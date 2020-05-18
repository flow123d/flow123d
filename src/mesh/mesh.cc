/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    mesh.cc
 * @ingroup mesh
 * @brief   Mesh construction
 */

#include <unistd.h>
#include <set>
#include <unordered_map>

#include "system/system.hh"
#include "system/exceptions.hh"
#include "system/index_types.hh"
#include "input/reader_to_storage.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"
#include "la/distribution.hh"

#include "mesh/mesh.h"
#include "mesh/bc_mesh.hh"
#include "mesh/ref_element.hh"
#include "mesh/region_set.hh"
#include "mesh/range_wrapper.hh"

// think about following dependencies
#include "mesh/accessors.hh"
#include "mesh/node_accessor.hh"
#include "mesh/partitioning.hh"
#include "mesh/neighbours.h"


#include "mesh/bih_tree.hh"
#include "mesh/duplicate_nodes.h"

#include "intersection/mixed_mesh_intersections.hh"



//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "mesh/region.hh"

#define NDEF  -1

namespace IT = Input::Type;

const Input::Type::Selection & Mesh::get_input_intersection_variant() {
    return Input::Type::Selection("Types of search algorithm for finding intersection candidates.")
        .add_value(Mesh::BIHsearch, "BIHsearch",
            "Use BIH for finding initial candidates, then continue by prolongation.")
        .add_value(Mesh::BIHonly, "BIHonly",
            "Use BIH for finding all candidates.")
        .add_value(Mesh::BBsearch, "BBsearch",
            "Use bounding boxes for finding initial candidates, then continue by prolongation.")
        .close();
}

const IT::Record & Mesh::get_input_type() {
	return IT::Record("Mesh","Record with mesh related data." )
	    .allow_auto_conversion("mesh_file")
		.declare_key("mesh_file", IT::FileName::input(), IT::Default::obligatory(),
				"Input file with mesh description.")
		.declare_key("regions", IT::Array( RegionSetBase::get_input_type() ), IT::Default::optional(),
				"List of additional region and region set definitions not contained in the mesh. "
				"There are three region sets implicitly defined:\n\n"
				"- ALL (all regions of the mesh)\n"
				"- .BOUNDARY (all boundary regions)\n"
				"- BULK (all bulk regions)")
		.declare_key("partitioning", Partitioning::get_input_type(), IT::Default("\"any_neighboring\""), "Parameters of mesh partitioning algorithms.\n" )
	    .declare_key("print_regions", IT::Bool(), IT::Default("true"), "If true, print table of all used regions.")
        .declare_key("intersection_search", Mesh::get_input_intersection_variant(), 
                     IT::Default("\"BIHsearch\""), "Search algorithm for element intersections.")
        .declare_key("global_snap_radius", IT::Double(0.0), IT::Default("1E-3"),
                     "Maximal snapping distance from the mesh in various search operations. In particular, it is used "
                     "to find the closest mesh element of an observe point; and in FieldFormula to find closest surface "
                     "element in plan view (Z projection).")
        .declare_key("raw_ngh_output", IT::FileName::output(), IT::Default::optional(),
                     "Output file with neighboring data from mesh.")
		.close();
}

const unsigned int Mesh::undef_idx;

Mesh::Mesh()
: tree(nullptr),
  bulk_size_(0),
  nodes_(3, 1, 0),
  row_4_el(nullptr),
  el_4_loc(nullptr),
  el_ds(nullptr),
  node_4_loc_(nullptr),
  node_ds_(nullptr),  
  bc_mesh_(nullptr)
  
{}



Mesh::Mesh(Input::Record in_record, MPI_Comm com)
: tree(nullptr),
  in_record_(in_record),
  comm_(com),
  bulk_size_(0),
  nodes_(3, 1, 0),
  row_4_el(nullptr),
  el_4_loc(nullptr),
  el_ds(nullptr),
  node_4_loc_(nullptr),
  node_ds_(nullptr),
  bc_mesh_(nullptr)
{
	// set in_record_, if input accessor is empty
	if (in_record_.is_empty()) {
		istringstream is("{mesh_file=\"\"}");
	    Input::ReaderToStorage reader;
	    IT::Record &in_rec = const_cast<IT::Record &>(Mesh::get_input_type());
	    in_rec.finish();
	    reader.read_stream(is, in_rec, Input::FileFormat::format_JSON);
	    in_record_ = reader.get_root_interface<Input::Record>();
	}

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            // optionally open raw output file
            FilePath raw_output_file_path;
            if (in_record_.opt_val("raw_ngh_output", raw_output_file_path)) {
            	MessageOut() << "Opening raw ngh output: " << raw_output_file_path << "\n";
            	try {
            		raw_output_file_path.open_stream(raw_ngh_output_file);
            	} INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, (in_record_))
            }

        }
	init();
}

Mesh::IntersectionSearch Mesh::get_intersection_search()
{
    return in_record_.val<Mesh::IntersectionSearch>("intersection_search");
}


void Mesh::init()
{

    n_insides = NDEF;
    n_exsides = NDEF;
    n_sides_ = NDEF;

    // number of element of particular dimension
    n_lines = 0;
    n_triangles = 0;
    n_tetrahedras = 0;

    for (int d=0; d<3; d++) max_edge_sides_[d] = 0;

    // Initialize numbering of nodes on sides.
    // This is temporary solution, until class Element is templated
    // by dimension. Then we can replace Mesh::side_nodes by
    // RefElement<dim>::side_nodes.

    // indices of side nodes in element node array
    // Currently this is made ad libitum
    // with some ordering here we can get sides with correct orientation.
    // This speedup normal calculation.

    side_nodes.resize(3); // three side dimensions
    for(int dim=0; dim < 3; dim++) {
        side_nodes[dim].resize(dim+2); // number of sides
        for(int i_side=0; i_side < dim+2; i_side++)
            side_nodes[dim][i_side].resize(dim+1);
    }

    for (unsigned int sid=0; sid<RefElement<1>::n_sides; sid++)
    	for (unsigned int nid=0; nid<RefElement<1>::n_nodes_per_side; nid++)
            side_nodes[0][sid][nid] = RefElement<1>::interact(Interaction<0,0>(sid))[nid];

    for (unsigned int sid=0; sid<RefElement<2>::n_sides; sid++)
        	for (unsigned int nid=0; nid<RefElement<2>::n_nodes_per_side; nid++)
                side_nodes[1][sid][nid] = RefElement<2>::interact(Interaction<0,1>(sid))[nid];

    for (unsigned int sid=0; sid<RefElement<3>::n_sides; sid++)
        	for (unsigned int nid=0; nid<RefElement<3>::n_nodes_per_side; nid++)
        		side_nodes[2][sid][nid] = RefElement<3>::interact(Interaction<0,2>(sid))[nid];
}


Mesh::~Mesh() {
    for(EdgeData &edg : this->edges)
        if (edg.side_) delete[] edg.side_;

    for (unsigned int idx=0; idx < bulk_size_; idx++) {
    	Element *ele=&(element_vec_[idx]);
        if (ele->boundary_idx_) delete[] ele->boundary_idx_;
        if (ele->neigh_vb) delete[] ele->neigh_vb;
    }

    for(unsigned int idx=bulk_size_; idx < element_vec_.size(); idx++) {
        Element *ele=&(element_vec_[idx]);
        if (ele->boundary_idx_) delete[] ele->boundary_idx_;
    }

    if (row_4_el != nullptr) delete[] row_4_el;
    if (el_4_loc != nullptr) delete[] el_4_loc;
    if (el_ds != nullptr) delete el_ds;
    if (node_4_loc_ != nullptr) delete[] node_4_loc_;
    if (node_ds_ != nullptr) delete node_ds_;
    if (bc_mesh_ != nullptr) delete bc_mesh_;
    if (tree != nullptr) delete tree;
}


unsigned int Mesh::n_sides() const
{
    if (n_sides_ == NDEF) {
        n_sides_=0;
        for (auto ele : this->elements_range()) n_sides_ += ele->n_sides();
    }
    return n_sides_;
}

unsigned int Mesh::n_vb_neighbours() const {
     return vb_neighbours_.size();
 }


unsigned int Mesh::n_corners() {
    unsigned int li, count = 0;
    for (auto ele : this->elements_range()) {
    	for (li=0; li<ele->n_nodes(); li++) {
            count++;
        }
    }
    return count;
}

Edge Mesh::edge(uint edge_idx) const
{
    ASSERT_LT_DBG(edge_idx, edges.size());
    return Edge(this, edge_idx);
}

Boundary Mesh::boundary(uint bc_idx) const
{
    ASSERT_LT_DBG(bc_idx, boundary_.size());
    return Boundary(&boundary_[bc_idx]);
}

Partitioning *Mesh::get_part() {
    return part_.get();
}

const LongIdx *Mesh::get_local_part() {
    return (LongIdx*)this->get_part()->get_loc_part();
}


//=============================================================================
// COUNT ELEMENT TYPES
//=============================================================================

void Mesh::count_element_types() {
	for (auto elm : this->elements_range())
    switch (elm->dim()) {
        case 1:
            n_lines++;
            break;
        case 2:
            n_triangles++;
            break;
        case 3:
            n_tetrahedras++;
            break;
        }
}



void Mesh::modify_element_ids(const RegionDB::MapElementIDToRegionID &map) {

    // get dim of the first element in the map, if it exists
    uint dim_to_check = RegionDB::undefined_dim;
    std::string reg_name = "UndefinedRegion";
    if(map.size() > 0){
        Element &ele = element_vec_[ elem_index(map.begin()->first) ];
        dim_to_check = ele.dim();
        reg_name = region_db_.find_id(map.begin()->second).label();
    }

	for (auto elem_to_region : map) {
		Element &ele = element_vec_[ elem_index(elem_to_region.first) ];
        
        if( ele.dim() != dim_to_check){
            xprintf(UsrErr, "User defined region '%s' (id %d) by 'From_Elements' cannot have elements of different dimensions.\n"
                            "Thrown due to: dim %d neq dim %d (ele id %d).\n"
                            "Split elements by dim, create separate regions and then possibly use Union.\n",
                    reg_name.c_str(), elem_to_region.second, dim_to_check, ele.dim(), elem_to_region.first);
        }

		ele.region_idx_ = region_db_.get_region( elem_to_region.second, ele.dim() );
		region_db_.mark_used_region(ele.region_idx_.idx());
	}
}


void Mesh::check_mesh_on_read() {
    std::vector<uint> nodes_new_idx( this->n_nodes(), Mesh::undef_idx );

    // check element quality and flag used nodes
    for (auto ele : this->elements_range()) {
        // element quality
    	double quality = ele.quality_measure_smooth(ele.side(0));
        if ( quality< 0.001)
            WarningOut().fmt("Bad quality (<0.001) of the element {}.\n", ele.idx());

        // flag used nodes
        for (uint ele_node=0; ele_node<ele->n_nodes(); ele_node++) {
            uint inode = ele->node_idx(ele_node);
            nodes_new_idx[inode] = inode;
        }
    }

    // possibly build new node ids map
    BidirectionalMap<int> new_node_ids_;
    new_node_ids_.reserve(node_ids_.size());

    // remove unused nodes from the mesh
    uint inode_new = 0;
    for(uint inode = 0; inode < nodes_new_idx.size(); inode++) {
        if(nodes_new_idx[inode] == Mesh::undef_idx){
            WarningOut().fmt("A node {} does not belong to any element "
                         " and will be removed.",
                         find_node_id(inode));
        }
        else{
            // map new node numbering
            nodes_new_idx[inode] = inode_new;
            
            // possibly move the nodes
            nodes_.vec<3>(inode_new) = nodes_.vec<3>(inode);
            new_node_ids_.add_item(node_ids_[inode]);

            inode_new++;
        }
    }

    uint n_nodes_new = inode_new;

    // if some node erased, update node ids in elements
    if(n_nodes_new < nodes_new_idx.size()){
        
        DebugOut() << "Updating node-element numbering due to unused nodes: "
            << print_var(n_nodes_new) << print_var(nodes_new_idx.size()) << "\n";

        // throw away unused nodes
        nodes_.resize(n_nodes_new);
        node_ids_ = new_node_ids_;

        // update node-element numbering
        for (auto ele : this->elements_range()) {
            for (uint ele_node=0; ele_node<ele->n_nodes(); ele_node++) {
                uint inode_orig = ele->node_idx(ele_node);
                uint inode = nodes_new_idx[inode_orig];
                ASSERT_DBG(inode != Mesh::undef_idx);
                const_cast<Element*>(ele.element())->nodes_[ele_node] = inode;
            }
        }
    }
}

void Mesh::setup_topology() {
    START_TIMER("MESH - setup topology");
    
    count_element_types();
    check_mesh_on_read();

    make_neighbours_and_edges();
    element_to_neigh_vb();
    make_edge_permutations();
    count_side_types();
    
    tree = new DuplicateNodes(this);

    part_ = std::make_shared<Partitioning>(this, in_record_.val<Input::Record>("partitioning") );

    // create parallel distribution and numbering of elements
    LongIdx *id_4_old = new LongIdx[n_elements()];
    int i = 0;
    for (auto ele : this->elements_range())
        id_4_old[i++] = ele.idx();
    part_->id_maps(n_elements(), id_4_old, el_ds, el_4_loc, row_4_el);

    delete[] id_4_old;
    
    this->distribute_nodes();

    output_internal_ngh_data();
}


//
void Mesh::count_side_types()
{

    n_insides = 0;
    n_exsides = 0;
	for (auto ele : this->elements_range())
        for(SideIter sde = ele.side(0); sde->side_idx() < ele->n_sides(); ++sde) {
            if (sde->is_external()) n_exsides++;
            else n_insides++;
        }
}



void Mesh::create_node_element_lists() {
    // for each node we make a list of elements that use this node
    node_elements_.resize( this->n_nodes() );

    for (auto ele : this->elements_range())
        for (unsigned int n=0; n<ele->n_nodes(); n++)
            node_elements_[ele.node(n).idx()].push_back(ele.idx());

    for (vector<vector<unsigned int> >::iterator n=node_elements_.begin(); n!=node_elements_.end(); n++)
        stable_sort(n->begin(), n->end());
}


void Mesh::intersect_element_lists(vector<unsigned int> const &nodes_list, vector<unsigned int> &intersection_element_list)
{
	if (node_elements_.size() == 0) {
		this->create_node_element_lists();
	}

	if (nodes_list.size() == 0) {
        intersection_element_list.clear();
    } else if (nodes_list.size() == 1) {
        intersection_element_list = node_elements_[ nodes_list[0] ];
	} else {
	    vector<unsigned int>::const_iterator it1=nodes_list.begin();
	    vector<unsigned int>::const_iterator it2=it1+1;
	    intersection_element_list.resize( node_elements_[*it1].size() ); // make enough space

	    it1=set_intersection(
                node_elements_[*it1].begin(), node_elements_[*it1].end(),
                node_elements_[*it2].begin(), node_elements_[*it2].end(),
                intersection_element_list.begin());
        intersection_element_list.resize(it1-intersection_element_list.begin()); // resize to true size

        for(;it2<nodes_list.end();++it2) {
            it1=set_intersection(
                    intersection_element_list.begin(), intersection_element_list.end(),
                    node_elements_[*it2].begin(), node_elements_[*it2].end(),
                    intersection_element_list.begin());
            intersection_element_list.resize(it1-intersection_element_list.begin()); // resize to true size
        }
    }
}


bool Mesh::find_lower_dim_element( vector<unsigned int> &element_list, unsigned int dim, unsigned int &element_idx) {
    bool is_neighbour = false;

    vector<unsigned int>::iterator e_dest=element_list.begin();
    for( vector<unsigned int>::iterator ele = element_list.begin(); ele!=element_list.end(); ++ele) {
        if (element_vec_[*ele].dim() == dim) { // keep only indexes of elements of same dimension
            *e_dest=*ele;
            ++e_dest;
        } else if (element_vec_[*ele].dim() == dim-1) { // get only first element of lower dimension
            if (is_neighbour) xprintf(UsrErr, "Too matching elements id: %d and id: %d in the same mesh.\n",
            		this->elem_index(*ele), this->elem_index(element_idx) );

            is_neighbour = true;
            element_idx = *ele;
        }
    }
    element_list.resize( e_dest - element_list.begin());
    return is_neighbour;
}

bool Mesh::same_sides(const SideIter &si, vector<unsigned int> &side_nodes) {
    // check if nodes lists match (this is slow and will be faster only when we convert whole mesh into hierarchical design like in deal.ii)
    unsigned int ni=0;
    while ( ni < si->n_nodes()
        && find(side_nodes.begin(), side_nodes.end(), si->node(ni).idx() ) != side_nodes.end() ) ni++;
    return ( ni == si->n_nodes() );
}

/**
 * TODO:
 * - use std::is_any for setting is_neigbour
 * - possibly make appropriate constructors for Edge and Neighbour
 * - check side!=-1 when searching neigbouring element
 * - process boundary elements first, there should be no Neigh, but check it
 *   set Edge and boundary there
 */

void Mesh::make_neighbours_and_edges()
{
	ASSERT(bc_element_tmp_.size()==0)
			.error("Temporary structure of boundary element data is not empty. Did you call create_boundary_elements?");

    Neighbour neighbour;
    EdgeData *edg = nullptr;
    unsigned int ngh_element_idx;
    unsigned int last_edge_idx = Mesh::undef_idx;

    neighbour.mesh_ = this;

    create_node_element_lists();

	// pointers to created edges
	//vector<Edge *> tmp_edges;
    edges.resize(0); // be sure that edges are empty

	vector<unsigned int> side_nodes;
	vector<unsigned int> intersection_list; // list of elements in intersection of node element lists

	for( unsigned int i=bulk_size_; i<element_vec_.size(); ++i) {

		ElementAccessor<3> bc_ele = this->element_accessor(i);
		ASSERT(bc_ele.region().is_boundary());
        // Find all elements that share this side.
        side_nodes.resize(bc_ele->n_nodes());
        for (unsigned n=0; n<bc_ele->n_nodes(); n++) side_nodes[n] = bc_ele->node_idx(n);
        intersect_element_lists(side_nodes, intersection_list);
        bool is_neighbour = find_lower_dim_element(intersection_list, bc_ele->dim() +1, ngh_element_idx);
        if (is_neighbour) {
            xprintf(UsrErr, "Boundary element (id: %d) match a regular element (id: %d) of lower dimension.\n",
                    bc_ele.idx(), this->elem_index(ngh_element_idx));
        } else {
            if (intersection_list.size() == 0) {
                // no matching dim+1 element found
            	WarningOut().fmt("Lonely boundary element, id: {}, region: {}, dimension {}.\n",
            			bc_ele.idx(), bc_ele.region().id(), bc_ele->dim());
                continue; // skip the boundary element
            }
            last_edge_idx=edges.size();
            edges.resize(last_edge_idx+1);
            edg = &( edges.back() );
            edg->n_sides = 0;
            edg->side_ = new struct SideIter[ intersection_list.size() ];

            // common boundary object
            unsigned int bdr_idx=boundary_.size();
            boundary_.resize(bdr_idx+1);
            BoundaryData &bdr=boundary_.back();
            bdr.bc_ele_idx_ = i;
            bdr.edge_idx_ = last_edge_idx;
            bdr.mesh_=this;

            // for 1d boundaries there can be more then one 1d elements connected to the boundary element
            // we do not detect this case later in the main search over bulk elements
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect)  {
                ElementAccessor<3> elem = this->element_accessor(*isect);
                for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++) {
                    SideIter si = elem.side(ecs);
                    if ( same_sides( si, side_nodes) ) {
                        if (elem->edge_idx(ecs) != Mesh::undef_idx) {
                        	OLD_ASSERT(elem->boundary_idx_!=nullptr, "Null boundary idx array.\n");
                            int last_bc_ele_idx=this->boundary_[elem->boundary_idx_[ecs]].bc_ele_idx_;
                            int new_bc_ele_idx=i;
                            THROW( ExcDuplicateBoundary()
                                    << EI_ElemLast(this->find_elem_id(last_bc_ele_idx))
                                    << EI_RegLast(this->element_accessor(last_bc_ele_idx).region().label())
                                    << EI_ElemNew(this->find_elem_id(new_bc_ele_idx))
                                    << EI_RegNew(this->element_accessor(new_bc_ele_idx).region().label())
                                    );
                        }
                        element_vec_[*isect].edge_idx_[ecs] = last_edge_idx;
                        edg->side_[ edg->n_sides++ ] = si;

                        if (elem->boundary_idx_ == NULL) {
                        	Element *el = &(element_vec_[*isect]);
                        	el->boundary_idx_ = new unsigned int [ el->n_sides() ];
                            std::fill( el->boundary_idx_, el->boundary_idx_ + el->n_sides(), Mesh::undef_idx);
                        }
                        elem->boundary_idx_[ecs] = bdr_idx;
                        break; // next element in intersection list
                    }
                }
            }

        }

	}
	// Now we go through all element sides and create edges and neighbours
	for (auto e : this->elements_range()) {
		for (unsigned int s=0; s<e->n_sides(); s++)
		{
			// skip sides that were already found
			if (e->edge_idx(s) != Mesh::undef_idx) continue;


			// Find all elements that share this side.
			side_nodes.resize(e.side(s)->n_nodes());
			for (unsigned n=0; n<e.side(s)->n_nodes(); n++) side_nodes[n] = e.side(s)->node(n).idx();
			intersect_element_lists(side_nodes, intersection_list);

			bool is_neighbour = find_lower_dim_element(intersection_list, e->dim(), ngh_element_idx);

			if (is_neighbour) { // edge connects elements of different dimensions
				// Initialize for the neighbour case.
			    neighbour.elem_idx_ = ngh_element_idx;
            } else { // edge connects only elements of the same dimension
                // Initialize for the edge case.
                last_edge_idx=edges.size();
                edges.resize(last_edge_idx+1);
                edg = &( edges.back() );
                edg->n_sides = 0;
                edg->side_ = new struct SideIter[ intersection_list.size() ];
                if (intersection_list.size() > max_edge_sides_[e->dim()-1])
                	max_edge_sides_[e->dim()-1] = intersection_list.size();

                if (intersection_list.size() == 1) {
                	// outer edge, create boundary object as well
                	Element &elm = element_vec_[e.idx()];
                    edg->n_sides=1;
                    edg->side_[0] = e.side(s);
                    element_vec_[e.idx()].edge_idx_[s] = last_edge_idx;

                    if (e->boundary_idx_ == NULL) {
                    	elm.boundary_idx_ = new unsigned int [ e->n_sides() ];
                        std::fill( elm.boundary_idx_, elm.boundary_idx_ + e->n_sides(), Mesh::undef_idx);
                    }

                    unsigned int bdr_idx=boundary_.size()+1; // need for VTK mesh that has no boundary elements
                                                             // and bulk elements are indexed from 0
                    boundary_.resize(bdr_idx+1);
                    BoundaryData &bdr=boundary_.back();
                    elm.boundary_idx_[s] = bdr_idx;

                    // fill boundary element
                    Element * bc_ele = add_element_to_vector(-bdr_idx);
                    bc_ele->init(e->dim()-1, region_db_.implicit_boundary_region() );
                    region_db_.mark_used_region( bc_ele->region_idx_.idx() );
                    for(unsigned int ni = 0; ni< side_nodes.size(); ni++) bc_ele->nodes_[ni] = side_nodes[ni];

                    // fill Boundary object
                    bdr.edge_idx_ = last_edge_idx;
                    bdr.bc_ele_idx_ = elem_index(-bdr_idx);
                    bdr.mesh_=this;

                    continue; // next side of element e
                }
			}

			// go through the elements connected to the edge or neighbour
			// setup neigbour or edge
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect) {
            	ElementAccessor<3> elem = this->element_accessor(*isect);
                for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++) {
                    if (elem->edge_idx(ecs) != Mesh::undef_idx) continue; // ??? This should not happen.
                    SideIter si = elem.side(ecs);
                    if ( same_sides( si, side_nodes) ) {
                        if (is_neighbour) {
                            // create a new edge and neighbour for this side, and element to the edge
                            last_edge_idx=edges.size();
                            edges.resize(last_edge_idx+1);
                            edg = &( edges.back() );
                            edg->n_sides = 1;
                            edg->side_ = new struct SideIter[1];
                            edg->side_[0] = si;
                            element_vec_[elem.idx()].edge_idx_[ecs] = last_edge_idx;

                            neighbour.edge_idx_ = last_edge_idx;

                            vb_neighbours_.push_back(neighbour); // copy neighbour with this edge setting
                        } else {
                            // connect the side to the edge, and side to the edge
                            ASSERT_PTR_DBG(edg);
                            edg->side_[ edg->n_sides++ ] = si;
                            ASSERT_DBG(last_edge_idx != Mesh::undef_idx);
                            element_vec_[elem.idx()].edge_idx_[ecs] = last_edge_idx;
                        }
                        break; // next element from intersection list
                    }
                } // search for side of other connected element
            } // connected elements

            if (! is_neighbour)
				ASSERT_EQ( (unsigned int) edg->n_sides, intersection_list.size())(e.index())(s).error("Missing edge sides.");
		} // for element sides
	}   // for elements

	MessageOut().fmt( "Created {} edges and {} neighbours.\n", edges.size(), vb_neighbours_.size() );
}



void Mesh::make_edge_permutations()
{
    // node numbers is the local index of the node on the last side
    // this maps the side nodes to the nodes of the reference side(0)
    std::unordered_map<unsigned int,unsigned int> node_numbers;

    for (auto edg : edge_range())
	{
        unsigned int n_side_nodes = edg.side(0)->n_nodes();
		// side 0 is reference, so its permutation is 0
		edg.side(0)->element()->permutation_idx_[edg.side(0)->side_idx()] = 0;

		if (edg.n_sides() > 1)
		{
		    // For every node on the reference side(0) give its local idx on the current side.
		    unsigned int permutation[n_side_nodes];


		    node_numbers.clear();
			for (uint i=0; i<n_side_nodes; i++)
				node_numbers[edg.side(0)->node(i).idx()] = i;

			for (uint sid=1; sid<edg.n_sides(); sid++)
			{
				for (uint i=0; i<n_side_nodes; i++)
					permutation[node_numbers[edg.side(sid)->node(i).idx()]] = i;

				switch (edg.side(0)->dim())
				{
				case 0:
					edg.side(sid)->element()->permutation_idx_[edg.side(sid)->side_idx()] = RefElement<1>::permutation_index(permutation);
					break;
				case 1:
					edg.side(sid)->element()->permutation_idx_[edg.side(sid)->side_idx()] = RefElement<2>::permutation_index(permutation);
					break;
				case 2:
					edg.side(sid)->element()->permutation_idx_[edg.side(sid)->side_idx()] = RefElement<3>::permutation_index(permutation);
					break;
				}
			}
		}
	}

	for (vector<Neighbour>::iterator nb=vb_neighbours_.begin(); nb!=vb_neighbours_.end(); nb++)
	{
        // node numbers is the local index of the node on the last side
        // this maps the side nodes to the nodes of the side(0)
	    unsigned int n_side_nodes = nb->element()->n_nodes();
	    unsigned int permutation[n_side_nodes];
        node_numbers.clear();

		// element of lower dimension is reference, so
		// we calculate permutation for the adjacent side
		for (unsigned int i=0; i<n_side_nodes; i++)
			node_numbers[nb->element().node(i).idx()] = i;

		for (unsigned int i=0; i<n_side_nodes; i++)
			permutation[node_numbers[nb->side()->node(i).idx()]] = i;

		switch (nb->side()->dim())
		{
		case 0:
			nb->side()->element()->permutation_idx_[nb->side()->side_idx()] = RefElement<1>::permutation_index(permutation);
			break;
		case 1:
			nb->side()->element()->permutation_idx_[nb->side()->side_idx()] = RefElement<2>::permutation_index(permutation);
			break;
		case 2:
			nb->side()->element()->permutation_idx_[nb->side()->side_idx()] = RefElement<3>::permutation_index(permutation);
			break;
		}
	}
}





//=============================================================================
//
//=============================================================================
void Mesh::element_to_neigh_vb()
{

	//MessageOut() << "Element to neighbours of vb2 type... "/*orig verb 5*/;

	for (vector<Element>::iterator ele = element_vec_.begin(); ele!= element_vec_.begin()+bulk_size_; ++ele)
		ele->n_neighs_vb_ =0;

    // count vb neighs per element
    for (auto & ngh : this->vb_neighbours_)  ngh.element()->n_neighs_vb_++;

    // Allocation of the array per element
    for (vector<Element>::iterator ele = element_vec_.begin(); ele!= element_vec_.begin()+bulk_size_; ++ele)
        if( ele->n_neighs_vb() > 0 ) {
            ele->neigh_vb = new struct Neighbour* [ele->n_neighs_vb()];
            ele->n_neighs_vb_=0;
        }

    // fill
    ElementAccessor<3> ele;
    for (auto & ngh : this->vb_neighbours_) {
        ele = ngh.element();
        ele->neigh_vb[ ele->n_neighs_vb_++ ] = &ngh;
    }

    //MessageOut() << "... O.K.\n"/*orig verb 6*/;
}






MixedMeshIntersections & Mesh::mixed_intersections() {
	/* Algorithm:
	 *
	 * 1) create BIH tree
	 * 2) for every 1D, find list of candidates
	 * 3) compute intersections for 1d, store it to master_elements
	 *
	 */
    if (! intersections) {
        intersections = std::make_shared<MixedMeshIntersections>(this);
        intersections->compute_intersections();
    }
    return *intersections;
}



ElementAccessor<3> Mesh::element_accessor(unsigned int idx) const {
    return ElementAccessor<3>(this, idx);
}



NodeAccessor<3> Mesh::node(unsigned int idx) const {
    return NodeAccessor<3>(this, idx);
}



void Mesh::elements_id_maps( vector<LongIdx> & bulk_elements_id, vector<LongIdx> & boundary_elements_id) const
{
    if (bulk_elements_id.size() ==0) {
        std::vector<LongIdx>::iterator map_it;
        LongIdx last_id;

        bulk_elements_id.resize(n_elements());
        map_it = bulk_elements_id.begin();
        last_id = -1;
        for(unsigned int idx=0; idx < n_elements(); idx++, ++map_it) {
        	LongIdx id = this->find_elem_id(idx);
            if (last_id >= id) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", id);
            last_id=*map_it = id;
        }

        boundary_elements_id.resize(n_elements(true));
        map_it = boundary_elements_id.begin();
        last_id = -1;
        for(unsigned int idx=bulk_size_; idx<element_vec_.size(); idx++, ++map_it) {
        	LongIdx id = this->find_elem_id(idx);
            // We set ID for boundary elements created by the mesh itself to "-1"
            // this force gmsh reader to skip all remaining entries in boundary_elements_id
            // and thus report error for any remaining data lines
            if (id < 0) last_id=*map_it=-1;
            else {
                if (last_id >= id) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", id);
                last_id=*map_it = id;
            }
        }
    }
}


bool compare_points(const arma::vec3 &p1, const arma::vec3 &p2) {
    static const double point_tolerance = 1E-10;
	return fabs(p1[0]-p2[0]) < point_tolerance
		&& fabs(p1[1]-p2[1]) < point_tolerance
		&& fabs(p1[2]-p2[2]) < point_tolerance;
}


bool Mesh::check_compatible_mesh( Mesh & mesh, vector<LongIdx> & bulk_elements_id, vector<LongIdx> & boundary_elements_id )
{
	std::vector<unsigned int> node_ids; // allow mapping ids of nodes from source mesh to target mesh
	std::vector<unsigned int> node_list;
	std::vector<unsigned int> candidate_list; // returned by intersect_element_lists
	std::vector<unsigned int> result_list; // list of elements with same dimension as vtk element
	unsigned int i; // counter over vectors

    {
        // iterates over node vector of \p this object
        // to each node must be found just only one node in target \p mesh
        // store orders (mapping between source and target meshes) into node_ids vector
        std::vector<unsigned int> searched_elements; // for BIH tree
        unsigned int i_node, i_elm_node;
        const BIHTree &bih_tree=mesh.get_bih_tree();

    	// create nodes of mesh
        node_ids.resize( this->n_nodes() );
        i=0;
        for (auto nod : this->node_range()) {
            uint found_i_node = Mesh::undef_idx;
            bih_tree.find_point(*nod, searched_elements);

            for (std::vector<unsigned int>::iterator it = searched_elements.begin(); it!=searched_elements.end(); it++) {
                ElementAccessor<3> ele = mesh.element_accessor( *it );
                for (i_node=0; i_node<ele->n_nodes(); i_node++)
                {
                    static const double point_tolerance = 1E-10;
                    if ( arma::norm(*ele.node(i_node) - *nod, 1) < point_tolerance) {
                        i_elm_node = ele.node(i_node).idx();
                        if (found_i_node == Mesh::undef_idx) found_i_node = i_elm_node;
                        else if (found_i_node != i_elm_node) {
                            // duplicate nodes in target mesh
                        	this->elements_id_maps(bulk_elements_id, boundary_elements_id);
                            return false;
                        }
                    }
                }
            }
            if (found_i_node == Mesh::undef_idx) {
                // no node found in target mesh
            	this->elements_id_maps(bulk_elements_id, boundary_elements_id);
            	return false;
            }
            node_ids[i] = found_i_node;
            searched_elements.clear();
            i++;
        }
    }

    {
        // iterates over bulk elements of \p this object
        // elements in both meshes must be in ratio 1:1
        // store orders (mapping between both mesh files) into bulk_elements_id vector
        bulk_elements_id.clear();
        bulk_elements_id.resize(this->n_elements());
        // iterate trough bulk part of element vector, to each element in source mesh must exist only one element in target mesh
        // fill bulk_elements_id vector
        i=0;
        for (auto elm : this->elements_range()) {
            for (unsigned int j=0; j<elm->n_nodes(); j++) { // iterate trough all nodes of any element
                node_list.push_back( node_ids[ elm->node_idx(j) ] );
            }
            mesh.intersect_element_lists(node_list, candidate_list);
            for (auto i_elm : candidate_list) {
            	if ( mesh.element_accessor(i_elm)->dim() == elm.dim() ) result_list.push_back( elm.index() );
            }
            if (result_list.size() != 1) {
            	// intersect_element_lists must produce one element
            	this->elements_id_maps(bulk_elements_id, boundary_elements_id);
            	return false;
            }
            bulk_elements_id[i] = (LongIdx)result_list[0];
            node_list.clear();
            result_list.clear();
        	i++;
        }
    }

    {
        // iterates over boundary elements of \p this object
        // elements in both meshes must be in ratio 1:1
        // store orders (mapping between both mesh files) into boundary_elements_id vector
    	auto bc_mesh = this->get_bc_mesh();
        boundary_elements_id.clear();
        boundary_elements_id.resize(bc_mesh->n_elements());
        // iterate trough boundary part of element vector, to each element in source mesh must exist only one element in target mesh
        // fill boundary_elements_id vector
        i=0;
        for (auto elm : bc_mesh->elements_range()) {
            for (unsigned int j=0; j<elm->n_nodes(); j++) { // iterate trough all nodes of any element
                node_list.push_back( node_ids[ elm->node_idx(j) ] );
            }
            mesh.get_bc_mesh()->intersect_element_lists(node_list, candidate_list);
            for (auto i_elm : candidate_list) {
            	if ( mesh.get_bc_mesh()->element_accessor(i_elm)->dim() == elm.dim() ) result_list.push_back( elm.index() );
            }
            if (result_list.size() != 1) {
            	// intersect_element_lists must produce one element
            	this->elements_id_maps(bulk_elements_id, boundary_elements_id);
            	return false;
            }
            boundary_elements_id[i] = (LongIdx)result_list[0];
            node_list.clear();
            result_list.clear();
        	i++;
        }
    }

    return true;
}

void Mesh::read_regions_from_input(Input::Array region_list)
{
	for (Input::Iterator<Input::AbstractRecord> it = region_list.begin<Input::AbstractRecord>();
				it != region_list.end();
				++it) {
		// constructor has side effect in the mesh - create new region or set and store them to Mesh::region_db_
		(*it).factory< RegionSetBase, const Input::Record &, Mesh * >(*it, this);
	}
}

void Mesh::check_and_finish()
{
	modify_element_ids(region_db_.el_to_reg_map_);
	region_db_.el_to_reg_map_.clear();
	region_db_.close();
	region_db_.check_regions();

	if ( in_record_.val<bool>("print_regions") ) {
		stringstream ss;
		region_db_.print_region_table(ss);
		MessageOut() << ss.str();
	}
}


std::vector<BoundingBox> Mesh::get_element_boxes() {
    START_TIMER("Mesh::compute_element_boxes");
    std::vector<BoundingBox> boxes;

    // make element boxes
    unsigned int i=0;
    boxes.resize(this->n_elements());
    for (auto element : this->elements_range()) {
        boxes[i] = element.bounding_box();
    	i++;
    }

    return boxes;
}

const BIHTree &Mesh::get_bih_tree() {
    if (! this->bih_tree_) {
        bih_tree_ = std::make_shared<BIHTree>();
        bih_tree_->add_boxes( this->get_element_boxes() );
        bih_tree_->construct();
	}
    return *bih_tree_;
}

double Mesh::global_snap_radius() const {
	return in_record_.val<double>("global_snap_radius");
}

void Mesh::add_physical_name(unsigned int dim, unsigned int id, std::string name) {
	region_db_.add_region(id, name, dim, "$PhysicalNames");
}


void Mesh::add_node(unsigned int node_id, arma::vec3 coords) {

    nodes_.append(coords);
    node_ids_.add_item(node_id);
}


void Mesh::add_element(unsigned int elm_id, unsigned int dim, unsigned int region_id, unsigned int partition_id,
		std::vector<unsigned int> node_ids) {
	RegionIdx region_idx = region_db_.get_region( region_id, dim );
	if ( !region_idx.is_valid() ) {
		region_idx = region_db_.add_region( region_id, region_db_.create_label_from_id(region_id), dim, "$Element" );
	}
	region_db_.mark_used_region(region_idx.idx());

	if (region_idx.is_boundary()) {
		bc_element_tmp_.push_back( ElementTmpData(elm_id, dim, region_idx, partition_id, node_ids) );
	} else {
		if(dim == 0 ) {
			WarningOut().fmt("Bulk elements of zero size(dim=0) are not supported. Element ID: {}.\n", elm_id);
		}
		else {
			Element *ele = add_element_to_vector(elm_id);
			bulk_size_++;
			this->init_element(ele, elm_id, dim, region_idx, partition_id, node_ids);
		}
	}
}


void Mesh::init_element(Element *ele, unsigned int elm_id, unsigned int dim, RegionIdx region_idx, unsigned int partition_id,
		std::vector<unsigned int> node_ids) {
	ele->init(dim, region_idx);
	ele->pid_ = partition_id;

	for (unsigned int ni=0; ni<ele->n_nodes(); ni++) {
		ele->nodes_[ni] = this->node_index(node_ids[ni]);
	}

    // check that tetrahedron element is numbered correctly and is not degenerated
    if(ele->dim() == 3)
    {
        ElementAccessor<3> ea = this->element_accessor( this->elem_index(elm_id) );
        double jac = ea.tetrahedron_jacobian();
        if( ! (jac > 0) ) {
            WarningOut().fmt("Tetrahedron element with id {} has wrong numbering or is degenerated (Jacobian = {}).",elm_id, jac);
        }
    }
}


vector<vector<unsigned int> > const & Mesh::node_elements() {
	if (node_elements_.size() == 0) {
		this->create_node_element_lists();
	}
	return node_elements_;
}


void Mesh::init_element_vector(unsigned int size) {
	element_vec_.clear();
    element_ids_.clear();
	element_vec_.reserve(size);
    element_ids_.reserve(size);
	bc_element_tmp_.clear();
	bc_element_tmp_.reserve(size);
	bulk_size_ = 0;
	boundary_loaded_size_ = 0;
}


void Mesh::init_node_vector(unsigned int size) {
	nodes_.reinit(size);
	node_ids_.clear();
	node_ids_.reserve(size);
}


Element * Mesh::add_element_to_vector(int id) {
    element_vec_.push_back( Element() );
    Element * elem = &element_vec_.back(); //[element_vec_.size()-1];
    element_ids_.add_item((unsigned int)(id));
	return elem;
}

Range<ElementAccessor<3>> Mesh::elements_range() const {
	auto bgn_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(this, 0) );
	auto end_it = make_iter<ElementAccessor<3>>( ElementAccessor<3>(this, bulk_size_) );
	return Range<ElementAccessor<3>>(bgn_it, end_it);
}

Range<NodeAccessor<3>> Mesh::node_range() const {
	auto bgn_it = make_iter<NodeAccessor<3>>( NodeAccessor<3>(this, 0) );
	auto end_it = make_iter<NodeAccessor<3>>( NodeAccessor<3>(this, n_nodes()) );
    return Range<NodeAccessor<3>>(bgn_it, end_it);
}

Range<Edge> Mesh::edge_range() const {
	auto bgn_it = make_iter<Edge>( Edge(this, 0) );
	auto end_it = make_iter<Edge>( Edge(this, edges.size()) );
    return Range<Edge>(bgn_it, end_it);
}

inline void Mesh::check_element_size(unsigned int elem_idx) const
{
    ASSERT(elem_idx < element_vec_.size())(elem_idx)(element_vec_.size()).error("Index of element is out of bound of element vector!");
}

/*
 * Output of internal flow data.
 */
void Mesh::output_internal_ngh_data()
{
    START_TIMER("Mesh::output_internal_ngh_data");

    if (! raw_ngh_output_file.is_open()) return;
    
    // header
    raw_ngh_output_file <<  "// fields:\n//ele_id    n_sides    ns_side_neighbors[n]    neighbors[n*ns]    n_vb_neighbors    vb_neighbors[n_vb]\n";
    raw_ngh_output_file <<  fmt::format("{}\n" , n_elements());

    int cit = 0;
    
    // map from higher dim elements to its lower dim neighbors, using gmsh IDs: ele->id()
    unsigned int undefined_ele_id = -1;
    std::map<unsigned int, std::vector<unsigned int>> neigh_vb_map;
    for (auto ele : this->elements_range()) {
        if(ele->n_neighs_vb() > 0){
            for (unsigned int i = 0; i < ele->n_neighs_vb(); i++){
                ElementAccessor<3> higher_ele = ele->neigh_vb[i]->side()->element();
                
                auto search = neigh_vb_map.find(higher_ele.idx());
                if(search != neigh_vb_map.end()){
                    // if found, add id to correct local side idx
                    search->second[ele->neigh_vb[i]->side()->side_idx()] = ele.idx();
                }
                else{
                    // if not found, create new vector, each side can have one vb neighbour
                    std::vector<unsigned int> higher_ele_side_ngh(higher_ele->n_sides(), undefined_ele_id);
                    higher_ele_side_ngh[ele->neigh_vb[i]->side()->side_idx()] = ele.idx();
                    neigh_vb_map[higher_ele.idx()] = higher_ele_side_ngh;
                }
            }
        }
    }
    
    for (auto ele : this->elements_range()) {
        raw_ngh_output_file << ele.idx() << " ";
        raw_ngh_output_file << ele->n_sides() << " ";
        
        auto search_neigh = neigh_vb_map.end();
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            unsigned int n_side_neighs = ele.side(i)->edge().n_sides()-1;  //n_sides - the current one
            // check vb neighbors (lower dimension)
            if(n_side_neighs == 0){
                //update search
                if(search_neigh == neigh_vb_map.end())
                    search_neigh = neigh_vb_map.find(ele.idx());
                
                if(search_neigh != neigh_vb_map.end())
                    if(search_neigh->second[i] != undefined_ele_id)
                        n_side_neighs = 1;
            }
            raw_ngh_output_file << n_side_neighs << " ";
        }
        
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            Edge edge = ele.side(i)->edge();
            if(edge.n_sides() > 1){
                for (uint j = 0; j < edge.n_sides(); j++) {
                    if(edge.side(j) != ele.side(i))
                        raw_ngh_output_file << edge.side(j)->element().idx() << " ";
                }
            }
            //check vb neighbour
            else if(search_neigh != neigh_vb_map.end()
                    && search_neigh->second[i] != undefined_ele_id){
                raw_ngh_output_file << search_neigh->second[i] << " ";
            }
        }
        
        // list higher dim neighbours
        raw_ngh_output_file << ele->n_neighs_vb() << " ";
        for (unsigned int i = 0; i < ele->n_neighs_vb(); i++)
            raw_ngh_output_file << ele->neigh_vb[i]->side()->element().idx() << " ";
        
        raw_ngh_output_file << endl;
        cit ++;
    }
    raw_ngh_output_file << "$EndFlowField\n" << endl;
}


void Mesh::create_boundary_elements() {
    // Copy boundary elements in temporary storage to the second part of the element vector
	for(ElementTmpData &e_data : bc_element_tmp_) {
	    Element *ele = add_element_to_vector(e_data.elm_id);
		this->init_element(ele, e_data.elm_id, e_data.dim, e_data.region_idx,
				e_data.partition_id, e_data.node_ids);
	}
	// release memory
	vector<ElementTmpData>().swap(bc_element_tmp_);
}


void Mesh::permute_tetrahedron(unsigned int elm_idx, std::vector<unsigned int> permutation_vec)
{
	ASSERT_LT_DBG(elm_idx, element_vec_.size());
    ASSERT_EQ_DBG(permutation_vec.size(), 4);

    std::array<unsigned int, 4> tmp_nodes;
    Element &elem = element_vec_[elm_idx];
    ASSERT_EQ_DBG(elem.dim(), 3);

    for(unsigned int i=0; i<elem.n_nodes(); i++)
    {
       	tmp_nodes[i] = elem.nodes_[permutation_vec[i]];
    }
    elem.nodes_ = tmp_nodes;
}


void Mesh::permute_triangle(unsigned int elm_idx, std::vector<unsigned int> permutation_vec)
{
	ASSERT_LT_DBG(elm_idx, element_vec_.size());
	ASSERT_EQ_DBG(permutation_vec.size(), 3);

    std::array<unsigned int, 4> tmp_nodes;
    Element &elem = element_vec_[elm_idx];
    ASSERT_EQ_DBG(elem.dim(), 2);

    for(unsigned int i=0; i<elem.n_nodes(); i++)
    {
       	tmp_nodes[i] = elem.nodes_[permutation_vec[i]];
    }
    elem.nodes_ = tmp_nodes;
}


BCMesh *Mesh::get_bc_mesh() {
	if (bc_mesh_ == nullptr) bc_mesh_ = new BCMesh(this);
	return bc_mesh_;
}


void Mesh::distribute_nodes() {
    ASSERT_PTR(el_4_loc).error("Array 'el_4_loc' is not initialized. Did you call Partitioning::id_maps?\n");

    unsigned int i_proc, i_node, i_ghost_node, elm_node;
    unsigned int my_proc = el_ds->myp();
    unsigned int n_proc = el_ds->np();

    // distribute nodes between processes, every node is assigned to minimal process of elements that own node
    // fill node_proc vector with same values on all processes
    std::vector<unsigned int> node_proc( this->n_nodes(), n_proc );
    std::vector<bool> local_node_flag( this->n_nodes(), false );

    for ( auto elm : this->elements_range() ) {
        i_proc = elm.proc();
        for (elm_node=0; elm_node<elm->n_nodes(); elm_node++) {
            i_node = elm->node_idx(elm_node);
            if (i_proc == my_proc) local_node_flag[i_node] = true;
            if (i_proc < node_proc[i_node]) node_proc[i_node] = i_proc;
        }
    }

    unsigned int n_own_nodes=0, n_local_nodes=0; // number of own and ghost nodes
    for(uint loc_flag : local_node_flag) if (loc_flag) n_local_nodes++;
    for(uint i_proc : node_proc) {
        if (i_proc == my_proc)
            n_own_nodes++;
        else if (i_proc == n_proc)
            ASSERT(0)(find_node_id(n_own_nodes)).error("A node does not belong to any element!");
    }

    //DebugOut() << print_var(n_own_nodes) << print_var(n_local_nodes) << this->n_nodes();
    // create and fill node_4_loc_ (mapping local to global indexes)
    node_4_loc_ = new LongIdx [ n_local_nodes ];
    i_node=0;
    i_ghost_node = n_own_nodes;
    for (unsigned int i=0; i<this->n_nodes(); ++i) {
        if (local_node_flag[i]) {
            if (node_proc[i]==my_proc)
                node_4_loc_[i_node++] = i;
            else
                node_4_loc_[i_ghost_node++] = i;
        }
    }

    // Construct node distribution object, set number of local nodes (own+ghost)
    node_ds_ = new Distribution(n_own_nodes, PETSC_COMM_WORLD);
    node_ds_->get_lsizes_array(); // need to initialize lsizes data member
    n_local_nodes_ = n_local_nodes;

}

//-----------------------------------------------------------------------------
// vim: set cindent:
