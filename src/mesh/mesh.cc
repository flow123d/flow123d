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


#include "system/system.hh"
#include "system/exceptions.hh"
#include "input/reader_to_storage.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "system/sys_profiler.hh"
#include "la/distribution.hh"

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"
#include "mesh/region_set.hh"

// think about following dependencies
#include "mesh/boundaries.h"
#include "mesh/accessors.hh"
#include "mesh/partitioning.hh"


#include "mesh/bih_tree.hh"

#include "mesh/ngh/include/triangle.h"
#include "mesh/ngh/include/abscissa.h"
#include "mesh/ngh/include/intersection.h"

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
		.declare_key("global_observe_search_radius", IT::Double(0.0), IT::Default("1E-3"),
					 "Maximal distance of observe point from Mesh relative to its size (bounding box). "
					 "Value is global and it can be rewrite at arbitrary ObservePoint by setting the key search_radius.")
		.close();
}

const unsigned int Mesh::undef_idx;

Mesh::Mesh()
: row_4_el(nullptr),
  el_4_loc(nullptr),
  el_ds(nullptr)
{}



Mesh::Mesh(Input::Record in_record, MPI_Comm com)
: in_record_(in_record),
  comm_(com),
  row_4_el(nullptr),
  el_4_loc(nullptr),
  el_ds(nullptr)
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

	reinit(in_record_);
}

Mesh::IntersectionSearch Mesh::get_intersection_search()
{
    return in_record_.val<Mesh::IntersectionSearch>("intersection_search");
}


void Mesh::reinit(Input::Record in_record)
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
    for(int i=0; i < 3; i++) {
        side_nodes[i].resize(i+2); // number of sides
        for(int j=0; j < i+2; j++)
            side_nodes[i][j].resize(i+1);
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
    for(Edge &edg : this->edges)
        if (edg.side_) delete[] edg.side_;

    FOR_ELEMENTS( this, ele ) {
        if (ele->node) delete[] ele->node;
        if (ele->edge_idx_) delete[] ele->edge_idx_;
        if (ele->permutation_idx_) delete[] ele->permutation_idx_;
        if (ele->boundary_idx_) delete[] ele->boundary_idx_;
        if (ele->neigh_vb) delete[] ele->neigh_vb;
    }

    for(unsigned int idx=0; idx < element_vec_.size(); idx++) {
        Element *ele=&(element_vec_[idx]);
        if (ele->node) delete[] ele->node;
        if (ele->edge_idx_) delete[] ele->edge_idx_;
        if (ele->permutation_idx_) delete[] ele->permutation_idx_;
        if (ele->boundary_idx_) delete[] ele->boundary_idx_;
        if (ele->neigh_vb) delete[] ele->neigh_vb;
    }

    if (row_4_el != nullptr) delete[] row_4_el;
    if (el_4_loc != nullptr) delete[] el_4_loc;
    if (el_ds != nullptr) delete el_ds;
}


unsigned int Mesh::n_sides()
{
    if (n_sides_ == NDEF) {
        n_sides_=0;
        FOR_ELEMENTS(this, ele) n_sides_ += ele->n_sides();
    }
    return n_sides_;
}

unsigned int Mesh::n_corners() {
    unsigned int li, count = 0;
    FOR_ELEMENTS(this, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            count++;
        }
    }
    return count;
}

Partitioning *Mesh::get_part() {
    return part_.get();
}


//=============================================================================
// COUNT ELEMENT TYPES
//=============================================================================

void Mesh::count_element_types() {
    FOR_ELEMENTS(this, elm)
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
	for (auto elem_to_region : map) {
		vector<Element>::iterator ele = element_vec_.begin() + elem_index(elem_to_region.first);
		ele->region_idx_ = region_db_.get_region( elem_to_region.second, ele->dim() );
		region_db_.mark_used_region(ele->region_idx_.idx());
	}
}



void Mesh::setup_topology() {
    START_TIMER("MESH - setup topology");
    
    count_element_types();

    // check mesh quality
    FOR_ELEMENTS(this, ele)
        if (ele->quality_measure_smooth() < 0.001) WarningOut().fmt("Bad quality (<0.001) of the element {}.\n", ele->id());

    make_neighbours_and_edges();
    element_to_neigh_vb();
    make_edge_permutations();
    count_side_types();
    
    part_ = std::make_shared<Partitioning>(this, in_record_.val<Input::Record>("partitioning") );

    // create parallel distribution and numbering of elements
    IdxInt *id_4_old = new IdxInt[n_elements()];
    int i = 0;
    FOR_ELEMENTS(this, ele)
        id_4_old[i++] = elem_index( ele->id() );
    part_->id_maps(n_elements(), id_4_old, el_ds, el_4_loc, row_4_el);

    delete[] id_4_old;
}


//
void Mesh::count_side_types()
{

    n_insides = 0;
    n_exsides = 0;
    FOR_SIDES(this,  sde ) {
        if (sde->is_external()) n_exsides++;
        else n_insides++;
    }
}



void Mesh::create_node_element_lists() {
    // for each node we make a list of elements that use this node
    node_elements_.resize(node_vector.size());

    FOR_ELEMENTS( this, e )
        for (unsigned int n=0; n<e->n_nodes(); n++)
            node_elements_[node_vector.index(e->node[n])].push_back(e->index());

    for (vector<vector<unsigned int> >::iterator n=node_elements_.begin(); n!=node_elements_.end(); n++)
        stable_sort(n->begin(), n->end());
}


void Mesh::intersect_element_lists(vector<unsigned int> const &nodes_list, vector<unsigned int> &intersection_element_list)
{
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
    for( vector<unsigned int>::iterator ele = element_list.begin(); ele!=element_list.end(); ++ele)
        if (element_vec_[*ele].dim() == dim) { // keep only indexes of elements of same dimension
            *e_dest=*ele;
            ++e_dest;
        } else if (element_vec_[*ele].dim() == dim-1) { // get only first element of lower dimension
            if (is_neighbour) xprintf(UsrErr, "Too matching elements id: %d and id: %d in the same mesh.\n",
            		element_vec_[*ele].id(), element_vec_[element_idx].id() );

            is_neighbour = true;
            element_idx = *ele;
        }
    element_list.resize( e_dest - element_list.begin());
    return is_neighbour;
}

bool Mesh::same_sides(const SideIter &si, vector<unsigned int> &side_nodes) {
    // check if nodes lists match (this is slow and will be faster only when we convert whole mesh into hierarchical design like in deal.ii)
    unsigned int ni=0;
    while ( ni < si->n_nodes()
        && find(side_nodes.begin(), side_nodes.end(), node_vector.index( si->node(ni) ) ) != side_nodes.end() ) ni++;
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
    Neighbour neighbour;
    Edge *edg;
    unsigned int ngh_element_idx, last_edge_idx;

    create_node_element_lists();

	// pointers to created edges
	//vector<Edge *> tmp_edges;
    edges.resize(0); // be sure that edges are empty

	vector<unsigned int> side_nodes;
	vector<unsigned int> intersection_list; // list of elements in intersection of node element lists

	auto bc_range = get_element_range(true);
	for( unsigned int i=bc_range.first; i<bc_range.second; ++i) {
		Element &bc_ele = element_vec_[i];
        // Find all elements that share this side.
        side_nodes.resize(bc_ele.n_nodes());
        for (unsigned n=0; n<bc_ele.n_nodes(); n++) side_nodes[n] = node_vector.index(bc_ele.node[n]);
        intersect_element_lists(side_nodes, intersection_list);
        bool is_neighbour = find_lower_dim_element(intersection_list, bc_ele.dim() +1, ngh_element_idx);
        if (is_neighbour) {
            xprintf(UsrErr, "Boundary element (id: %d) match a regular element (id: %d) of lower dimension.\n",
                    bc_ele.id(), element_vec_[ngh_element_idx].id());
        } else {
            if (intersection_list.size() == 0) {
                // no matching dim+1 element found
            	WarningOut().fmt("Lonely boundary element, id: {}, region: {}, dimension {}.\n",
            			bc_ele.id(), bc_ele.region().id(), bc_ele.dim());
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
            Boundary &bdr=boundary_.back();
            bdr.bc_ele_idx_ = i;
            bdr.edge_idx_ = last_edge_idx;
            bdr.mesh_=this;

            // for 1d boundaries there can be more then one 1d elements connected to the boundary element
            // we do not detect this case later in the main search over bulk elements
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect)  {
                Element *elem = &(element_vec_[*isect]);
                for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++) {
                    SideIter si = elem->side(ecs);
                    if ( same_sides( si, side_nodes) ) {
                        if (elem->edge_idx_[ecs] != Mesh::undef_idx) {
                        	OLD_ASSERT(elem->boundary_idx_!=nullptr, "Null boundary idx array.\n");
                            int last_bc_ele_idx=this->boundary_[elem->boundary_idx_[ecs]].bc_ele_idx_;
                            int new_bc_ele_idx=i;
                            THROW( ExcDuplicateBoundary()
                                    << EI_ElemLast(this->element_vec_[last_bc_ele_idx].id())
                                    << EI_RegLast(this->element_vec_[last_bc_ele_idx].region().label())
                                    << EI_ElemNew(this->element_vec_[new_bc_ele_idx].id())
                                    << EI_RegNew(this->element_vec_[new_bc_ele_idx].region().label())
                                    );
                        }
                        elem->edge_idx_[ecs] = last_edge_idx;
                        edg->side_[ edg->n_sides++ ] = si;

                        if (elem->boundary_idx_ == NULL) {
                            elem->boundary_idx_ = new unsigned int [ elem->n_sides() ];
                            std::fill( elem->boundary_idx_, elem->boundary_idx_ + elem->n_sides(), Mesh::undef_idx);
                        }
                        elem->boundary_idx_[ecs] = bdr_idx;
                        break; // next element in intersection list
                    }
                }
            }

        }

	}
	// Now we go through all element sides and create edges and neighbours
	FOR_ELEMENTS( this, e )
	{
		for (unsigned int s=0; s<e->n_sides(); s++)
		{
			// skip sides that were already found
			if (e->edge_idx_[s] != Mesh::undef_idx) continue;


			// Find all elements that share this side.
			side_nodes.resize(e->side(s)->n_nodes());
			for (unsigned n=0; n<e->side(s)->n_nodes(); n++) side_nodes[n] = node_vector.index(e->side(s)->node(n));
			intersect_element_lists(side_nodes, intersection_list);

			bool is_neighbour = find_lower_dim_element(intersection_list, e->dim(), ngh_element_idx);

			if (is_neighbour) { // edge connects elements of different dimensions
			    neighbour.element_ = element_vec_.begin() + ngh_element_idx;
            } else { // edge connects only elements of the same dimension
                // Allocate the array of sides.
                last_edge_idx=edges.size();
                edges.resize(last_edge_idx+1);
                edg = &( edges.back() );
                edg->n_sides = 0;
                edg->side_ = new struct SideIter[ intersection_list.size() ];
                if (intersection_list.size() > max_edge_sides_[e->dim()-1])
                	max_edge_sides_[e->dim()-1] = intersection_list.size();

                if (intersection_list.size() == 1) { // outer edge, create boundary object as well
                    edg->n_sides=1;
                    edg->side_[0] = e->side(s);
                    e->edge_idx_[s] = last_edge_idx;

                    if (e->boundary_idx_ == NULL) {
                        e->boundary_idx_ = new unsigned int [ e->n_sides() ];
                        std::fill( e->boundary_idx_, e->boundary_idx_ + e->n_sides(), Mesh::undef_idx);
                    }

                    unsigned int bdr_idx=boundary_.size()+1; // need for VTK mesh that has no boundary elements
                                                             // and bulk elements are indexed from 0
                    boundary_.resize(bdr_idx+1);
                    Boundary &bdr=boundary_.back();
                    e->boundary_idx_[s] = bdr_idx;

                    // fill boundary element
                    vector<Element>::iterator bc_ele_iter = add_element_to_vector(-bdr_idx, true);
                    bc_ele_iter->init(e->dim()-1, -bdr_idx, this, region_db_.implicit_boundary_region() );
                    region_db_.mark_used_region( bc_ele_iter->region_idx_.idx() );
                    for(unsigned int ni = 0; ni< side_nodes.size(); ni++) bc_ele_iter->node[ni] = &( node_vector[side_nodes[ni]] );

                    // fill Boundary object
                    bdr.edge_idx_ = last_edge_idx;
                    bdr.bc_ele_idx_ = elem_index(-bdr_idx);
                    bdr.mesh_=this;

                    continue; // next side of element e
                }
			}

			// go through the elements connected to the edge or neighbour
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect) {
                Element *elem = &(element_vec_[*isect]);
                for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++) {
                    if (elem->edge_idx_[ecs] != Mesh::undef_idx) continue;
                    SideIter si = elem->side(ecs);
                    if ( same_sides( si, side_nodes) ) {
                        if (is_neighbour) {
                            // create a new edge and neighbour for this side, and element to the edge
                            last_edge_idx=edges.size();
                            edges.resize(last_edge_idx+1);
                            edg = &( edges.back() );
                            edg->n_sides = 1;
                            edg->side_ = new struct SideIter[1];
                            edg->side_[0] = si;
                            elem->edge_idx_[ecs] = last_edge_idx;

                            neighbour.edge_idx_ = last_edge_idx;

                            vb_neighbours_.push_back(neighbour); // copy neighbour with this edge setting
                        } else {
                            // connect the side to the edge, and side to the edge
                            edg->side_[ edg->n_sides++ ] = si;
                            elem->edge_idx_[ecs] = last_edge_idx;
                        }
                        break; // next element from intersection list
                    }
                } // search for side of other connected element
            } // connected elements
            OLD_ASSERT( is_neighbour || ( (unsigned int) edg->n_sides ) == intersection_list.size(), "Some connected sides were not found.\n");
		} // for element sides
	}   // for elements

	MessageOut().fmt( "Created {} edges and {} neighbours.\n", edges.size(), vb_neighbours_.size() );
}



void Mesh::make_edge_permutations()
{
	for (EdgeVector::iterator edg=edges.begin(); edg!=edges.end(); edg++)
	{
		// side 0 is reference, so its permutation is 0
		edg->side(0)->element()->permutation_idx_[edg->side(0)->el_idx()] = 0;

		if (edg->n_sides > 1)
		{
			map<const Node*,unsigned int> node_numbers;
			unsigned int permutation[edg->side(0)->n_nodes()];

			for (unsigned int i=0; i<edg->side(0)->n_nodes(); i++)
				node_numbers[edg->side(0)->node(i)] = i;

			for (int sid=1; sid<edg->n_sides; sid++)
			{
				for (unsigned int i=0; i<edg->side(0)->n_nodes(); i++)
					permutation[node_numbers[edg->side(sid)->node(i)]] = i;

				switch (edg->side(0)->dim())
				{
				case 0:
					edg->side(sid)->element()->permutation_idx_[edg->side(sid)->el_idx()] = RefElement<1>::permutation_index(permutation);
					break;
				case 1:
					edg->side(sid)->element()->permutation_idx_[edg->side(sid)->el_idx()] = RefElement<2>::permutation_index(permutation);
					break;
				case 2:
					edg->side(sid)->element()->permutation_idx_[edg->side(sid)->el_idx()] = RefElement<3>::permutation_index(permutation);
					break;
				}
			}
		}
	}

	for (vector<Neighbour>::iterator nb=vb_neighbours_.begin(); nb!=vb_neighbours_.end(); nb++)
	{
		map<const Node*,unsigned int> node_numbers;
		unsigned int permutation[nb->element()->n_nodes()];

		// element of lower dimension is reference, so
		// we calculate permutation for the adjacent side
		for (unsigned int i=0; i<nb->element()->n_nodes(); i++)
			node_numbers[nb->element()->node[i]] = i;

		for (unsigned int i=0; i<nb->side()->n_nodes(); i++)
			permutation[node_numbers[nb->side()->node(i)]] = i;

		switch (nb->side()->dim())
		{
		case 0:
			nb->side()->element()->permutation_idx_[nb->side()->el_idx()] = RefElement<1>::permutation_index(permutation);
			break;
		case 1:
			nb->side()->element()->permutation_idx_[nb->side()->el_idx()] = RefElement<2>::permutation_index(permutation);
			break;
		case 2:
			nb->side()->element()->permutation_idx_[nb->side()->el_idx()] = RefElement<3>::permutation_index(permutation);
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

    FOR_ELEMENTS(this,ele) ele->n_neighs_vb =0;

    // count vb neighs per element
    FOR_NEIGHBOURS(this,  ngh )  ngh->element_->n_neighs_vb++;

    // Allocation of the array per element
    FOR_ELEMENTS(this,  ele )
        if( ele->n_neighs_vb > 0 ) {
            ele->neigh_vb = new struct Neighbour* [ele->n_neighs_vb];
            ele->n_neighs_vb=0;
        }

    // fill
    ElementIterator ele;
    FOR_NEIGHBOURS(this,  ngh ) {
        ele = ngh->element();
        ele->neigh_vb[ ele->n_neighs_vb++ ] = &( *ngh );
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



ElementAccessor<3> Mesh::element_accessor(unsigned int idx, bool boundary) {
    return ElementAccessor<3>(this, idx, boundary);
}



void Mesh::elements_id_maps( vector<IdxInt> & bulk_elements_id, vector<IdxInt> & boundary_elements_id) const
{
    if (bulk_elements_id.size() ==0) {
        std::vector<IdxInt>::iterator map_it;
        IdxInt last_id;

        bulk_elements_id.resize(n_elements());
        map_it = bulk_elements_id.begin();
        last_id = -1;
        for(unsigned int idx=0; idx < n_elements(); idx++, ++map_it) {
        	IdxInt id = element_vec_[idx].id();
            if (last_id >= id) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", id);
            last_id=*map_it = id;
        }

        boundary_elements_id.resize(n_elements(true));
        map_it = boundary_elements_id.begin();
        last_id = -1;
        auto bc_range = get_element_range(true);
        for(unsigned int idx=bc_range.first; idx<bc_range.second; idx++, ++map_it) {
        	IdxInt id = element_vec_[idx].id();
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


void Mesh::compute_element_boxes() {
    START_TIMER("Mesh::compute_element_boxes");
    if (element_box_.size() > 0) return;

    // make element boxes
    element_box_.resize(this->n_elements());
    unsigned int i=0;
    FOR_ELEMENTS(this, element) {
         element_box_[i] = element->bounding_box();
         i++;
    }

    // make mesh box
    Node* node = this->node_vector.begin();
    mesh_box_ = BoundingBox(node->point(), node->point());
    FOR_NODES(this, node ) {
        mesh_box_.expand( node->point() );
    }

}

const BIHTree &Mesh::get_bih_tree() {
    if (! this->bih_tree_)
        bih_tree_ = std::make_shared<BIHTree>(this);
    return *bih_tree_;
}

double Mesh::global_observe_radius() const {
	return in_record_.val<double>("global_observe_search_radius");
}

void Mesh::add_physical_name(unsigned int dim, unsigned int id, std::string name) {
	region_db_.add_region(id, name, dim, "$PhysicalNames");
}


void Mesh::add_node(unsigned int node_id, arma::vec3 coords) {
	NodeFullIter node = node_vector.add_item(node_id);
	node->point() = coords;
}


void Mesh::add_element(unsigned int elm_id, unsigned int dim, unsigned int region_id, unsigned int partition_id,
		std::vector<unsigned int> node_ids) {
	vector<Element>::iterator ele_iter;
	RegionIdx region_idx = region_db_.get_region( region_id, dim );
	if ( !region_idx.is_valid() ) {
		region_idx = region_db_.add_region( region_id, region_db_.create_label_from_id(region_id), dim, "$Element" );
	}
	region_db_.mark_used_region(region_idx.idx());

	if (region_idx.is_boundary()) {
		ele_iter = add_element_to_vector(elm_id, true);
	} else {
		if(dim == 0 ) {
			WarningOut().fmt("Bulk elements of zero size(dim=0) are not supported. Element ID: {}.\n", elm_id);
			return;
		}
		else {
			ele_iter = add_element_to_vector(elm_id);
		}
	}
	ele_iter->init(dim, elm_id, this, region_idx);
	ele_iter->pid = partition_id;

	unsigned int ni;
	FOR_ELEMENT_NODES(&*ele_iter, ni) {
		unsigned int node_id = node_ids[ni];
		NodeIter node = node_vector.find_id( node_id );
		INPUT_CHECK( node != node_vector.end(),
				"Unknown node id %d in specification of element with id=%d.\n", node_id, elm_id);
		ele_iter->node[ni] = node;
	}

    // check that tetrahedron element is numbered correctly and is not degenerated
    if(ele_iter->dim() == 3)
    {
        double jac = ele_iter->tetrahedron_jacobian();
        if( ! (jac > 0) )
            WarningOut().fmt("Tetrahedron element with id {} has wrong numbering or is degenerated (Jacobian = {}).",ele_iter->id(),jac);
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
	element_vec_.resize(size);
	element_id_map_.clear();
	bulk_size_ = 0;
	boundary_size_ = 0;
}


vector<Element>::iterator Mesh::add_element_to_vector(int id, bool boundary) {
	ASSERT( element_id_map_.find(id) == element_id_map_.end())(id).error("Can not add item with given id since it already exists.");

	vector<Element>::iterator ele_iter;
	if (boundary) {
		if (id>=0) {
			unsigned int boundary_pos = element_vec_.size() - boundary_size_ - 1;
			ele_iter = element_vec_.begin()+boundary_pos;
			element_id_map_[id]=boundary_pos;
		} else {
            element_vec_.push_back( Element() );
            ele_iter = element_vec_.end()-1;
			element_id_map_[id]=element_vec_.size()-1;
		}
		boundary_size_++;
	} else {
		ele_iter = element_vec_.begin()+bulk_size_;
		element_id_map_[id]=bulk_size_;
		bulk_size_++;
	}

	return ele_iter;
}

//-----------------------------------------------------------------------------
// vim: set cindent:
