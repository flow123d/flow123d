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
#include "system/xio.h"
#include "input/reader_to_storage.hh"
#include "input/input_type.hh"
#include "system/sys_profiler.hh"
#include "la/distribution.hh"

#include <boost/tokenizer.hpp>
#include "boost/lexical_cast.hpp"
#include <boost/make_shared.hpp>

#include "mesh/mesh.h"
#include "mesh/ref_element.hh"

// think about following dependencies
#include "mesh/boundaries.h"
#include "mesh/accessors.hh"
#include "mesh/partitioning.hh"

#include "mesh/bih_tree.hh"


//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "mesh/msh_gmshreader.h"
#include "mesh/region.hh"

#define NDEF  -1

namespace IT = Input::Type;


const IT::Record & Mesh::get_input_type() {
	return IT::Record("Mesh","Record with mesh related data." )
		.declare_key("mesh_file", IT::FileName::input(), IT::Default::obligatory(),
				"Input file with mesh description.")
		.declare_key("regions", IT::Array( RegionDB::get_region_input_type() ), IT::Default::optional(),
				"List of additional region definitions not contained in the mesh.")
		.declare_key("sets", IT::Array( RegionDB::get_region_set_input_type()), IT::Default::optional(),
				"List of region set definitions. There are three region sets implicitly defined:\n\n"
				" - ALL (all regions of the mesh)\n - BOUNDARY (all boundary regions)\n - and BULK (all bulk regions)")
		.declare_key("partitioning", Partitioning::get_input_type(), IT::Default("\"any_neighboring\""), "Parameters of mesh partitioning algorithms.\n" )
		.close();
}



const unsigned int Mesh::undef_idx;

Mesh::Mesh(const std::string &input_str, MPI_Comm comm)
:comm_(comm),
 row_4_el(nullptr),
 el_ds(nullptr),
 el_4_loc(nullptr)
{

    Input::ReaderToStorage reader( input_str, Mesh::get_input_type(), Input::FileFormat::format_JSON );
    in_record_ = reader.get_root_interface<Input::Record>();

    reinit(in_record_);
}



Mesh::Mesh(Input::Record in_record, MPI_Comm com)
: in_record_(in_record),
  comm_(com),
  row_4_el(nullptr),
  el_ds(nullptr),
  el_4_loc(nullptr)
{
    reinit(in_record_);
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
    		side_nodes[0][sid][nid] = RefElement<1>::side_nodes[sid][nid];

    for (unsigned int sid=0; sid<RefElement<2>::n_sides; sid++)
        	for (unsigned int nid=0; nid<RefElement<2>::n_nodes_per_side; nid++)
        		side_nodes[1][sid][nid] = RefElement<2>::side_nodes[sid][nid];

    for (unsigned int sid=0; sid<RefElement<3>::n_sides; sid++)
        	for (unsigned int nid=0; nid<RefElement<3>::n_nodes_per_side; nid++)
        		side_nodes[2][sid][nid] = RefElement<3>::side_nodes[sid][nid];
}


Mesh::~Mesh() {
    for(Edge &edg : this->edges)
        if (edg.side_) delete[] edg.side_;

    FOR_ELEMENTS( this, ele ) {
        if (ele->node) delete[] ele->node;
        if (ele->edge_idx_) delete[] ele->edge_idx_;
        if (ele->permutation_idx_) delete[] ele->permutation_idx_;
        if (ele->boundary_idx_) delete[] ele->boundary_idx_;
    }

    for(unsigned int idx=0; idx < this->bc_elements.size(); idx++) {
        Element *ele=&(bc_elements[idx]);
        if (ele->node) delete[] ele->node;
        if (ele->edge_idx_) delete[] ele->edge_idx_;
        if (ele->permutation_idx_) delete[] ele->permutation_idx_;
        if (ele->boundary_idx_) delete[] ele->boundary_idx_;
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


void Mesh::read_gmsh_from_stream(istream &in) {
  
    START_TIMER("Reading mesh - from_stream");
    
    GmshMeshReader reader(in);
    reader.read_mesh(this);
    setup_topology();
    //close region_db_.
    region_db_.close();
}



void Mesh::init_from_input() {
    START_TIMER("Reading mesh - init_from_input");
    
    Input::Array region_list;
    RegionDB::MapElementIDToRegionID el_to_reg_map;

    // read raw mesh, add regions from GMSH file
    GmshMeshReader reader( in_record_.val<FilePath>("mesh_file") );
    reader.read_mesh(this);
    // possibly add implicit_boundary region.
    setup_topology();
    // create regions from our input
    if (in_record_.opt_val("regions", region_list)) {
        region_db_.read_regions_from_input(region_list, el_to_reg_map);
    }
    modify_element_ids(el_to_reg_map);
    //close region_db_.
    region_db_.close();
    // create sets
    Input::Array set_list;
    if (in_record_.opt_val("sets", set_list)) {
        region_db_.read_sets_from_input(set_list);
    }
}




void Mesh::modify_element_ids(const RegionDB::MapElementIDToRegionID &map) {
	for (auto elem_to_region : map) {
		ElementIter ele = this->element.find_id(elem_to_region.first);
		ele->region_idx_ = region_db_.add_region( elem_to_region.second, ele->dim() );
	}
}




void Mesh::setup_topology() {
    START_TIMER("MESH - setup topology");
    
    count_element_types();

    // check mesh quality
    FOR_ELEMENTS(this, ele)
        if (ele->quality_measure_smooth() < 0.001) xprintf(Warn, "Bad quality (<0.001) of the element %u.\n", ele.id());

    make_neighbours_and_edges();
    element_to_neigh_vb();
    make_edge_permutations();
    count_side_types();

    part_ = boost::make_shared<Partitioning>(this, in_record_.val<Input::Record>("partitioning") );

    // create parallel distribution and numbering of elements
    int *id_4_old = new int[element.size()];
    int i = 0;
    FOR_ELEMENTS(this, ele) id_4_old[i++] = ele.index();
    part_->id_maps(element.size(), id_4_old, el_ds, el_4_loc, row_4_el);
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
    node_elements.resize(node_vector.size());

    FOR_ELEMENTS( this, e )
        for (unsigned int n=0; n<e->n_nodes(); n++)
            node_elements[node_vector.index(e->node[n])].push_back(e->index());

    for (vector<vector<unsigned int> >::iterator n=node_elements.begin(); n!=node_elements.end(); n++)
        stable_sort(n->begin(), n->end());
}


void Mesh::intersect_element_lists(vector<unsigned int> const &nodes_list, vector<unsigned int> &intersection_element_list)
{
    if (nodes_list.size() == 0) {
        intersection_element_list.clear();
    } else if (nodes_list.size() == 1) {
        intersection_element_list = node_elements[ nodes_list[0] ];
	} else {
	    vector<unsigned int>::const_iterator it1=nodes_list.begin();
	    vector<unsigned int>::const_iterator it2=it1+1;
	    intersection_element_list.resize( node_elements[*it1].size() ); // make enough space

	    it1=set_intersection(
                node_elements[*it1].begin(), node_elements[*it1].end(),
                node_elements[*it2].begin(), node_elements[*it2].end(),
                intersection_element_list.begin());
        intersection_element_list.resize(it1-intersection_element_list.begin()); // resize to true size

        for(;it2<nodes_list.end();++it2) {
            it1=set_intersection(
                    intersection_element_list.begin(), intersection_element_list.end(),
                    node_elements[*it2].begin(), node_elements[*it2].end(),
                    intersection_element_list.begin());
            intersection_element_list.resize(it1-intersection_element_list.begin()); // resize to true size
        }
    }
}


bool Mesh::find_lower_dim_element( ElementVector &elements, vector<unsigned int> &element_list,
        unsigned int dim, unsigned int &element_idx) {
    bool is_neighbour = false;

    vector<unsigned int>::iterator e_dest=element_list.begin();
    for( vector<unsigned int>::iterator ele = element_list.begin(); ele!=element_list.end(); ++ele)
        if (elements[*ele].dim() == dim) { // keep only indexes of elements of same dimension
            *e_dest=*ele;
            ++e_dest;
        } else if (elements[*ele].dim() == dim-1) { // get only first element of lower dimension
            if (is_neighbour) xprintf(UsrErr, "Too matching elements id: %d and id: %d in the same mesh.\n",
                    elements(*ele).id(), elements(element_idx).id() );

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
 * - process bc_elements first, there should be no Neigh, but check it
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

	for( ElementFullIter bc_ele = bc_elements.begin(); bc_ele != bc_elements.end(); ++bc_ele) {
        // Find all elements that share this side.
        side_nodes.resize(bc_ele->n_nodes());
        for (unsigned n=0; n<bc_ele->n_nodes(); n++) side_nodes[n] = node_vector.index(bc_ele->node[n]);
        intersect_element_lists(side_nodes, intersection_list);
        bool is_neighbour = find_lower_dim_element(element, intersection_list, bc_ele->dim() +1, ngh_element_idx);
        if (is_neighbour) {
            xprintf(UsrErr, "Boundary element (id: %d) match a regular element (id: %d) of lower dimension.\n",
                    bc_ele.id(), element(ngh_element_idx).id());
        } else {
            if (intersection_list.size() == 0) {
                // no matching dim+1 element found
                xprintf(Warn, "Lonely boundary element, id: %d, region: %d, dimension %d.\n", bc_ele.id(), bc_ele->region().id(), bc_ele->dim());
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
            bdr.bc_ele_idx_ = bc_ele.index();
            bdr.edge_idx_ = last_edge_idx;
            bdr.mesh_=this;

            // for 1d boundaries there can be more then one 1d elements connected to the boundary element
            // we do not detect this case later in the main search over bulk elements
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect)  {
                Element *elem = &(element[*isect]);
                for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++) {
                    SideIter si = elem->side(ecs);
                    if ( same_sides( si, side_nodes) ) {
                        if (elem->edge_idx_[ecs] != Mesh::undef_idx) {
                            ASSERT(elem->boundary_idx_!=nullptr, "Null boundary idx array.\n");
                            int last_bc_ele_idx=this->boundary_[elem->boundary_idx_[ecs]].bc_ele_idx_;
                            int new_bc_ele_idx=bc_ele.index();
                            THROW( ExcDuplicateBoundary()
                                    << EI_ElemLast(this->bc_elements.get_id(last_bc_ele_idx))
                                    << EI_RegLast(this->bc_elements[last_bc_ele_idx].region().label())
                                    << EI_ElemNew(this->bc_elements.get_id(new_bc_ele_idx))
                                    << EI_RegNew(this->bc_elements[new_bc_ele_idx].region().label())
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

			bool is_neighbour = find_lower_dim_element(element, intersection_list, e->dim(), ngh_element_idx);

			if (is_neighbour) { // edge connects elements of different dimensions
			    neighbour.element_ = &(element[ngh_element_idx]);
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

                    unsigned int bdr_idx=boundary_.size();
                    boundary_.resize(bdr_idx+1);
                    Boundary &bdr=boundary_.back();
                    e->boundary_idx_[s] = bdr_idx;

                    // fill boundary element
                    ElementFullIter bc_ele = bc_elements.add_item( -bdr_idx ); // use negative bcd index as ID,
                    bc_ele->init(e->dim()-1, this, region_db_.implicit_boundary_region() );
                    for(unsigned int ni = 0; ni< side_nodes.size(); ni++) bc_ele->node[ni] = &( node_vector[side_nodes[ni]] );

                    // fill Boundary object
                    bdr.edge_idx_ = last_edge_idx;
                    bdr.bc_ele_idx_ = bc_ele.index();
                    bdr.mesh_=this;

                    continue; // next side of element e
                }
			}

			// go through the elements connected to the edge or neighbour
            for( vector<unsigned int>::iterator isect = intersection_list.begin(); isect!=intersection_list.end(); ++isect) {
                Element *elem = &(element[*isect]);
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
            ASSERT( is_neighbour || ( (unsigned int) edg->n_sides ) == intersection_list.size(), "Some connected sides were not found.\n");
		} // for element sides
	}   // for elements

	xprintf( Msg, "Created %d edges and %d neighbours.\n", edges.size(), vb_neighbours_.size() );
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

    xprintf( MsgVerb, "   Element to neighbours of vb2 type... ")/*orig verb 5*/;

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
    ElementIter ele;
    FOR_NEIGHBOURS(this,  ngh ) {
        ele = ngh->element();
        ele->neigh_vb[ ele->n_neighs_vb++ ] = &( *ngh );
    }

    xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}




#include "mesh/ngh/include/triangle.h"
#include "mesh/ngh/include/abscissa.h"
#include "mesh/ngh/include/intersection.h"


void Mesh::make_intersec_elements() {
	/* Algorithm:
	 *
	 * 1) create BIH tree
	 * 2) for every 1D, find list of candidates
	 * 3) compute intersections for 1d, store it to master_elements
	 *
	 */
	BIHTree bih_tree( this );
	master_elements.resize(n_elements());

	for(unsigned int i_ele=0; i_ele<n_elements(); i_ele++) {
		Element &ele = this->element[i_ele];

		if (ele.dim() == 1) {
			vector<unsigned int> candidate_list;
                        bih_tree.find_bounding_box(ele.bounding_box(), candidate_list);
                        
			//for(unsigned int i_elm=0; i_elm<n_elements(); i_elm++) {
                        for(unsigned int i_elm : candidate_list) {
				ElementFullIter elm = this->element( i_elm );
				if (elm->dim() == 2) {
					IntersectionLocal *intersection;
					GetIntersection( TAbscissa(ele), TTriangle(*elm), intersection);
					if (intersection && intersection->get_type() == IntersectionLocal::line) {

						master_elements[i_ele].push_back( intersections.size() );
						intersections.push_back( Intersection(this->element(i_ele), elm, intersection) );
				    }
				}

			}
		}
	}

}



ElementAccessor<3> Mesh::element_accessor(unsigned int idx, bool boundary) {
    return ElementAccessor<3>(this, idx, boundary);
}



vector<int> const & Mesh::elements_id_maps( bool boundary_domain) const
{
    if (bulk_elements_id_.size() ==0) {
        std::vector<int>::iterator map_it;
        int last_id;

        bulk_elements_id_.resize(n_elements());
        map_it = bulk_elements_id_.begin();
        last_id = -1;
        for(unsigned int idx=0; idx < element.size(); idx++, ++map_it) {
        	int id = element.get_id(idx);
            if (last_id >= id) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", id);
            last_id=*map_it = id;
        }

        boundary_elements_id_.resize(bc_elements.size());
        map_it = boundary_elements_id_.begin();
        last_id = -1;
        for(unsigned int idx=0; idx < bc_elements.size(); idx++, ++map_it) {
        	int id = bc_elements.get_id(idx);
            // We set ID for boundary elements created by the mesh itself to "-1"
            // this force gmsh reader to skip all remaining entries in boundary_elements_id_
            // and thus report error for any remaining data lines
            if (id < 0) last_id=*map_it=-1;
            else {
                if (last_id >= id) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", id);
                last_id=*map_it = id;
            }
        }
    }

    if (boundary_domain) return boundary_elements_id_;
    return bulk_elements_id_;
}

//-----------------------------------------------------------------------------
// vim: set cindent:
