/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @ingroup mesh
 * @brief  Mesh construction
 *
 */

#include <unistd.h>
#include <set>


#include "system/system.hh"
#include "system/xio.h"
#include "input/input_type.hh"
#include "system/sys_profiler.hh"

#include <boost/tokenizer.hpp>
#include "boost/lexical_cast.hpp"

#include "mesh/mesh.h"
#include "fem/ref_element.hh"

// think about following dependencies
#include "mesh/boundaries.h"
#include "mesh/accessors.hh"


//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "mesh/msh_reader.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/region.hh"


using namespace Input::Type;


Record Mesh::input_type
	= Record("Mesh","Record with mesh related data." )
	.declare_key("mesh_file", FileName::input(), Default::obligatory(),
			"Input file with mesh description.")
	.declare_key("regions", Array( RegionDB::region_input_type ), Default::optional(),
	        "List of additional region definitions not contained in the mesh.")
	.declare_key("sets", Array( RegionDB::region_set_input_type), Default::optional(),
	        "List of region set definitions. There are three region sets implicitly defined:\n"
	        "ALL (all regions of the mesh), BOUNDARY (all boundary regions), and BULK (all bulk regions)")
	.close();



const unsigned int Mesh::undef_idx;

Mesh::Mesh()
{
    reinit(in_record_);
}



Mesh::Mesh(Input::Record in_record)
: in_record_(in_record) {
    reinit(in_record_);
}



void Mesh::reinit(Input::Record in_record)
{

    n_materials = NDEF;

    n_insides = NDEF;
    n_exsides = NDEF;
    n_sides_ = NDEF;

    // number of element of particular dimension
    n_lines = 0;
    n_triangles = 0;
    n_tetrahedras = 0;


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


unsigned int Mesh::n_sides()
{
    if (n_sides_ == NDEF) {
        n_sides_=0;
        FOR_ELEMENTS(this, ele) n_sides_ += ele->n_sides();
    }
    return n_sides_;
}


//=============================================================================
// COUNT ELEMENT TYPES
//=============================================================================

void Mesh::count_element_types() {
    F_ENTRY;

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
  
    START_TIMER("READING MESH - from_stream");
    
    GmshMeshReader reader(in);
    reader.read_mesh(this);
    setup_topology();
}



void Mesh::init_from_input() {
    F_ENTRY;

    START_TIMER("READING MESH - init_from_input");
    
    Input::Array region_list;
    RegionDB::MapElementIDToRegionID el_to_reg_map;

    // create regions from our input
    if (in_record_.opt_val("regions", region_list)) {
        region_db_.read_regions_from_input(region_list, el_to_reg_map);
    }
    // read raw mesh, add regions from GMSH file
    GmshMeshReader reader( in_record_.val<FilePath>("mesh_file") );
    reader.read_mesh(this, &el_to_reg_map);
    // possibly add implicit_boundary region, close region_db_.
    setup_topology();
    // create sets
    Input::Array set_list;
    if (in_record_.opt_val("sets", set_list)) {
        region_db_.read_sets_from_input(set_list);
    }
}




void Mesh::setup_topology() {
    F_ENTRY;

    START_TIMER("setup topology");
    
    count_element_types();
    make_neighbours_and_edges();
    element_to_neigh_vb();
    make_edge_permutations();
    count_side_types();

    region_db_.close();
}


//
void Mesh::count_side_types()
{
    struct Side *sde;

    n_insides = 0;
    n_exsides = 0;
    //FOR_SIDES(this,  sde ) {
    //    DBGMSG( "ele: %d edge: %d\n", sde->element().index(), sde->edge_idx());
    //}
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
                xprintf(Warn, "Lonely boundary element, id: %d, dimension %d.\n", bc_ele.id(), bc_ele->dim());
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
                        edg->side_[ edg->n_sides++ ] = si;
                        elem->edge_idx_[ecs] = last_edge_idx;

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
                neighbour.sigma = 1;
            } else { // edge connects only elements of the same dimension
                // Allocate the array of sides.
                last_edge_idx=edges.size();
                edges.resize(last_edge_idx+1);
                edg = &( edges.back() );
                edg->n_sides = 0;
                edg->side_ = new struct SideIter[ intersection_list.size() ];

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
            ASSERT( is_neighbour || edg->n_sides == intersection_list.size(), "Some connected sides were not found.\n");
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
			const unsigned int dim = edg->side(0)->dim()+1;

			for (int i=0; i<edg->side(0)->n_nodes(); i++)
				node_numbers[edg->side(0)->node(i)] = i;

			for (int sid=1; sid<edg->n_sides; sid++)
			{
				for (int i=0; i<edg->side(0)->n_nodes(); i++)
					permutation[node_numbers[edg->side(sid)->node(i)]] = i;

				switch (dim)
				{
				case 1:
					edg->side(sid)->element()->permutation_idx_[edg->side(sid)->el_idx()] = RefElement<1>::permutation_index(permutation);
					break;
				case 2:
					edg->side(sid)->element()->permutation_idx_[edg->side(sid)->el_idx()] = RefElement<2>::permutation_index(permutation);
					break;
				case 3:
					edg->side(sid)->element()->permutation_idx_[edg->side(sid)->el_idx()] = RefElement<3>::permutation_index(permutation);
					break;
				}
			}
		}
	}
}





/**
 * Set Element->boundaries_ for all external sides. (temporary solution)
 */
void Mesh::create_external_boundary()
{
    /*
    // set to non zero all pointers including boundary connected to lower dim elements
    // these have only one side per edge
    Boundary empty_boundary;


    FOR_ELEMENTS(this, ele) {
        // is there any outer side
        bool outer=false;
        FOR_ELEMENT_SIDES(ele, si)
        {
            if ( ele->side(si)->edge()->n_sides == 1) {
                outer=true;
                break;
            }
        }
       if (outer) {
           // for elements on the boundary set boundaries_
           FOR_ELEMENT_SIDES(ele,si)
                if ( ele->side(si)->edge()->n_sides == 1)
                    ele->boundaries_[si] = &empty_boundary;
                else
                    ele->boundaries_[si] = NULL;

       } else {
           // can delete boundaries on internal elements !!
            delete ele->boundaries_;
            ele->boundaries_=NULL;
       }
    }

    int count=0;
    // pass through neighbours and set to NULL internal interfaces
    FOR_NEIGHBOURS(this,  ngh ) {
        SideIter s = ngh->side();
        if (s->element()->boundaries_ == NULL) continue;
        s->element()->boundaries_[ s->el_idx() ] = NULL;
    }

    // count remaining
    unsigned int n_boundaries=0;
    FOR_ELEMENTS(this, ele) {
        if (ele->boundaries_ == NULL) continue;
        FOR_ELEMENT_SIDES(ele, si)
            if (ele->boundaries_[si]) n_boundaries ++;
    }

    // fill boundaries
    BoundaryFullIter bcd(boundary);
    unsigned int ni;

    boundary.reserve(n_boundaries);
    DBGMSG("bc_elements size after read: %d\n", bc_elements.size());
    bc_elements.reserve(n_boundaries);
    FOR_ELEMENTS(this, ele) {
         if (ele->boundaries_ == NULL) continue;
         FOR_ELEMENT_SIDES(ele, si)
             if (ele->boundaries_[si]) {
                 // add boundary object
                 bcd = boundary.add_item();


                 // fill boundary element
                 Element * bc_ele = bc_elements.add_item( -bcd.index() ); // use negative bcd index as ID,
                 bc_ele->dim_ = ele->dim()-1;
                 bc_ele->node = new Node * [bc_ele->n_nodes()];
                 FOR_ELEMENT_NODES(bc_ele, ni) {
                     bc_ele->node[ni] = (Node *)ele->side(si)->node(ni);
                 }

                 // fill Boudary object
                 bcd->side = ele->side(si);
                 bcd->bc_element_ = bc_ele;
                 ele->boundaries_[si] = bcd;

             }
    }*/
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


/*
void Mesh::setup_materials( MaterialDatabase &base)
{
  
    xprintf( MsgVerb, "   Element to material... ");//orig verb 5;
    FOR_ELEMENTS(this, ele ) {
        ele->material=base.find_id(ele->region().id());
        INPUT_CHECK( ele->material != base.end(),
                "Reference to undefined material %d in element %d\n", ele->region().id(), ele.id() );
    }
    xprintf( MsgVerb, "O.K.\n");//orig verb 6
  
}
*/

void Mesh::read_intersections() {

    using namespace boost;

    ElementFullIter master(element), slave(element);

    char tmp_line[LINE_SIZE];
    string file_name = in_record_.val<FilePath>("neighbouring");
    FILE *in = xfopen( file_name , "rt" );

    tokenizer<boost::char_separator<char> >::iterator tok;

    xprintf( Msg, "Reading intersections...")/*orig verb 2*/;
    skip_to(in, "$Intersections");
    xfgets(tmp_line, LINE_SIZE - 2, in);
    int n_intersect = atoi(xstrtok(tmp_line));
    INPUT_CHECK( n_intersect >= 0 ,"Negative number of neighbours!\n");

    intersections.reserve(n_intersect);

    for (int i = 0; i < n_intersect; i++) {
        xfgets(tmp_line, LINE_SIZE - 2, in);
        string line = tmp_line;
        tokenizer<boost::char_separator<char> > line_tokenizer(line, boost::char_separator<char>("\t \n"));

        tok = line_tokenizer.begin();

        try {
            ++tok; // skip id token
            int type = lexical_cast<int> (*tok);
            ++tok;
            int master_id = lexical_cast<int> (*tok);
            ++tok;
            int slave_id = lexical_cast<int> (*tok);
            ++tok;
            double sigma = lexical_cast<double> (*tok);
            ++tok;

            int n_intersect_points = lexical_cast<int> (*tok);
            ++tok;
            master = element.find_id(master_id);
            slave = element.find_id(slave_id);

            intersections.push_back(Intersection(n_intersect_points - 1, master, slave, tok));
        } catch (bad_lexical_cast &) {
            xprintf(UsrErr, "Wrong number format at line %d in file %s x%sx\n",i, file_name.c_str(),(*tok).c_str());
        }

    }

    xprintf( Msg, "O.K.\n")/*orig verb 2*/;

}


void Mesh::make_intersec_elements() {

     // calculate sizes and make allocations
     vector<int >sizes(n_elements(),0);
     for( vector<Intersection>::iterator i=intersections.begin(); i != intersections.end(); ++i )
     sizes[i->master_iter().index()]++;
     master_elements.resize(n_elements());
     for(unsigned int i=0;i<n_elements(); ++i ) master_elements[i].reserve(sizes[i]);

     // fill intersec_elements
     for( vector<Intersection>::iterator i=intersections.begin(); i != intersections.end(); ++i )
     master_elements[i->master_iter().index()].push_back( i-intersections.begin() );

}



ElementAccessor<3> Mesh::element_accessor(unsigned int idx, bool boundary) {
    return ElementAccessor<3>(this, idx, boundary);
}



vector<int> const & Mesh::all_elements_id() {
    if (all_elements_id_.size() ==0) {

        all_elements_id_.resize(n_all_input_elements_);
        std::vector<int>::iterator all_it = all_elements_id_.begin();
        unsigned int last_id = element.begin().id();

        for(ElementFullIter it=element.begin(); it!=element.end(); ++it, ++all_it) {
            if (last_id > it.id()) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", it.id());
            last_id=*all_it = it.id();
        }
        for(ElementFullIter it=bc_elements.begin(); all_it!=all_elements_id_.end(); ++it, ++all_it) {
            if (last_id > it.id()) xprintf(UsrErr, "Element IDs in non-increasing order, ID: %d\n", it.id());
            last_id=*all_it = it.id();
        }
    }
    return all_elements_id_;
}


/*
void Mesh::make_edge_list_from_neigh() {
    int edi;
    Mesh *mesh = this;
    struct Neighbour *ngh;

    xprintf( Msg, "Creating edges from neigbours... ");

    int n_edges = mesh->n_sides;
    FOR_NEIGHBOURS( ngh )
        if (ngh->type == BB_E || ngh->type == BB_EL)
            n_edges-=( ngh->n_elements - 1 );

    mesh->edge.resize(n_edges);

    xprintf( MsgVerb, " O.K. %d edges created.", mesh->n_edges());

    EdgeFullIter edg = mesh->edge.begin();
    n_edges=0;
    FOR_NEIGHBOURS( ngh )
        if (ngh->type == BB_E || ngh->type == BB_EL) {
            ngh->edge = edg;
            edg->neigh_bb = ngh;
            ++edg;
            n_edges++;
        }
    xprintf( MsgVerb, "O.K. %d\n");

}*/
//-----------------------------------------------------------------------------
// vim: set cindent:
