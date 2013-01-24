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

#include <boost/tokenizer.hpp>
#include "boost/lexical_cast.hpp"

#include "mesh/mesh.h"

// think about following dependencies
#include "mesh/boundaries.h"
#include "mesh/accessors.hh"


//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "mesh/topology.cc"
#include "mesh/msh_reader.h"
#include "mesh/msh_gmshreader.h"
#include "mesh/region.hh"

void count_element_types(Mesh*);
void read_node_list(Mesh*);


using namespace Input::Type;

Record BoundarySegment::input_type
    = Record("BoundarySegment","Record with specification of boundary segments,\n"
            "i.e. subsets of the domain boundary where we prescribe one boundary condition.\n"
            "Currently very GMSH oriented. (NOT IMPLEMENTED YET)")
    .declare_key("physical_domains", Array(Integer(0)),
                "Numbers of physical domains (submeshes) that forms a segment and that "
                "will be removed from the computational mesh.")
    .declare_key("elements", Array(Integer(0)),
                        "Numbers of elements that forms a segment and that "
                        "will be removed from the computational mesh.")
    .declare_key("sides", Array( Array(Integer(0), 2,2) ),
                        "Pairs [ element, local_side] specifying sides that are part of the boundary segment."
                        "Sides are NOT removed from computation.");


Record Mesh::input_type
	= Record("Mesh","Record with mesh related data." )
	.declare_key("mesh_file", FileName::input(), Default::obligatory(),
			"Input file with mesh description.")
	.declare_key("regions", Array( RegionDB::region_input_type ), Default::optional(),
	        "List of additional region definitions not contained in the mesh.")
	.declare_key("sets", Array( RegionDB::region_set_input_type), Default::optional(),
	        "List of region set definitions.")
    .declare_key("boundary",
            Record( "BoundaryLists", "Lists of boundary regions and sets.\n "
                    "All regions with name starting with period '.' are automatically treated as boundary regions.")
            .declare_key("ids", Array( Integer(0) ), Default::optional(),
                    "List with ID numbers of boundary regions.")
            .declare_key("labels", Array( String() ), Default::optional(),
                    "List of labels of boundary regions.")
            .declare_key("sets", Array( String() ), Default::optional(),
                    "List of regions sets whose regions will be marked as boundary.")
            .close(),
            Default::optional(),
            "")
    .declare_key("neighbouring", FileName::input(), Default::obligatory(),
    		"File with mesh connectivity data.");



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

    side_nodes[0][0][0] = 0;
    side_nodes[0][1][0] = 1;


    side_nodes[1][0][0] = 0;
    side_nodes[1][0][1] = 1;

    side_nodes[1][1][0] = 1;
    side_nodes[1][1][1] = 2;

    side_nodes[1][2][0] = 2;
    side_nodes[1][2][1] = 0;


    side_nodes[2][0][0] = 1;
    side_nodes[2][0][1] = 2;
    side_nodes[2][0][2] = 3;

    side_nodes[2][1][0] = 0;
    side_nodes[2][1][1] = 2;
    side_nodes[2][1][2] = 3;

    side_nodes[2][2][0] = 0;
    side_nodes[2][2][1] = 1;
    side_nodes[2][2][2] = 3;

    side_nodes[2][3][0] = 0;
    side_nodes[2][3][1] = 1;
    side_nodes[2][3][2] = 2;
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

/**
 *  Setup whole topology for read mesh.
 */
void Mesh::setup_topology(istream *in) {
    F_ENTRY;
    Mesh *mesh=this;


    count_element_types();

    // topology
    make_neighbours_and_edges();
    element_to_neigh_vb();
    create_external_boundary();

    count_side_types();



    xprintf(MsgVerb, "Topology O.K.\n")/*orig verb 4*/;


    // cleanup
    {
        vector<Neighbour_both> empty_vec;
        neighbours_.swap(empty_vec);
    }

}



//
void Mesh::count_side_types()
{
    struct Side *sde;

    n_insides = 0;
    n_exsides = 0;
    FOR_SIDES(this,  sde )
        if (sde->is_external()) n_exsides++;
        else n_insides++;
}



void Mesh::make_neighbours_and_edges()
{
	// for each node we make a list of elements that use this node
	map<const Node*,vector<unsigned int> > node_elements;

	// create the node_elements map
	FOR_ELEMENTS( this, e )
		for (unsigned int n=0; n<e->n_nodes(); n++)
			node_elements[e->node[n]].push_back(e->index());

	// pointers to created edges
	vector<Edge *> tmp_edges;
	// Now we go through all element sides and create edges and neighbours
	FOR_ELEMENTS( this, e )
	{
		for (unsigned int s=0; s<e->n_sides(); s++)
		{
			// skip sides that were already found
			if (e->edges_[s] != NULL) continue;

			Neighbour neighbour;
			bool is_neighbour = false;
			unsigned int n_edg_sides = 0;

			// Find all elements that share this side.
			// element_count contains number of nodes that are used by the respective element.
			map<unsigned int,unsigned int> element_count;
			set<const Node *> set_of_nodes;
			for (unsigned n=0; n<e->side(s)->n_nodes(); n++)
			{
				set_of_nodes.insert(e->side(s)->node(n));
				for (vector<unsigned int>::iterator eit=node_elements[e->side(s)->node(n)].begin();
						eit!=node_elements[e->side(s)->node(n)].end(); eit++)
					element_count[*eit]++;
			}

			// Count the number of elements that share the whole side/edge.
			for (map<unsigned int, unsigned int>::iterator ec=element_count.begin(); ec!=element_count.end(); ec++)
				if (ec->second == e->side(s)->n_nodes())
				{
					if (element[ec->first].dim_ == e->dim_)
						n_edg_sides++;
					else if (element[ec->first].dim_ == e->dim_-1)
					{
						is_neighbour = true;
						neighbour.element_ = &(element[ec->first]);
						neighbour.sigma = 1;
					}
				}

			if (is_neighbour)
			{ // edge connects elements of different dimensions
				for (map<unsigned int, unsigned int>::iterator ec=element_count.begin(); ec!=element_count.end(); ec++)
					if (ec->second == e->side(s)->n_nodes())
					{
						Element *elem = &(element[ec->first]);
						if (elem->dim_ == e->dim_)
						{ // element with the same dimension: update edge
							// find local side index on element ec->first
							for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++)
							{
								SideIter si = elem->side(ecs);
								int ni=0;
								while (ni<si->n_nodes() && set_of_nodes.find(si->node(ni)) != set_of_nodes.end()) ni++;
								if (ni>=si->n_nodes())
								{
									// create a new edge and neighbour for this side
									Edge *edg = new Edge;
									tmp_edges.push_back(edg);
									edg->n_sides = 1;
									edg->side_ = new struct SideIter[1];
									edg->side_[0] = si;
									elem->edges_[ecs] = edg;
									neighbour.edge_ = edg;
									vb_neighbours_.push_back(neighbour);
									break;
								}
							}
						}
					}
			} else { // edge connects only elements of the same dimension
				// Allocate the array of sides.
				Edge *edg = new Edge;
				tmp_edges.push_back(edg);
				edg->n_sides = n_edg_sides;
				edg->side_ = new struct SideIter[edg->n_sides];
				unsigned int i=0;
				// initialize edge data
				for (map<unsigned int, unsigned int>::iterator ec=element_count.begin(); ec!=element_count.end(); ec++)
					if (ec->second == e->side(s)->n_nodes())
					{
						Element *elem = &(element[ec->first]);

						if (elem->dim_ != e->dim_) continue;

						// find local side index on element ec->first
						for (unsigned int ecs=0; ecs<elem->n_sides(); ecs++)
						{
							SideIter si = elem->side(ecs);
							int ni=0;
							while (ni<si->n_nodes() && set_of_nodes.find(si->node(ni)) != set_of_nodes.end()) ni++;
							if (ni>=si->n_nodes())
							{
								edg->side_[i++] = si;
								elem->edges_[ecs] = edg;
								break;
							}
						}
					}
			}
		}
	}

	// Now we can create the vector of edges, after we know its size
	edge.resize(tmp_edges.size());
	map<Edge*,Edge*> edge_map;
	// update pointers to edges and free the temporary objects
	for (unsigned int i=0; i<tmp_edges.size(); i++)
	{
		edge[i] = *(tmp_edges[i]);
		edge_map[tmp_edges[i]] = &(edge[i]);
		for (int s=0; s<edge[i].n_sides; s++)
			edge[i].side_[s]->element()->edges_[edge[i].side_[s]->el_idx()] = &(edge[i]);
		delete tmp_edges[i];
	}
	FOR_NEIGHBOURS(this, ngh)
		ngh->edge_ = edge_map[ngh->edge_];

	xprintf( Msg, "Created %d edges and %d neighbours.\n", edge.size(), vb_neighbours_.size() );
}





/**
 * Set Element->boundaries_ for all external sides. (temporary solution)
 */
void Mesh::create_external_boundary()
{
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



void Mesh::setup_materials( MaterialDatabase &base)
{
    xprintf( MsgVerb, "   Element to material... ")/*orig verb 5*/;
    FOR_ELEMENTS(this, ele ) {
        ele->material=base.find_id(ele->region().id());
        INPUT_CHECK( ele->material != base.end(),
                "Reference to undefined material %d in element %d\n", ele->region().id(), ele.id() );
    }
    xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;
}


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
     for(int i=0;i<n_elements(); ++i ) master_elements[i].reserve(sizes[i]);

     // fill intersec_elements
     for( vector<Intersection>::iterator i=intersections.begin(); i != intersections.end(); ++i )
     master_elements[i->master_iter().index()].push_back( i-intersections.begin() );

}

ElementAccessor<3> Mesh::element_accessor(unsigned int idx, bool boundary) {
    return ElementAccessor<3>(this, idx, boundary);
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
