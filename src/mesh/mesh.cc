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


#include "system/system.hh"
#include "input/input_type.hh"

#include "mesh/mesh.h"

// think about following dependencies
#include "mesh/boundaries.h"


//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "mesh/topology.cc"
#include "mesh/msh_reader.h"
#include "mesh/msh_gmshreader.h"

void count_element_types(Mesh*);
void read_node_list(Mesh*);


Input::Type::Record BoundarySegment::get_input_type() {
    using namespace Input::Type;

    static Record rec("BoundarySegment","Record with specification of boundary segments,\n"
            "i.e. subsets of the domain boundary where we prescribe one boundary condition.\n"
            "Currently very GMSH oriented. (NOT IMPLEMENTED YET)");
    if (!rec.is_finished()) {
        rec.declare_key("physical_domains", Array(Integer(0)),
                "Numbers of physical domains (submeshes) that forms a segment and that "
                "will be removed from the computational mesh.");
        rec.declare_key("elements", Array(Integer(0)),
                        "Numbers of elements that forms a segment and that "
                        "will be removed from the computational mesh.");
        rec.declare_key("sides", Array( Array(Integer(0), 2,2) ),
                        "Pairs [ element, local_side] specifying sides that are part of the boundary segment."
                        "Sides are NOT removed from computation.");
        rec.finish();
    }
    return rec;
}

Mesh::Mesh(Input::Record in_record)
: in_record_(in_record) {

    n_materials = NDEF;
    n_sides = NDEF;
    side = NULL;
    l_side = NULL;
    n_insides = NDEF;
    n_exsides = NDEF;
    n_neighs = NDEF;
    neighbour = NULL;
    l_neighbour = NULL;

    // Hashes
    n_lines = 0;
    n_triangles = 0;
    n_tetrahedras = 0;



}


Input::Type::Record Mesh::get_input_type() {
    using namespace Input::Type;

    static Record rec("Mesh","Record with mesh related data." );
    if (!rec.is_finished()) {
        rec.declare_key("mesh_file", FileName::input(), Default::obligatory(),
                "Input file with mesh description.");
        rec.declare_key("boundary_segmants", Array( BoundarySegment::get_input_type() ),
                "Array with specification of boundary segments");
        rec.declare_key("neighbouring", FileName::input(), Default::obligatory(),
                "File with mesh connectivity data.");
        rec.finish();
    }
    return rec;
}


//=============================================================================
// COUNT ELEMENT TYPES
//=============================================================================

void Mesh::count_element_types() {
    Mesh *mesh = this;

    FOR_ELEMENTS(this, elm)
    switch (elm->type) {
        case 1:
            n_lines++;
            break;
        case 2:
            n_triangles++;
            break;
        case 4:
            n_tetrahedras++;
            break;
    }
}
//=============================================================================
// RETURN MAX NUMBER OF ENTRIES IN THE ROW
//=============================================================================
/*
int *max_entry() {

    int *max_size, size, i;
//    ElementIter elm;

    max_size = (int*) xmalloc(2 * sizeof (int));

    max_size[0] = 0; // entries count
    max_size[1] = 0; // row

    FOR_ELEMENTS(elm) {
        size = 0; // uloha se zapornymi zdroji =1

        FOR_ELEMENT_SIDES(elm, si) { //same dim
            if (elm->side[si]->cond != NULL) size++;
            else
                size += elm->side[ si ]->edge->n_sides - 1;
            if (elm->side[si]->neigh_bv != NULL) size++; // comp model
        } // end same dim


        //   printf("SD id:%d,size:%d\n",elm->id,size);
        //

        FOR_ELM_NEIGHS_VB(elm, i) { // comp model
            size += elm->neigh_vb[i]->n_elements - 1;
            //     printf("VB id:%d,size:%d\n",elm->id,size);
        } // end comp model


        /*
        if (elm->dim > 1)
          FOR_NEIGHBOURS(ngh)
            FOR_NEIGH_ELEMENTS(ngh,n)
              if (ngh->element[n]->id == elm->id && n == 1)
                size++;
         */
/*
        size += elm->n_neighs_vv; // non-comp model

        max_size[0] += size;
        if (max_size[1] < size) max_size[1] = size;
    }
    // getchar();
    return max_size;
}*/

void Mesh::setup_topology() {
    Mesh *mesh=this;

    /// initialize mesh topology (should be handled inside mesh object)
    read_neighbour_list(mesh, in_record_.val<FilePath>("neighbouring") );

    make_side_list( mesh);
    make_edge_list(mesh);

    //    make_hashes(problem);
    count_element_types();

    // topology
    node_to_element(mesh);
    element_to_side_both(mesh);
    neigh_vv_to_element(mesh);
    //element_to_neigh_vv(mesh);
    neigh_vb_to_element_and_side(mesh);
    neigh_bv_to_side(mesh);
    element_to_neigh_vb(mesh);
    side_shape_specific(mesh);
    side_to_node(mesh);
    neigh_bb_topology(mesh);
    neigh_bb_to_edge_both(mesh);
    edge_to_side_both(mesh);
    neigh_vb_to_edge_both(mesh);
    side_types(mesh);
    count_side_types(mesh);
    xprintf(MsgVerb, "Topology O.K.\n")/*orig verb 4*/;



}

void Mesh::setup_materials( MaterialDatabase &base)
{
    xprintf( MsgVerb, "   Element to material... ")/*orig verb 5*/;
    FOR_ELEMENTS(this, ele ) {
        ele->material=base.find_id(ele->mid);
        INPUT_CHECK( ele->material != base.end(),
                "Reference to undefined material %d in element %d\n", ele->mid, ele.id() );
    }
    xprintf( MsgVerb, "O.K.\n")/*orig verb 6*/;

    make_element_geometry();
}

/**
 * CALCULATE PROPERTIES OF ALL ELEMENTS OF THE MESH
 */
void Mesh::make_element_geometry() {

    xprintf(Msg, "Calculating properties of elements... ")/*orig verb 2*/;

    ASSERT(element.size() > 0, "Empty mesh.\n");

    FOR_ELEMENTS(this, ele) {
        //DBGMSG("\n ele: %d \n",ele.id());
        //FOR_ELEMENTS(ele1) {
        //    printf("%d(%d) ",ele1.id(),ele1->type);
        //    ele1->bas_alfa[0]=1.0;
       // }
        ele->calc_metrics();
        ele->calc_volume();
        ele->calc_centre();
    }

    xprintf(Msg, "O.K.\n")/*orig verb 2*/;
}

//=============================================================================
// ID-POS TRANSLATOR
//=============================================================================
// TODO: should be method of water_linsys
/*
int id2pos(Mesh* mesh, int id, int* list, int type) {
    int i, limit;

    switch (type) {
        case ELM:
            limit = mesh->n_elements();
            break;
        case BC:
            limit = mesh->n_boundaries();
            break;
        case NODE:
            limit = mesh->node_vector.size();
            break;
    }

    for (i = 0; i < limit; i++)
        if (list[i] == id)
            break;

    if (type != BC)
        return i;
    else
        return i + mesh->n_elements();
}*/
/*
  for(i=0;i< ((type == ELM) ? mesh->n_elements() : mesh->n_boundaries );i++)
        if(list[i] == id)
                return (type == ELM) ? i : (i + mesh->n_elements());
 */
//=============================================================================
// MAKE ID-POS ELEMENT & BOUNDARY LIST
//=============================================================================
// TODO: should be private method of water_linsys
/*
void make_id2pos_list() {
    F_ENTRY;

    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

    ElementIter elm;
    NodeIter node;
    int i, j, s;

    mesh->epos_id = (int*) xmalloc(mesh->n_elements() * sizeof (int));
    mesh->spos_id = (int*) xmalloc(mesh->n_boundaries() * sizeof (int));
    mesh->npos_id = (int*) xmalloc(mesh->node_vector.size() * sizeof (int));

    j = i = 0;

    FOR_ELEMENTS(elm) {
        mesh->epos_id[i++] = elm.id();
        FOR_ELEMENT_SIDES(elm, s)
        if (elm->side[s]->cond != NULL)
            mesh->spos_id[j++] = elm->side[s]->id;
    }
    i = 0;

    FOR_NODES( node ) {
        mesh->npos_id[i++] = node->id;
    }
}
*/
//-----------------------------------------------------------------------------
// vim: set cindent:
