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
 * @brief  Mesh construction
 *
 */

#include <unistd.h>

#include "mesh/ini_constants_mesh.hh"
#include "constantdb.h"

#include "system.hh"
#include "problem.h"
#include "mesh.h"
#include "hashes.h"
// think about following dependencies
#include "boundaries.h"
#include "initials.h"
#include "sources.h"
#include "concentrations.h"
#include "transport_bcd.h"
#include "transport.h"

//TODO: sources, concentrations, initial condition  and similarly boundary conditions should be
// instances of a Element valued field
// concentrations is in fact reimplemented in transport REMOVE it HERE

// After removing non-geometrical things from mesh, this should be part of mash initializing.
#include "topology.cc"
#include "msh_reader.h"
#include "msh_gmshreader.h"

void count_element_types(Mesh*);
void read_node_list(Mesh*);

Mesh::Mesh() {
    xprintf(Msg, " - Mesh()     - version with node_vector\n");

    n_materials = NDEF;
    concentration = NULL;
    l_concentration = NULL;
    transport_bcd = NULL;
    l_transport_bcd = NULL;
    //	n_sources        = NDEF;
    //	source           = NULL;
    //	l_source         = NULL;
    n_sides = NDEF;
    side = NULL;
    l_side = NULL;
    n_insides = NDEF;
    n_exsides = NDEF;
    n_edges = NDEF;
    edge = NULL;
    l_edge = NULL;
    n_neighs = NDEF;
    neighbour = NULL;
    l_neighbour = NULL;

    // Hashes
    n_lines = 0;
    n_triangles = 0;
    n_tetrahedras = 0;
    max_bou_id = NDEF;
    max_con_id = NDEF;
    max_tbc_id = NDEF;
    max_ngh_id = NDEF;

    concentration_hash = NULL;
    transport_bcd_hash = NULL;
    neighbour_hash = NULL;
}

//=============================================================================
// MAKE AND FILL ALL LISTS IN STRUCT MESH
//
// DF - method make_mesh() will be removed outside from this file
//=============================================================================

void make_mesh(struct Problem *problem) {
    F_ENTRY;

    ASSERT(!(problem == NULL), "NULL pointer as argument of function make_mesh()\n");

    const char* meshFileName = OptGetStr("Input", "Mesh", NULL);
    if (access(meshFileName, R_OK) != 0)
        xprintf(UsrErr, "Cannot read from file %s\n", meshFileName);

    Mesh* mesh = new Mesh();

    /* Test of object storage */
    ConstantDB::getInstance()->setObject(MESH::MAIN_INSTANCE, mesh);

    // get all file names
    // DF  - problem, it is not sure, why this is happening
    // ConstantDB::getInstance()->setChar("Mesh_geometry_fname", xstrcpy(meshFileName));

    // DF - Move to ConstantDB
    // mesh->material_fname = xstrcpy(problem->material_fname);
    // mesh->boundary_fname = xstrcpy(problem->boundary_fname);

    // if (ConstantDB::getInstance()->getInt("Problem_type") == UNSTEADY_SATURATED)
    // mesh->initial_fname = xstrcpy(problem->initial_fname);

    //if (OptGetBool("Transport", "Transport_on", "no") == true) {
    // mesh->concentration_fname = xstrcpy(problem->concentration_fname);
    // mesh->transport_bcd_fname = xstrcpy(problem->transport_bcd_fname);
    //}
    // mesh->neighbours_fname = xstrcpy(problem->neighbours_fname);
    // if (problem->sources_fname != NULL)
    //    mesh->sources_fname = xstrcpy(problem->sources_fname);
    //else
    //    mesh->sources_fname = NULL;


    // read all mesh files - this is work for MeshReader
    // read_node_list(mesh);
    // DF - elements are read by MeshReader
    // read_element_list(mesh);
    // --------------------- MeshReader testing - Begin
    MeshReader* meshReader = new GmshMeshReader();
    meshReader->read(OptGetStr("Input", "Mesh", NULL), mesh);
    // --------------------- MeshReader testing - End

    read_neighbour_list(mesh);
    //  if( mesh->sources_fname != NULL ) {
    //      read_source_list( mesh );
    //  }
    //  read_element_properties( mesh );
    //  make_element_geometry( mesh );

    make_side_list(mesh);
    make_edge_list(mesh);
    make_hashes(problem);
    count_element_types(mesh);

    // topology
    element_to_material(mesh, *(problem->material_database));
    node_to_element(mesh);
    element_to_side_both(mesh);
    neigh_vv_to_element(mesh);
    element_to_neigh_vv(mesh);
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


    read_boundary(mesh);

    if (OptGetBool("Transport", "Transport_on", "no") == true) {
        mesh->n_substances = problem->transport->n_substances;
        read_concentration_list(mesh);
        read_transport_bcd_list(mesh);
    }
    source_to_element_both(mesh);
    if (mesh->concentration != NULL) {
        concentration_to_element(mesh);
        transport_bcd_to_boundary(mesh);
    }



}

//=============================================================================
// COUNT ELEMENT TYPES
//=============================================================================

void count_element_types(Mesh* mesh) {
    //ElementIter elm;

    FOR_ELEMENTS(elm)
    switch (elm->type) {
        case 1:
            mesh->n_lines++;
            break;
        case 2:
            mesh->n_triangles++;
            break;
        case 4:
            mesh->n_tetrahedras++;
            break;
    }
}
//=============================================================================
// RETURN MAX NUMBER OF ENTRIES IN THE ROW
//=============================================================================

int *max_entry() {
    Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);

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

        size += elm->n_neighs_vv; // non-comp model

        max_size[0] += size;
        if (max_size[1] < size) max_size[1] = size;
    }
    // getchar();
    return max_size;
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
