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
 * @brief ???
 *
 */

#ifndef INITIALS_H
#define INITIALS_H

class Mesh;

//=============================================================================
// STRUCTURE OF THE INITIAL CONDITION
//=============================================================================
struct Initial
{
    // Data readed from initial conditions files
    int      id;        // Id number of condition
    int      eid;           // Id number of element where prescribed
    double   epress;    // Press in the element's center
    int      n_sides;       // # of sides of the element readed from file
    double   *spress;       // Presses on the sides of the element
    double   *sflux;    // Fluxes via sides of the element
    int      n_neighs_vv;   // # of neighbours, V-V type (noncomp.)
    int  *hid;      // ID numbers of neighbour
    // Topology of the mesh
    // List
    struct Initial *prev;   // Previous initial in the list
    struct Initial *next;   // Next initial in the list
    // Misc
    int  aux;       // Auxiliary flag
    double   faux;      // Auxiliary number
};

#define FOR_INITIALS(i)   for((i)=mesh->initial;(i)!=NULL;(i)=(i)->next)

void read_initial_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
