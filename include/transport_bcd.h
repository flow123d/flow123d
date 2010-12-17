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

#ifndef TRANSPORT_BCD_H
#define TRANSPORT_BCD_H

class Mesh;

//=============================================================================
// STRUCTURE OF THE TRANSPORT_BCD
//=============================================================================
struct Transport_bcd {
    // Data readed from transport bounary conditions files
    int id; // Id number of condition
    int bid; // ID number of boundary condition where prescribed
    int n_subst; // # of substances
    double *conc; // Values of concentrations

    // List
    struct Transport_bcd *prev; // Previous transport boundary in the list
    struct Transport_bcd *next; // Next transport boundary in the list

    // Misc
    int aux; // Auxiliary flag
    double faux; // Auxiliary number
};

#define FOR_TRANSPORT_BCDS(i)   for((i)=mesh->transport_bcd;(i)!=NULL;(i)=(i)->next)

void read_transport_bcd_list(Mesh*);

#endif
//-----------------------------------------------------------------------------
// vim: set cindent:
