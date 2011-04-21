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
 * @brief  Classes that model discrete mathematical fields.
 */

#ifndef FIELD_P0_HH_
#define FIELD_P0_HH_

#include "system/system.hh"
#include "xio.h"
#include <string>
#include "mesh/mesh_types.hh"
#include "mesh/mesh.h"


/**
 * Class can read scalar field - piecewise constant on elements.
 * Since we want to template it by dimension of returned vector or tensor, we keep implementation in the header file.
 *
 * TODO:
 * - use PETSC vectors for storing data (paralelization), field as accessor to a PETSC vector
 *   problem jak to pouzit v transportu: potrebuju mit pristup po jednotlivych letkach i po jednotlivych elementech
 *
 * - parallelization by given partitioning of mesh
 *   should provide: distribution, iterator over local elements, access to mesh elements,
 *   access through old mesh element indexes
 *
 * - P0 fields given by values for materials
 * - need output through streams - need something like xopen for streams distinguish output and input
 *   until that we can introduce type cast to string for : double, vector, tensor - for output
 *   - and vice versa for input
 * - fields for RT0 and Edge P0
 */

template <class Value>
class FieldP0 {
public:

    FieldP0(Mesh *mesh);
    void read_field(const string f_name, const string section);
    inline Value element_value(const unsigned int el_idx) const
        { return data[el_idx]; }
    void integrate_material_subdomains(vector<Value> &integrals);

    ~FieldP0() {};

private:
    std::vector<Value> data;
    Mesh *mesh;

};

// template implementation --------------------------------------------------

template <class Value>
FieldP0<Value>::FieldP0(Mesh *mesh)
: mesh(mesh), data(0)
{}


template <class Value>
void FieldP0<Value>::read_field(const string f_name,const string section)
{
    FILE    *in;   // input file
    char     line[ LINE_SIZE ];   // line of data file
    std::istringstream line_stream;
    ElementFullIter elm(mesh->element);
    int eid,n_sources;
    double value;

    F_ENTRY;

    in = xfopen( f_name.c_str(), "rt" );
    INPUT_CHECK( skip_to( in, section.c_str() ) , "Can not find $Sources section in the file %s.\n", f_name.c_str());

    xfgets( line, LINE_SIZE - 2, in );
    line_stream.str(string(line));
    line_stream >> n_sources;
    INPUT_CHECK( line_stream.good(), "Can not convert input to int.\n");

    data.resize(mesh->element.size());

    for(int i_line=0; i_line<n_sources; i_line++) {
        xfgets( line, LINE_SIZE - 2, in );
        line_stream.str(string(line));
        line_stream >> eid;
        INPUT_CHECK( line_stream.good(), "Can not convert input to int.\n");
        // get index from ID
        elm=mesh->element.find_id(eid);
        INPUT_CHECK( elm != mesh->element.end(), "Wrong element id on line: %d file: %s\n", i_line, f_name.c_str());
        // make value from rest of line
        line_stream >> data[elm.index()];
    }
    xfclose( in );
}

#endif /* FIELD_P0_HH_ */
