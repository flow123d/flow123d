/*
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a refer to Flow123d on your project site if you use the program for any purpose.
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of
 * the GNU General Public License version 3 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program;
 * if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 */

/**
 * $Id: field_p0.hh 807 2010-12-13 19:29:27Z jan.brezina $
 * $Revision: 807 $  $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2010-12-13 20:29:27 +0100 (Mon, 13 Dec 2010) $
 *
 * @file
 * @brief  Simple mathematical field piecewise constant on the elements of the mesh.
 */

#ifndef FIELD_P0_HH_
#define FIELD_P0_HH_

#include "system.hh"
#include "xio.h"
#include <string>
#include "mesh_types.hh"
#include "mesh.h"


/**
 * Class can read scalar field - piecewise constant on elements.
 * Since we want to template it by dimension of returned vector or tensor, we keep implementation in the header file.
 */

class FieldP0 {
public:
    FieldP0(const string f_name, const string section, ElementVector &el_vec);
    inline double element_value(const unsigned int el_idx) { return data[el_idx]; }

    ~FieldP0() {};

private:
    std::vector<double> data;

};


FieldP0::FieldP0(const string f_name,const string section, ElementVector &el_vec)
: data(el_vec.size(),0)
{
    FILE    *in;   // input file
    char     line[ LINE_SIZE ];   // line of data file
    ElementFullIter elm(el_vec);
    int eid,n_sources;
    double value;

    F_ENTRY;

    in = xfopen( f_name.c_str(), "rt" );
    INPUT_CHECK( skip_to( in, section.c_str() ) , "Can not find $Sources section in the file %s.\n", f_name.c_str());

    xfgets( line, LINE_SIZE - 2, in );
    n_sources = atoi( xstrtok( line) );

    for(int i_line=0; i_line<n_sources; i_line++) {
        xfgets( line, LINE_SIZE - 2, in );
        eid = atoi( xstrtok( line) );
        value = atof( xstrtok(NULL) );

        elm=el_vec.find_id(eid);
        INPUT_CHECK( elm != el_vec.end(), "Wrong element id on line: %d file: %s\n", i_line, f_name.c_str());

        data[elm.index()] = value;
    }
    xfclose( in );
}

#endif /* FIELD_P0_HH_ */
