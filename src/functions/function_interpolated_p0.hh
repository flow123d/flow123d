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
 * $Id: function_interpolated_p0.hh 1567 2012-02-28 13:24:58Z jan.brezina $
 * $Revision: 1567 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
 *
 *
 */

#ifndef FUNCTION_INTERPOLATED_P0_HH_
#define FUNCTION_INTERPOLATED_P0_HH_

#include "function_base.hh"
#include "mesh/mesh.h"
#include "mesh/mesh_types.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
#include "new_mesh/bounding_interval_hierarchy.hh"
#include "new_mesh/ngh/include/triangle.h"
#include "new_mesh/ngh/include/tetrahedron.h"

//enum TIntersectionType;

class FunctionInterpolatedP0: public FunctionBase<3> {
public:

	/**
	 * Constructor
	 */
	FunctionInterpolatedP0();

    /**
     * Set element for interpolation
     */
	void set_element(ElementFullIter &element);

	/**
	 * Set sources files of interpolation
	 *
	 * @param mesh_file file contained data of mesh
	 * @param raw_output file contained output
	 * @param ngh_file file with specification of adjacency between dimensions
	 * @param bcd_file file with boundary conditions
	 */
	void set_source_of_interpolation(const std::string & mesh_file,
									 const std::string & raw_output,
									 const std::string & ngh_file,
									 const std::string & bcd_file);

    /**
     * Returns one scalar value in one given point.
     */
    virtual double value(const Point &p, const unsigned int  component = 0) const;
    /**
     * Returns one vector value in one given point.
     */
    virtual void   vector_value(const Point &p, std::vector<double>     &value) const;

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void   value_list (const std::vector< Point >  &point_list,
                       std::vector<double>         &value_list,
                       const unsigned int  component = 0) const;

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void   vector_value_list (const std::vector< Point >    &point_list,
                              std::vector< std::vector<double> >      &value_list) const;

protected:
	/// element for interpolation
	ElementIter element_;

	/// mesh
	Mesh* mesh_;

	/// vector of pressures in nodes
	std::vector<double> pressures_;

	/// tree of mesh elements
	BoundingIntevalHierachy* bihTree_;

	/**
	 * Read pressures from file and put them to vector pressures_
	 *
	 * @param raw_output file contained output
	 */
	void read_pressures(FILE* raw_output);

	/// calculate pressures in inner mesh
	void calculate_interpolation();

	/**
	 * Calculate pressures in element
	 */
	double calculate_element(TTriangle &element, std::vector<int> &searchedElements);

	/**
	 * Create tetrahedron from element
	 */
	void createTetrahedron(ElementFullIter ele, TTetrahedron &te);

};



#endif /* FUNCTION_INTERPOLATED_P0_HH_ */
