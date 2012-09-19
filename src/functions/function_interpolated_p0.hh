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


template <int dim>
class FunctionInterpolatedP0: public FunctionBase<dim> {
public:
	typedef arma::vec::fixed<dim> Point;

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
	 */
	void set_source_of_interpolation(const std::string & mesh_file,
			const std::string & raw_output);

    /**
     * Returns one scalar value in one given point.
     */
	double value(const Point &p, const unsigned int component);

    /**
     * Returns one vector value in one given point.
     */
    void vector_value(const Point &p, std::vector<double> &value);

    /**
	 * Returns std::vector of scalar values in several points at once.
	 */
    void value_list (const std::vector<Point>  &point_list,
				  std::vector<double>         &value_list,
				  const unsigned int  component);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    void vector_value_list (const std::vector<Point> &point_list,
                            std::vector< std::vector<double> > &value_list);

protected:
	/// element for interpolation
	ElementFullIter element_;

	/// mesh
	Mesh* mesh_;

};

#endif /* FUNCTION_INTERPOLATED_P0_HH_ */
