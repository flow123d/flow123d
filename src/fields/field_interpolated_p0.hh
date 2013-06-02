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


#ifndef FIELD_INTERPOLATED_P0_HH_
#define FIELD_INTERPOLATED_P0_HH_

#include "field_base.hh"
#include "mesh/mesh.h"
#include "mesh/mesh_types.hh"
#include "system/system.hh"
#include "mesh/msh_gmshreader.h"
#include "mesh/bih_tree.hh"
#include "mesh/ngh/include/abscissa.h"
#include "mesh/ngh/include/triangle.h"
#include "mesh/ngh/include/tetrahedron.h"


template <int spacedim, class Value>
class FieldInterpolatedP0: public FieldBase<spacedim, Value> {
public:

	/**
	 * Constructor
	 */
	FieldInterpolatedP0(unsigned int n_comp=0);


	/**
	 * Declare Input type.
	 */
	static Input::Type::Record input_type;

	static Input::Type::Record get_input_type(Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit);

	/**
	 * Initialization from the input interface.
	 */
	virtual void init_from_input(const Input::Record &rec);


    /**
     * Update time and possibly update data from GMSH file.
     */
    virtual bool set_time(double time);


    /**
	 * Set sources files of interpolation
	 *
	 * @param mesh_file file contained data of mesh
	 * @param raw_output file contained output
	 *
	 * TODO: use streams instead of filenames (better testing)
	 * TODO: use just one GMSH file for both mesh and data (consistency)
	 */
	//void set_source_of_interpolation(const FilePath & mesh_file,
	//								 const FilePath & raw_output);


    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list(const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);

protected:
	/// mesh, which is interpolated
	Mesh* source_mesh_;

	/// mesh reader
	GmshMeshReader *reader_;

    /// Raw buffer of n_entities rows each containing Value::size() doubles.
    double *data_;

	/// vector stored suspect elements in calculating the intersection
	std::vector<unsigned int> searched_elements_;

	/// field name read from input
	std::string field_name_;

	/// tree of mesh elements
	BIHTree* bih_tree_;

	/// stored reference to last computed element
	const ElementAccessor<spacedim> *computed_elm_;

	/**
	 * Read scalar element data with name @p field_name using tokenizer @p tok initialized
	 * over a GMSH file or stream.
	 *
	 * This needs lot of work to be general enough to be outside of this class
	 * TODO:
	 * - we need concept of Fields so that we can fill corresponding vectors
	 * - then we should support scalar as well as vector or even tensor data
	 * - support for time dependent data
	 * - selective reading on submesh (parallelism - subdomains, or boundary data)
	 *
	 */
	//void read_element_data_from_gmsh(Tokenizer &tok, const  string &field_name);


	/**
	 * Calculate pressures in triangle element
	 */
	void calculate_triangle_value(TTriangle &element, unsigned int idx);
	/**
	 * Calculate pressures in abscissa element
	 */
	void calculate_abscissa_value(TAbscissa &element, unsigned int idx);

public:
	/**
	 * Create tetrahedron from element
	 */
	static void createTetrahedron(Element *ele, TTetrahedron &te);

	/**
	 * Create triangle from element
	 */
	static void createTriangle(Element *ele, TTriangle &tr);

	/**
	 * Create abscissa from element
	 */
	static void createAbscissa(Element *ele, TAbscissa &ab);
};



#endif /* FUNCTION_INTERPOLATED_P0_HH_ */
