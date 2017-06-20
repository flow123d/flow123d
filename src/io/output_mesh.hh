/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    output_mesh.hh
 * @brief   Classes for auxiliary output mesh.
 */

#ifndef OUTPUT_MESH_HH_
#define OUTPUT_MESH_HH_

#include <string>

#include "system/sys_profiler.hh"
#include "input/accessors.hh"

#include "fields/field_values.hh"
#include "fields/field_set.hh"

#include "tools/general_iterator.hh"

class Mesh;
template<int, class Value> class Field;
template<class T> class ElementDataCache;

class OutputElement;
typedef GeneralIterator<OutputElement> OutputElementIterator;

class OutputMeshBase;
class OutputMesh;
class OutputMeshDiscontinuous;


/// @brief Base class for Output mesh.
/**
 * Defines common members for OutputMesh classes.
 */
class OutputMeshBase : public std::enable_shared_from_this<OutputMeshBase>
{
public:
    DECLARE_EXCEPTION(ExcFieldNotScalar, << "Field '" << FieldCommon::EI_Field::qval
                                         << "' is not scalar in spacedim 3.");
    
    /// Shortcut instead of spacedim template. We suppose only spacedim=3 at the moment. 
    static const unsigned int spacedim = 3;
    
    /// Constructor. Takes computational mesh as a parameter.
    OutputMeshBase(Mesh &mesh);
    /// Constructor. Takes computational mesh and input record as a parameters.
    OutputMeshBase(Mesh &mesh, const Input::Record &in_rec);
    virtual ~OutputMeshBase();
    
    /**
     * @brief The specification of output mesh.
     * @return record for output mesh
     */
    static const Input::Type::Record & get_input_type();
    
    /// Gives iterator to the FIRST element of the output mesh.
    OutputElementIterator begin();
    /// Gives iterator to the LAST element of the output mesh.
    OutputElementIterator end();
    
    /// Selects the error control field out of output field set according to input record.
    void select_error_control_field(FieldSet &output_fields);
    
    /// Vector of element indices in the computational mesh. (Important when refining.)
    std::shared_ptr<std::vector<unsigned int>> orig_element_indices_;
    
    /// Vector of node coordinates. [spacedim x n_nodes]
    std::shared_ptr<ElementDataCache<double>> nodes_;
    /// Vector maps the nodes to their coordinates in vector @p nodes_.
    std::shared_ptr<ElementDataCache<unsigned int>> connectivity_;
    /// Vector of offsets of node indices of elements. Maps elements to their nodes in connectivity_.
    std::shared_ptr<ElementDataCache<unsigned int>> offsets_;
    
    /// Returns number of nodes.
    unsigned int n_nodes();
    /// Returns number of element.
    unsigned int n_elements();
    
protected:
    /// Input record for output mesh.
    Input::Record input_record_;
    
    /// Pointer to the computational mesh.
    Mesh *orig_mesh_;
    
    /// Maximal level of refinement.
    const unsigned int max_level_;
    
    /// Refinement error control field.
    Field<3, FieldValue<3>::Scalar> *error_control_field_;
    
    /// Friend provides access to vectors for element accessor class.
    friend class OutputElement;
};


/// @brief Class represents output mesh with continuous elements.
class OutputMesh : public OutputMeshBase
{
public:
    OutputMesh(Mesh &mesh);
    OutputMesh(Mesh &mesh, const Input::Record &in_rec);
    ~OutputMesh();
    
    /// Creates the output mesh identical to the computational one.
    void create_identical_mesh();
    
    /// Creates refined mesh.
    void create_refined_mesh();
    
protected:
    bool refinement_criterion();
    
    /// Friend provides access to vectors for discontinous output mesh.
    friend class OutputMeshDiscontinuous;
};


class OutputMeshDiscontinuous : public OutputMeshBase
{
public:
    OutputMeshDiscontinuous(Mesh &mesh);
    OutputMeshDiscontinuous(Mesh &mesh, const Input::Record& in_rec);
    ~OutputMeshDiscontinuous();
    
    /// Creates output mesh from the given continuous one.
    void create_mesh(std::shared_ptr<OutputMesh> output_mesh);
    
    /// Creates discontinuous refined mesh.
    void create_refined_mesh();
    
protected:
    bool refinement_criterion();
};

#endif  // OUTPUT_MESH_HH_

