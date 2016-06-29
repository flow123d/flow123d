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

#include "io/output_data_base.hh"

#include "fields/field_values.hh"

#include "tools/general_iterator.hh"

class Mesh;
template<int, class Value> class Field;

class OutputElement;
typedef GeneralIterator<OutputElement> OutputElementIterator;

/// Class representing data vector of geometry and topology information (especially for VTK).
/// Filling the vector is the users responsibility.
template <typename T>
class MeshData : public OutputDataBase {
public:
    /// Constructor. @p name is the possible name of the output vector.
    MeshData(std::string name, NumCompValueType n_elem = N_SCALAR)
    {
        output_field_name = name;
        n_elem_ = n_elem;
    }
    
    ~MeshData() override 
    {};
    
    /// Prints @p idx element of data vector into stream.
    void print(std::ostream& out_stream, unsigned int idx) override {
        ASSERT_LE(idx, this->n_values);
        out_stream << data_[idx] ;
    }
    
    /// Prints the whole data vector into stream.
    void print_all(std::ostream& out_stream) override {
        for(auto &d : data_)
            out_stream << d << " ";
    }
    
    /// Access i-th element in the data vector.
    T& operator[](unsigned int i){
        ASSERT(i < data_.size());
        return data_[i];
    }
    
    /// Data vector.
    std::vector<T> data_;
};


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
    /// Shortcut instead of spacedim template. We suppose only spacedim=3 at the moment. 
    static const unsigned int spacedim = 3;
    
    /// Constructor. Takes computational mesh as a parameter.
    OutputMeshBase(Mesh* mesh);
    virtual ~OutputMeshBase();
    
    /// Gives iterator to the FIRST element of the output mesh.
    OutputElementIterator begin();
    /// Gives iterator to the LAST element of the output mesh.
    OutputElementIterator end();
    
    /// Vector of element indices in the computational mesh. (Important when refining.)
    std::shared_ptr<std::vector<unsigned int>> orig_element_indices_;
    
    /// Vector of node coordinates. [spacedim x n_nodes]
    std::shared_ptr<MeshData<double>> nodes_;
    /// Vector maps the nodes to their coordinates in vector @p nodes_.
    std::shared_ptr<MeshData<unsigned int>> connectivity_;
    /// Vector of offsets of node indices of elements. Maps elements to their nodes in connectivity_.
    std::shared_ptr<MeshData<unsigned int>> offsets_;
    
    /// Returns number of nodes.
    unsigned int n_nodes();
    /// Returns number of element.
    unsigned int n_elements();
    
protected:
    /// Pointer to the computational mesh.
    Mesh *orig_mesh_;
    
    /// Friend provides access to vectors for element accessor class.
    friend class OutputElement;
};


/// @brief Class represents output mesh with continuous elements.
class OutputMesh : public OutputMeshBase
{
public:
    OutputMesh(Mesh* mesh);
    ~OutputMesh();
    
    /// Creates the output mesh identical to the computational one.
    void create_identical_mesh();
    
    /// Creates refined mesh.
    void create_refined_mesh(Field<3, FieldValue<3>::Scalar> *error_control_field);
    
protected:
    bool refinement_criterion();
    
    const unsigned int max_level = 2;
    
    /// Friend provides access to vectors for discontinous output mesh.
    friend class OutputMeshDiscontinuous;
};


class OutputMeshDiscontinuous : public OutputMeshBase
{
public:
    OutputMeshDiscontinuous(Mesh* mesh);
    ~OutputMeshDiscontinuous();
    
    /// Creates output mesh from the given continuous one.
    void create_mesh(std::shared_ptr<OutputMesh> output_mesh);
    
    /// Creates discontinuous refined mesh.
    void create_refined_mesh(Field<3, FieldValue<3>::Scalar> *error_control_field);
    
protected:
    bool refinement_criterion();
    
    const unsigned int max_level = 2;
};

#endif  // OUTPUT_MESH_HH_

