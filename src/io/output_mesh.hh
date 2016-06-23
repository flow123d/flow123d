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
#include "mesh/point.hh"


class Mesh;
template<int, class Value> class Field;
class OutputElementIterator;

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
    
    T& operator[](unsigned int i){
        ASSERT(i < data_.size());
        return data_[i];
    }
    
    /// Data vector.
    std::vector<T> data_;
};


class OutputMesh
{
public:
    static const unsigned int spacedim = 3;
    
    OutputMesh(Mesh* mesh, unsigned int max_refinement_level = 2);
    ~OutputMesh();
    
    void create_identical_mesh();
    void create_refined_mesh(Field<3, FieldValue<3>::Scalar> *error_control_field);
    
    OutputElementIterator begin();
    OutputElementIterator end();
    
    std::shared_ptr<std::vector<unsigned int>> orig_element_indices_;
    std::shared_ptr<std::vector<double>> local_nodes_;
    
    std::shared_ptr<MeshData<double>> nodes_;
    std::shared_ptr<MeshData<unsigned int>> connectivity_;
    std::shared_ptr<MeshData<unsigned int>> offsets_;
    
    void compute_discontinuous_data();
    std::shared_ptr<MeshData<double>> discont_nodes_;
    std::shared_ptr<MeshData<unsigned int>> discont_connectivity_;
    
    unsigned int n_nodes();
    unsigned int n_nodes_disc();
    unsigned int n_elements();
    
protected:
    
    struct AuxElement{
        std::vector<Space<spacedim>::Point> nodes;
        std::vector<unsigned int> connectivity;
        unsigned int level;
    };
    
    void fill_vectors();
    bool refinement_criterion(const AuxElement& ele);
    
    template<int dim>
    void refine_aux_element(const OutputMesh::AuxElement& aux_element,
                            std::vector< OutputMesh::AuxElement >& refinement,
                            unsigned int& last_node_idx);
    
    Mesh *orig_mesh_;
    bool discont_data_computed_;
    
    friend class OutputElement;
    
    std::vector<std::vector<AuxElement>> unit_refinement_;
    const unsigned int max_refinement_level_;
};


#endif  // OUTPUT_MESH_HH_

