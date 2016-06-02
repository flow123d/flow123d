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

#include <string>

#include "system/sys_profiler.hh"
#include "input/accessors.hh"

#include "io/output_data_base.hh"

class Mesh;

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
    
    /// Data vector.
    std::vector<T> data_;
};



class OutputMesh
{
public:
    OutputMesh(Mesh *mesh);
    ~OutputMesh();
    
    std::shared_ptr<MeshData<double>> nodes_;
    std::shared_ptr<MeshData<unsigned int>> connectivity_;
    std::shared_ptr<MeshData<unsigned int>> offsets_;
    
private:
    void fill_vectors(Mesh *mesh);
};