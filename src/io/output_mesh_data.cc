/*
 * output_mesh_data.cc
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#include "io/output_mesh_data.hh"


/// Constructor. @p name is the possible name of the output vector.
template <class T>
MeshData<T>::MeshData(std::string name, NumCompValueType n_elem)
{
    output_field_name = name;
    n_elem_ = n_elem;
}

template <class T>
MeshData<T>::~MeshData()
{};

/// Prints @p idx element of data vector into stream.
template <class T>
void MeshData<T>::print(std::ostream& out_stream, unsigned int idx)
{
    ASSERT_LE(idx, this->n_values);
    out_stream << data_[idx] ;
}

/// Prints the whole data vector into stream.
template <class T>
void MeshData<T>::print_all(std::ostream& out_stream)
{
    for(auto &d : data_)
        out_stream << d << " ";
}


/// Prints the whole data vector into stream.
template <class T>
void MeshData<T>::print_all_yaml(std::ostream& out_stream, unsigned int precision)
{
    ASSERT(false).error("Unsupported output of the mesh data to YAML format.");
}


/// Access i-th element in the data vector.
template <class T>
T& MeshData<T>::operator[](unsigned int i)
{
    ASSERT(i < data_.size());
    return data_[i];
}

template class MeshData<unsigned int>;
template class MeshData<double>;

