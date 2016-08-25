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
	this->set_vtk_type<T>();
    output_field_name = name;
    n_elem_ = n_elem;
}

template <class T>
MeshData<T>::~MeshData()
{};

/// Prints @p idx element of data vector into stream.
template <class T>
void MeshData<T>::print_ascii(std::ostream& out_stream, unsigned int idx)
{
    ASSERT_LT(idx, this->n_values);
    out_stream << data_[idx] ;
}

/// Prints the whole data vector into stream.
template <class T>
void MeshData<T>::print_ascii_all(std::ostream& out_stream)
{
    for(auto &d : data_)
        out_stream << d << " ";
}


/// Prints the whole data vector into stream.
template <class T>
void MeshData<T>::print_binary_all(ostream &append_stream, ostream &out_stream)
{
    print_binary_inner(append_stream, out_stream);
}


/// Prints the whole data vector into stream.
template <class T>
void MeshData<T>::print_all_yaml(std::ostream& out_stream, unsigned int precision)
{
    ASSERT(false).error("Unsupported output of the mesh data to YAML format.");
}


template<class T>
template<class Q>
void MeshData<T>::print_binary_inner(ostream &append_stream, ostream &out_stream,
		typename std::enable_if<std::is_floating_point<Q>::value >::type* )
{
	// variables storing range of data array
	Q range_min = std::numeric_limits<Q>::max();
	Q range_max = std::numeric_limits<Q>::min();
	// write size of data
	unsigned long long int data_byte_size = data_.size() * sizeof(Q);
	append_stream.write(reinterpret_cast<const char*>(&data_byte_size), sizeof(unsigned long long int));
	// write data
    for(auto &d : data_) {
    	append_stream.write(reinterpret_cast<const char*>(&d), sizeof(Q));
       	if (d < range_min) range_min = d;
       	if (d > range_max) range_max = d;
    }
    out_stream << "RangeMin=\"" << range_min << "\" RangeMax=\"" << range_max << "\" ";
}


template<class T>
template<class Q>
void MeshData<T>::print_binary_inner(ostream &append_stream, ostream &out_stream,
		typename std::enable_if<!std::is_floating_point<Q>::value >::type* )
{
	// write size of data
	unsigned long long int data_byte_size = data_.size() * sizeof(Q);
	append_stream.write(reinterpret_cast<const char*>(&data_byte_size), sizeof(unsigned long long int));
	// write data
    for(auto &d : data_)
    	append_stream.write(reinterpret_cast<const char*>(&d), sizeof(Q));
    out_stream << "RangeMin=\"\" RangeMax=\"\" ";
}


/// Access i-th element in the data vector.
template <class T>
T& MeshData<T>::operator[](unsigned int i)
{
    ASSERT_DBG(i < data_.size());
    return data_[i];
}

template class MeshData<unsigned int>;
template class MeshData<double>;

