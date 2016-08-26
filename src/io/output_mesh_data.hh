/*
 * output_mesh_data.hh
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#ifndef SRC_IO_OUTPUT_MESH_DATA_HH_
#define SRC_IO_OUTPUT_MESH_DATA_HH_

#include <ostream>
#include <string>
#include "system/asserts.hh"
#include "io/output_data_base.hh"
#include <type_traits>

/// Class representing data vector of geometry and topology information (especially for VTK).
/// Filling the vector is the users responsibility.
template <typename T>
class MeshData : public OutputDataBase {
public:
    /// Constructor. @p name is the possible name of the output vector.
    MeshData(std::string name, NumCompValueType n_elem = N_SCALAR);

    ~MeshData() override;

    /// Prints @p idx element of data vector into stream.
    void print_ascii(std::ostream& out_stream, unsigned int idx) override;

    /// Prints the whole data vector into stream.
    void print_ascii_all(std::ostream& out_stream) override;

    /**
     * Prints the whole data vector into stream in appended binary format.
     */
    void print_binary_all(ostream &out_stream) override;

    /// Prints the whole data vector into stream. UNSUPPORTED.
    void print_all_yaml(std::ostream& out_stream, unsigned int precision) override;

    /**
     * Find minimal and maximal range of stored data
     */
    void get_min_max_range(double &min, double &max) override;

    /// Access i-th element in the data vector.
    T& operator[](unsigned int i);

    /// Data vector.
    std::vector<T> data_;

};



#endif /* SRC_IO_OUTPUT_MESH_DATA_HH_ */
