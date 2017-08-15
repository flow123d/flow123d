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
 * @file    output_vtk.cc
 * @brief   The functions for outputs to VTK files.
 */

#include "output_vtk.hh"
#include "element_data_cache_base.hh"
#include "output_mesh.hh"

#include <limits.h>
#include "input/factory.hh"
#include "input/accessors_forward.hh"
#include "system/file_path.hh"

#include "config.h"
#include <zlib.h>

FLOW123D_FORCE_LINK_IN_CHILD(vtk)


using namespace Input::Type;

const Record & OutputVTK::get_input_type() {
    return Record("vtk", "Parameters of vtk output format.")
		// It is derived from abstract class
		.derive_from(OutputTime::get_input_format_type())
		.declare_key("variant", OutputVTK::get_input_type_variant(), Default("\"ascii\""),
			"Variant of output stream file format.")
		// The parallel or serial variant
		//.declare_key("parallel", Bool(), Default("false"),
		//	"Parallel or serial version of file format.")
		.close();
}


const Selection & OutputVTK::get_input_type_variant() {
    return Selection("VTK variant (ascii or binary)")
		.add_value(OutputVTK::VARIANT_ASCII, "ascii",
			"ASCII variant of VTK file format")
		.add_value(OutputVTK::VARIANT_BINARY_UNCOMPRESSED, "binary",
			"Uncompressed appended binary XML VTK format without usage of base64 encoding of appended data.")
#ifdef FLOW123D_HAVE_ZLIB
		.add_value(OutputVTK::VARIANT_BINARY_ZLIB, "binary_zlib",
			"Appended binary XML VTK format without usage of base64 encoding of appended data. Compressed with ZLib.")
#endif // FLOW123D_HAVE_ZLIB
		.close();
}


const int OutputVTK::registrar = Input::register_class< OutputVTK >("vtk") +
		OutputVTK::get_input_type().size();


const std::vector<std::string> OutputVTK::formats = { "ascii", "appended", "appended" };



OutputVTK::OutputVTK()
{
    this->enable_refinement_ = true;
}



OutputVTK::~OutputVTK()
{
    this->write_tail();
}



void OutputVTK::init_from_input(const std::string &equation_name, Mesh &mesh, const Input::Record &in_rec)
{
	OutputTime::init_from_input(equation_name, mesh, in_rec);

    if(this->rank == 0) {
        auto format_rec = (Input::Record)(input_record_.val<Input::AbstractRecord>("format"));
        variant_type_ = format_rec.val<VTKVariant>("variant");

        this->fix_main_file_extension(".pvd");
        try {
            this->_base_filename.open_stream( this->_base_file );
        } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, input_record_)

        LogOut() << "Writing flow output file: " << this->_base_filename << " ... ";

        this->make_subdirectory();
        this->write_head();
    }

}



int OutputVTK::write_data(void)
{
    ASSERT_PTR(output_mesh_).error();

    /* It's possible now to do output to the file only in the first process */
    if(this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    ASSERT(this->_base_file.is_open())(this->_base_filename).error();

    ostringstream ss;
    ss << main_output_basename_ << "-"
       << std::setw(6) << std::setfill('0') << this->current_step
       << ".vtu";


    std::string frame_file_name = ss.str();
    FilePath frame_file_path({main_output_dir_, main_output_basename_, frame_file_name}, FilePath::output_file);

    /* Set up data file */
    try {
        frame_file_path.open_stream(_data_file);
    } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, input_record_)


    LogOut() << __func__ << ": Writing output file " << this->_base_filename << " ... ";

    /* Set floating point precision to max */
    this->_base_file.precision(std::numeric_limits<double>::digits10);

    /* Strip out relative path and add "base/" string */
    std::string relative_frame_file = main_output_basename_ + "/" + frame_file_name;
    this->_base_file << scientific << "<DataSet timestep=\"" << (isfinite(this->time)?this->time:0)
            << "\" group=\"\" part=\"0\" file=\"" << relative_frame_file <<"\"/>" << endl;

    LogOut() << "O.K.";

    LogOut() << __func__ << ": Writing output (frame " << this->current_step << ") file " << relative_frame_file << " ... ";

    this->write_vtk_vtu();

    /* Close stream for file of current frame */
    _data_file.close();
    //delete data_file;
    //this->_data_file = NULL;

    LogOut() << "O.K.";

    return 1;
}




void OutputVTK::make_subdirectory()
{
	ASSERT_EQ(this->_base_filename.extension(), ".pvd").error();
	main_output_dir_ = this->_base_filename.parent_path();
	main_output_basename_ = this->_base_filename.stem();

    vector<string> sub_path = { main_output_dir_, main_output_basename_, "__tmp__" };
    FilePath fp(sub_path, FilePath::output_file);
    fp.create_output_dir();
}




void OutputVTK::write_vtk_vtu_head(void)
{
    ofstream &file = this->_data_file;

    file << "<?xml version=\"1.0\"?>" << endl;
    // TODO: test endianess of platform (this would be important, when raw
    // data will be saved to the VTK file)
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"";
    if ( this->variant_type_ != VTKVariant::VARIANT_ASCII ) {
    	file << " header_type=\"UInt64\"";
    }
    if ( this->variant_type_ == VTKVariant::VARIANT_BINARY_ZLIB ) {
    	file << " compressor=\"vtkZLibDataCompressor\"";
    }
    file << ">" << endl;
    file << "<UnstructuredGrid>" << endl;
}



std::shared_ptr<ElementDataCache<unsigned int>> OutputVTK::fill_element_types_data()
{    
    auto &offsets = *( output_mesh_->offsets_->get_component_data(0).get() );
    unsigned int n_elements = offsets.size();
    
    auto types = std::make_shared<ElementDataCache<unsigned int>>("types", (unsigned int)ElementDataCacheBase::N_SCALAR, 1, n_elements);
    std::vector< unsigned int >& data = *( types->get_component_data(0).get() );
    int n_nodes;
    
    n_nodes = offsets[0];
    switch(n_nodes) {
        case 2:
            data[0] = (unsigned int)VTK_LINE;
            break;
        case 3:
            data[0] = (unsigned int)VTK_TRIANGLE;
            break;
        case 4:
            data[0] = (unsigned int)VTK_TETRA;
            break;
        }
    
    for(unsigned int i=1; i < n_elements; i++)
    {
        n_nodes = offsets[i]-offsets[i-1];
        switch(n_nodes) {
        case 2:
            data[i] = (unsigned int)VTK_LINE;
            break;
        case 3:
            data[i] = (unsigned int)VTK_TRIANGLE;
            break;
        case 4:
            data[i] = (unsigned int)VTK_TETRA;
            break;
        }
    }

    return types;
}



void OutputVTK::write_vtk_data(OutputTime::OutputDataPtr output_data)
{
    // names of types in DataArray section
	static const std::vector<std::string> types = {
        "Int8", "UInt8", "Int16", "UInt16", "Int32", "UInt32", "Float32", "Float64" };

    ofstream &file = this->_data_file;

    file    << "<DataArray type=\"" << types[output_data->vtk_type()] << "\" ";
    // possibly write name
    if( ! output_data->field_input_name().empty())
        file << "Name=\"" << output_data->field_input_name() <<"\" ";
    // write number of components
    if (output_data->n_elem() > 1)
    {
        file
            << "NumberOfComponents=\"" << output_data->n_elem() << "\" ";
    }
    file    << "format=\"" << formats[this->variant_type_] << "\"";

    if ( this->variant_type_ == VTKVariant::VARIANT_ASCII ) {
    	// ascii output
    	file << ">" << endl;
    	file << std::fixed << std::setprecision(10); // Set precision to max
    	output_data->print_ascii_all(file);
    	file << "\n</DataArray>" << endl;
    } else {
    	// binary output is stored to appended_data_ stream
    	double range_min, range_max;
    	output_data->get_min_max_range(range_min, range_max);
    	file    << " offset=\"" << appended_data_.tellp() << "\" ";
    	file    << "RangeMin=\"" << range_min << "\" RangeMax=\"" << range_max << "\"/>" << endl;
    	if ( this->variant_type_ == VTKVariant::VARIANT_BINARY_UNCOMPRESSED ) {
    		output_data->print_binary_all( appended_data_ );
    	} else { // ZLib compression
    		stringstream uncompressed_data, compressed_data;
    		output_data->print_binary_all( uncompressed_data, false );
    		this->compress_data(uncompressed_data, compressed_data);
    		appended_data_ << compressed_data.str();
    	}
    }

}


void OutputVTK::compress_data(stringstream &uncompressed_stream, stringstream &compressed_stream) {
    // size of block of compressed data.
	static const size_t BUF_SIZE = 32 * 1024;

	string uncompressed_string = uncompressed_stream.str();  // full uncompressed string
	uLong uncompressed_size = uncompressed_string.size();    // size of uncompressed string
	stringstream compressed_data;                            // helper stream stores blocks of compress data

	uLong count_of_blocks = (uncompressed_size + BUF_SIZE - 1) / BUF_SIZE;
	uLong last_block_size = (uncompressed_size % BUF_SIZE);
	compressed_stream.write(reinterpret_cast<const char*>(&count_of_blocks), sizeof(unsigned long long int));
	compressed_stream.write(reinterpret_cast<const char*>(&BUF_SIZE), sizeof(unsigned long long int));
	compressed_stream.write(reinterpret_cast<const char*>(&last_block_size), sizeof(unsigned long long int));

	for (uLong i=0; i<count_of_blocks; ++i) {
		// get block of data for compression
		std::string data_block = uncompressed_string.substr(i*BUF_SIZE, BUF_SIZE);
		uLong data_block_size = data_block.size();

		std::vector<uint8_t> buffer;
		uint8_t temp_buffer[BUF_SIZE];

		// set zlib object
		z_stream strm;
		strm.zalloc = 0;
		strm.zfree = 0;
		strm.next_in = reinterpret_cast<uint8_t *>(&data_block[0]);
		strm.avail_in = data_block_size;
		strm.next_out = temp_buffer;
		strm.avail_out = BUF_SIZE;

		// compression of data
		deflateInit(&strm, Z_BEST_COMPRESSION);
		while (strm.avail_in != 0) {
			int res = deflate(&strm, Z_NO_FLUSH);
			ASSERT_EQ(res, Z_OK).error();
			if (strm.avail_out == 0) {
				buffer.insert(buffer.end(), temp_buffer, temp_buffer + BUF_SIZE);
				strm.next_out = temp_buffer;
				strm.avail_out = BUF_SIZE;
			}
		}
		int deflate_res = Z_OK;
		while (deflate_res == Z_OK) {
			if (strm.avail_out == 0) {
				buffer.insert(buffer.end(), temp_buffer, temp_buffer + BUF_SIZE);
				strm.next_out = temp_buffer;
				strm.avail_out = BUF_SIZE;
			}
			deflate_res = deflate(&strm, Z_FINISH);
		}
		ASSERT_EQ(deflate_res, Z_STREAM_END).error();
		buffer.insert(buffer.end(), temp_buffer, temp_buffer + BUF_SIZE - strm.avail_out);
		deflateEnd(&strm);

		// store compress data and its size to streams
		std::string str(buffer.begin(), buffer.end());
		uLong compressed_data_size = str.size();
		compressed_stream.write(reinterpret_cast<const char*>(&compressed_data_size), sizeof(unsigned long long int));
		compressed_data << str;
	}
	// push compress data to returned stream
	compressed_stream << compressed_data.str();
}


void OutputVTK::write_vtk_field_data(OutputDataFieldVec &output_data_vec)
{
    for(OutputDataPtr data :  output_data_vec)
        write_vtk_data(data);
}




void OutputVTK::write_vtk_data_names(ofstream &file,
        OutputDataFieldVec &output_data_vec)
{
    if (output_data_vec.empty()) return;

    file << "Scalars=\"";
    for(OutputDataPtr data :  output_data_vec )
		if (data->n_elem() == ElementDataCacheBase::N_SCALAR) file << data->field_input_name() << ",";
	file << "\" ";

    file << "Vectors=\"";
    for(OutputDataPtr data :  output_data_vec )
		if (data->n_elem() == ElementDataCacheBase::N_VECTOR) file << data->field_input_name() << ",";
	file << "\" ";

    file << "Tensors=\"";
    for(OutputDataPtr data :  output_data_vec )
		if (data->n_elem() == ElementDataCacheBase::N_TENSOR) file << data->field_input_name() << ",";
	file << "\"";
}


void OutputVTK::write_vtk_node_data(void)
{
    ofstream &file = this->_data_file;

    // merge node and corner data
    OutputDataFieldVec node_corner_data(output_data_vec_[NODE_DATA]);
    node_corner_data.insert(node_corner_data.end(),
            output_data_vec_[CORNER_DATA].begin(), output_data_vec_[CORNER_DATA].end());

    if( ! node_corner_data.empty() ) {
        /* Write <PointData begin */
        file << "<PointData ";
        write_vtk_data_names(file, node_corner_data);
        file << ">" << endl;

        /* Write data on nodes */
        this->write_vtk_field_data(output_data_vec_[NODE_DATA]);

        /* Write data in corners of elements */
        this->write_vtk_field_data(output_data_vec_[CORNER_DATA]);

        /* Write PointData end */
        file << "</PointData>" << endl;
    }
}


void OutputVTK::write_vtk_element_data(void)
{
    ofstream &file = this->_data_file;

    auto &data_map = this->output_data_vec_[ELEM_DATA];
    if (data_map.empty()) return;

    /* Write CellData begin */
    file << "<CellData ";
    write_vtk_data_names(file, data_map);
    file << ">" << endl;

    /* Write own data */
    this->write_vtk_field_data(data_map);

    /* Write PointData end */
    file << "</CellData>" << endl;
}


void OutputVTK::write_vtk_native_data(void)
{
    ofstream &file = this->_data_file;

    auto &data_map = this->output_data_vec_[NATIVE_DATA];
    if (data_map.empty()) return;

    /* Write Flow123dData begin */
    file << "<Flow123dData ";
    write_vtk_data_names(file, data_map);
    file << ">" << endl;

    /* Write own data */
    for(OutputDataPtr output_data : data_map) {
        file  << "<DataArray type=\"Float64\" ";
        file  << "Name=\"" << output_data->field_input_name() <<"\" ";
        file  << "format=\"" << formats[this->variant_type_] << "\" ";
        file  << "dof_handler_hash=\"" << output_data->dof_handler_hash() << "\" ";
        file  << "n_dofs_per_element=\"" << output_data->n_elem() << "\"";

        if ( this->variant_type_ == VTKVariant::VARIANT_ASCII ) {
        	// ascii output
        	file << ">" << endl;
        	file << std::fixed << std::setprecision(10); // Set precision to max
        	output_data->print_ascii_all(file);
        	file << "\n</DataArray>" << endl;
        } else {
        	// binary output is stored to appended_data_ stream
        	double range_min, range_max;
        	output_data->get_min_max_range(range_min, range_max);
        	file    << " offset=\"" << appended_data_.tellp() << "\" ";
        	file    << "RangeMin=\"" << range_min << "\" RangeMax=\"" << range_max << "\"/>" << endl;
        	if ( this->variant_type_ == VTKVariant::VARIANT_BINARY_UNCOMPRESSED ) {
        		output_data->print_binary_all( appended_data_ );
        	} else { // ZLib compression
        		stringstream uncompressed_data, compressed_data;
        		output_data->print_binary_all( uncompressed_data, false );
        		this->compress_data(uncompressed_data, compressed_data);
        		appended_data_ << compressed_data.str();
        	}
        }
    }

    /* Write Flow123dData end */
    file << "</Flow123dData>" << endl;
}


void OutputVTK::write_vtk_vtu_tail(void)
{
    ofstream &file = this->_data_file;

    file << "</UnstructuredGrid>" << endl;
    if ( this->variant_type_ != VTKVariant::VARIANT_ASCII ) {
    	// appended data of binary compressed output
    	file << "<AppendedData encoding=\"raw\">" << endl;
    	// appended data starts with '_' character
    	file << "_" << appended_data_.str() << endl;
    	file << "</AppendedData>" << endl;
    }
    file << "</VTKFile>" << endl;
}


void OutputVTK::write_vtk_vtu(void)
{
    ofstream &file = this->_data_file;

    /* Write header */
    this->write_vtk_vtu_head();

    /* When there is no discontinuous data, then write classical vtu */
    if ( this->output_data_vec_[CORNER_DATA].empty() )
    {
        /* Write Piece begin */
        file << "<Piece NumberOfPoints=\"" << output_mesh_->n_nodes()
                  << "\" NumberOfCells=\"" << output_mesh_->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        file << "<Points>" << endl;
            write_vtk_data(output_mesh_->nodes_);
        file << "</Points>" << endl;
    
        
        /* Write VTK Topology */
        file << "<Cells>" << endl;
            write_vtk_data(output_mesh_->connectivity_);
            write_vtk_data(output_mesh_->offsets_);
            auto types = fill_element_types_data();
           	write_vtk_data( types );
        file << "</Cells>" << endl;

        /* Write VTK scalar and vector data on nodes to the file */
        this->write_vtk_node_data();

        /* Write VTK data on elements */
        this->write_vtk_element_data();

        /* Write own VTK native data (skipped by Paraview) */
        this->write_vtk_native_data();

        /* Write Piece end */
        file << "</Piece>" << endl;

    } else {
        /* Write Piece begin */
        file << "<Piece NumberOfPoints=\"" << output_mesh_discont_->n_nodes()
                  << "\" NumberOfCells=\"" << output_mesh_->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        file << "<Points>" << endl;
            write_vtk_data(output_mesh_discont_->nodes_);
        file << "</Points>" << endl;

        /* Write VTK Topology */
        file << "<Cells>" << endl;
            write_vtk_data(output_mesh_discont_->connectivity_);
            write_vtk_data(output_mesh_discont_->offsets_);
            auto types = fill_element_types_data();
           	write_vtk_data( types );
        file << "</Cells>" << endl;

        /* Write VTK scalar and vector data on nodes to the file */
        this->write_vtk_node_data();

        /* Write VTK data on elements */
        this->write_vtk_element_data();

        /* Write own VTK native data (skipped by Paraview) */
        this->write_vtk_native_data();

        /* Write Piece end */
        file << "</Piece>" << endl;
    }

    /* Write tail */
    this->write_vtk_vtu_tail();
}



int OutputVTK::write_head(void)
{
    /* It's possible now to do output to the file only in the first process */
    if(this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    LogOut() << __func__ << ": Writing output file (head) " << this->_base_filename << " ... ";

    this->_base_file << "<?xml version=\"1.0\"?>" << endl;
    this->_base_file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    this->_base_file << "<Collection>" << endl;

    LogOut() << "O.K.";

    return 1;
}


int OutputVTK::write_tail(void)
{
    /* It's possible now to do output to the file only in the first process */
    if(this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    LogOut() << __func__ << ": Writing output file (tail) " << this->_base_filename << " ... ";

    this->_base_file << "</Collection>" << endl;
    this->_base_file << "</VTKFile>" << endl;

    LogOut() << "O.K.";

    return 1;
}






