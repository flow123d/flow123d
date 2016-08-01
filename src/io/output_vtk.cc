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
#include "output_data_base.hh"
#include "output_mesh_data.hh"
#include "output_mesh.hh"

#include <limits.h>
#include "input/factory.hh"
#include "input/accessors_forward.hh"

FLOW123D_FORCE_LINK_IN_CHILD(vtk)


using namespace Input::Type;

const Record & OutputVTK::get_input_type() {
    return Record("vtk", "Parameters of vtk output format.")
		// It is derived from abstract class
		.derive_from(OutputTime::get_input_format_type())
		.declare_key("variant", OutputVTK::get_input_type_variant(), Default("\"ascii\""),
			"Variant of output stream file format.")
		// The parallel or serial variant
		.declare_key("parallel", Bool(), Default("false"),
			"Parallel or serial version of file format.")
		.close();
}


const Selection & OutputVTK::get_input_type_variant() {
    return Selection("VTK variant (ascii or binary)")
		.add_value(OutputVTK::VARIANT_ASCII, "ascii",
			"ASCII variant of VTK file format")
		.add_value(OutputVTK::VARIANT_BINARY_UNCOMPRESSED, "binary",
			"Uncompressed binary variant of VTK file format (not supported yet)")
		.add_value(OutputVTK::VARIANT_BINARY_ZLIB, "binary_zlib",
			"Binary variant of VTK file format compressed with zlib (not supported yet)")
		.close();
}


const int OutputVTK::registrar = Input::register_class< OutputVTK >("vtk") +
		OutputVTK::get_input_type().size();



OutputVTK::OutputVTK()
{
    this->enable_refinement_ = true;
}



OutputVTK::~OutputVTK()
{
    this->write_tail();
}




int OutputVTK::write_data(void)
{
    ASSERT_PTR(output_mesh_).error();

    /* It's possible now to do output to the file only in the first process */
    if(this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    if (! this->_base_file.is_open()) {
        this->fix_main_file_extension(".pvd");

        this->_base_file.open(this->_base_filename.c_str());
        INPUT_CHECK( this->_base_file.is_open() , "Can not open output file: %s\n", this->_base_filename.c_str() );
        xprintf(MsgLog, "Writing flow output file: %s ... \n", this->_base_filename.c_str());

        this->make_subdirectory();
        this->write_head();
    }

    ostringstream ss;
    ss << main_output_basename_ << "-"
       << std::setw(6) << std::setfill('0') << this->current_step
       << ".vtu";


    std::string frame_file_name = ss.str();
    std::string frame_file_path = main_output_dir_ + "/" + main_output_basename_ + "/" + frame_file_name;

    _data_file.open(frame_file_path);
    if(_data_file.is_open() == false) {
        xprintf(Err, "Could not write output to the file: %s\n", frame_file_path.c_str());
        return 0;
    } else {
        /* Set up data file */

        xprintf(MsgLog, "%s: Writing output file %s ... ",
                __func__, this->_base_filename.c_str());


        /* Set floating point precision to max */
        this->_base_file.precision(std::numeric_limits<double>::digits10);

        /* Strip out relative path and add "base/" string */
        std::string relative_frame_file = main_output_basename_ + "/" + frame_file_name;
        this->_base_file << scientific << "<DataSet timestep=\"" << (isfinite(this->time)?this->time:0)
                << "\" group=\"\" part=\"0\" file=\"" << relative_frame_file <<"\"/>" << endl;

        xprintf(MsgLog, "O.K.\n");

        xprintf(MsgLog, "%s: Writing output (frame %d) file %s ... ", __func__,
                this->current_step, relative_frame_file.c_str());

        this->write_vtk_vtu();

        /* Close stream for file of current frame */
        _data_file.close();
        //delete data_file;
        //this->_data_file = NULL;

        xprintf(MsgLog, "O.K.\n");
    }

    return 1;
}




void OutputVTK::make_subdirectory()
{
    string main_file="./" + this->_base_filename; // guarantee that find_last_of succeeds
    ASSERT( main_file.substr( main_file.size() - 4) == ".pvd").error("Incorrect extension of VTK output file.");
    unsigned int last_sep_pos=main_file.find_last_of(DIR_DELIMITER);
    main_output_dir_=main_file.substr(2, last_sep_pos-2);
    main_output_basename_=main_file.substr(last_sep_pos+1);
    main_output_basename_=main_output_basename_.substr(0, main_output_basename_.size() - 4); // 5 = ".pvd".size() +1

    FilePath fp(main_output_dir_ + DIR_DELIMITER + main_output_basename_ + DIR_DELIMITER + "__tmp__", FilePath::output_file);
    fp.create_output_dir();
}




void OutputVTK::write_vtk_vtu_head(void)
{
    ofstream &file = this->_data_file;

    file << "<?xml version=\"1.0\"?>" << endl;
    // TODO: test endianess of platform (this would be important, when raw
    // data will be saved to the VTK file)
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    file << "<UnstructuredGrid>" << endl;
}



void OutputVTK::fill_element_types_vector(std::vector< unsigned int >& data)
{    
    auto offsets = output_mesh_->offsets_->data_;
    unsigned int n_elements = offsets.size();
    
    data.resize(n_elements);
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
}



void OutputVTK::write_vtk_data(OutputTime::OutputDataPtr output_data, VTKValueType type)
{
    ofstream &file = this->_data_file;

    file    << "<DataArray type=\"" << vtk_value_type_map(type) << "\" ";
    // possibly write name
    if( ! output_data->output_field_name.empty()) 
        file << "Name=\"" << output_data->output_field_name <<"\" ";
    // write number of components
    if (output_data->n_elem_ > 1)
    {
        file
            << "NumberOfComponents=\"" << output_data->n_elem_ << "\" ";
    }
    file    << "format=\"" << vtk_variant_map(this->variant_type_) << "\"";

    if ( this->variant_type_ == VTKVariant::VARIANT_BINARY_ZLIB ) {
    	// binary compressed output
    	file    << " offset=\"" << appended_data_.tellp() << "\"/>" << endl;
    	output_data->print_binary_all(appended_data_);
    } else {
    	file	<< ">" << endl;
    	if ( this->variant_type_ == VTKVariant::VARIANT_ASCII ) {
            // ascii output
            file.precision(std::numeric_limits<double>::digits10); // Set precision to max
            output_data->print_ascii_all(file);
    	} else {
    		// binary output
    		output_data->print_binary_all(file);
    	}
    	file << "\n</DataArray>" << endl;
    }

}


void OutputVTK::write_vtk_data(OutputDataFieldVec &output_data_vec)
{
    for(OutputDataPtr data :  output_data_vec)
        write_vtk_data(data, VTK_FLOAT64);
}




void OutputVTK::write_vtk_data_names(ofstream &file,
        OutputDataFieldVec &output_data_vec)
{
    if (output_data_vec.empty()) return;

    file << "Scalars=\"";
    for(OutputDataPtr data :  output_data_vec )
		if (data->n_elem_ == OutputDataBase::N_SCALAR) file << data->output_field_name << ",";
	file << "\" ";

    file << "Vectors=\"";
    for(OutputDataPtr data :  output_data_vec )
		if (data->n_elem_ == OutputDataBase::N_VECTOR) file << data->output_field_name << ",";
	file << "\" ";

    file << "Tensors=\"";
    for(OutputDataPtr data :  output_data_vec )
		if (data->n_elem_ == OutputDataBase::N_TENSOR) file << data->output_field_name << ",";
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
        this->write_vtk_data(output_data_vec_[NODE_DATA]);

        /* Write data in corners of elements */
        this->write_vtk_data(output_data_vec_[CORNER_DATA]);

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
    this->write_vtk_data(data_map);

    /* Write PointData end */
    file << "</CellData>" << endl;
}


void OutputVTK::write_vtk_vtu_tail(void)
{
    ofstream &file = this->_data_file;

    file << "</UnstructuredGrid>" << endl;
    if ( this->variant_type_ == VTKVariant::VARIANT_BINARY_ZLIB ) {
    	// appended data of binary compressed output
    	WarningOut() << "Zlib library is not supported yet. Appended output is not compressed." << endl;
    	file << "<UnstructuredGrid encoding=\"raw\">" << endl;
    	file << appended_data_.str() << endl;
    	file << "</UnstructuredGrid>" << endl;
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
            write_vtk_data(output_mesh_->nodes_, VTK_FLOAT64 );
        file << "</Points>" << endl;
    
        
        /* Write VTK Topology */
        file << "<Cells>" << endl;
            write_vtk_data(output_mesh_->connectivity_, VTK_INT32 );
            write_vtk_data(output_mesh_->offsets_, VTK_INT32 );
            auto types = std::make_shared<MeshData<unsigned int>>("types");
            fill_element_types_vector(types->data_);
           	write_vtk_data( types, VTK_UINT8 );
        file << "</Cells>" << endl;

        /* Write VTK scalar and vector data on nodes to the file */
        this->write_vtk_node_data();

        /* Write VTK data on elements */
        this->write_vtk_element_data();

        /* Write Piece end */
        file << "</Piece>" << endl;

    } else {
        /* Write Piece begin */
        file << "<Piece NumberOfPoints=\"" << output_mesh_discont_->n_nodes()
                  << "\" NumberOfCells=\"" << output_mesh_->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        file << "<Points>" << endl;
            write_vtk_data(output_mesh_discont_->nodes_, VTK_FLOAT64 );
        file << "</Points>" << endl;

        /* Write VTK Topology */
        file << "<Cells>" << endl;
            write_vtk_data(output_mesh_discont_->connectivity_, VTK_INT32 );
            write_vtk_data(output_mesh_discont_->offsets_, VTK_INT32 );
            auto types = std::make_shared<MeshData<unsigned int>>("types");
            fill_element_types_vector(types->data_);
           	write_vtk_data( types, VTK_UINT8 );
        file << "</Cells>" << endl;

        /* Write VTK scalar and vector data on nodes to the file */
        this->write_vtk_node_data();

        /* Write VTK data on elements */
        this->write_vtk_element_data();

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

    xprintf(MsgLog, "%s: Writing output file (head) %s ... ", __func__,
            this->_base_filename.c_str() );

    this->_base_file << "<?xml version=\"1.0\"?>" << endl;
    this->_base_file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    this->_base_file << "<Collection>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}


int OutputVTK::write_tail(void)
{
    /* It's possible now to do output to the file only in the first process */
    if(this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    xprintf(MsgLog, "%s: Writing output file (tail) %s ... ", __func__,
            this->_base_filename.c_str() );

    this->_base_file << "</Collection>" << endl;
    this->_base_file << "</VTKFile>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}






