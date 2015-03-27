/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
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
 * $Id: output_vtk.cc 2505 2013-09-13 14:52:27Z jiri.hnidek $
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file    output_vtk.cc
 * @brief   The functions for outputs to VTK files.
 *
 */

#include "output_vtk.hh"
#include <limits.h>
#include "mesh/mesh.h"
#include "output_data_base.hh"


using namespace Input::Type;

Record OutputVTK::input_type
    = Record("vtk", "Parameters of vtk output format.")
    // It is derived from abstract class
    .derive_from(OutputTime::input_format_type)
    .declare_key("variant", input_type_variant, Default("ascii"),
        "Variant of output stream file format.")
    // The parallel or serial variant
    .declare_key("parallel", Bool(), Default("false"),
        "Parallel or serial version of file format.")
    // Type of compression
    .declare_key("compression", input_type_compression, Default("none"),
        "Compression used in output stream file format.");


Selection OutputVTK::input_type_variant
    = Selection("VTK variant (ascii or binary)")
    .add_value(OutputVTK::VARIANT_ASCII, "ascii",
        "ASCII variant of VTK file format")
    .add_value(OutputVTK::VARIANT_BINARY, "binary",
        "Binary variant of VTK file format (not supported yet)");


Selection OutputVTK::input_type_compression
    = Selection("Type of compression of VTK file format")
    .add_value(OutputVTK::COMPRESSION_NONE, "none",
        "Data in VTK file format are not compressed")
    .add_value(OutputVTK::COMPRESSION_GZIP, "zlib",
        "Data in VTK file format are compressed using zlib (not supported yet)");


OutputVTK::OutputVTK(const Input::Record &in_rec) : OutputTime(in_rec)
{
    this->fix_main_file_extension(".pvd");

    if(this->rank == 0) {
        this->_base_file.open(this->_base_filename.c_str());
        INPUT_CHECK( this->_base_file.is_open() , "Can not open output file: %s\n", this->_base_filename.c_str() );
        xprintf(MsgLog, "Writing flow output file: %s ... \n", this->_base_filename.c_str());
    }

    this->make_subdirectory();
    this->write_head();
}



OutputVTK::OutputVTK()
{}



OutputVTK::~OutputVTK()
{
    this->write_tail();
}




int OutputVTK::write_data(void)
{
    ASSERT(_mesh != nullptr, "Null mesh.\n");

    /* It's possible now to do output to the file only in the first process */
    if(this->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
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
    ASSERT( main_file.substr( main_file.size() - 4) == ".pvd" , "none");
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

void OutputVTK::write_vtk_geometry(void)
{
    Mesh *mesh = this->_mesh;
    ofstream &file = this->_data_file;

    int tmp;

    /* Write Points begin*/
    file << "<Points>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    tmp = 0;
    /* Set floating point precision */
    file.precision(std::numeric_limits<double>::digits10);
    FOR_NODES(mesh, node) {
        node->aux = tmp;   /* store index in the auxiliary variable */

        file << scientific << node->getX() << " ";
        file << scientific << node->getY() << " ";
        file << scientific << node->getZ() << " ";

        tmp++;
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;
    /* Write Points end*/
    file << "</Points>" << endl;
}

void OutputVTK::write_vtk_topology(void)
{
    Mesh *mesh = this->_mesh;
    ofstream &file = this->_data_file;

    Node* node;
    //ElementIter ele;
    unsigned int li;
    int tmp;

    /* Write Cells begin*/
    file << "<Cells>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            file << node->aux << " ";   /* Write connectivity */
        }
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    /* Write number of nodes for each element */
    tmp = 0;
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->dim()) {
        case 1:
            tmp += VTK_LINE_SIZE;
            break;
        case 2:
            tmp += VTK_TRIANGLE_SIZE;
            break;
        case 3:
            tmp += VTK_TETRA_SIZE;
            break;
        }
        file << tmp << " ";
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    /* Write type of nodes for each element */
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->dim()) {
        case 1:
            file << (int)VTK_LINE << " ";
            break;
        case 2:
            file << (int)VTK_TRIANGLE << " ";
            break;
        case 3:
            file << (int)VTK_TETRA << " ";
            break;
        }
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write Cells end*/
    file << "</Cells>" << endl;
}

void OutputVTK::write_vtk_discont_geometry(void)
{
    Mesh *mesh = this->_mesh;
    ofstream &file = this->_data_file;

    NodeIter node;
    unsigned int li;

    /* Write Points begin*/
    file << "<Points>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Set floating point precision */
    file.precision(std::numeric_limits<double>::digits10);
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];

            file << scientific << node->getX() << " ";
            file << scientific << node->getY() << " ";
            file << scientific << node->getZ() << " ";
        }
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;
    /* Write Points end*/
    file << "</Points>" << endl;
}

void OutputVTK::write_vtk_discont_topology(void)
{
    Mesh *mesh = this->_mesh;
    ofstream &file = this->_data_file;

    //Node* node;
    //ElementIter ele;
    unsigned int li, tmp;

    /* Write Cells begin*/
    file << "<Cells>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    tmp = 0;
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            file << tmp << " ";   /* Write connectivity */
            tmp++;
        }
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    /* Write number of nodes for each element */
    tmp = 0;
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->dim()) {
        case 1:
            tmp += VTK_LINE_SIZE;
            break;
        case 2:
            tmp += VTK_TRIANGLE_SIZE;
            break;
        case 3:
            tmp += VTK_TETRA_SIZE;
            break;
        }
        file << tmp << " ";
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    /* Write type of nodes for each element */
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->dim()) {
        case 1:
            file << (int)VTK_LINE << " ";
            break;
        case 2:
            file << (int)VTK_TRIANGLE << " ";
            break;
        case 3:
            file << (int)VTK_TETRA << " ";
            break;
        }
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write Cells end*/
    file << "</Cells>" << endl;
}




void OutputVTK::write_vtk_data_ascii(OutputDataFieldVec &output_data_vec)
{
    ofstream &file = this->_data_file;

    for(OutputDataPtr data :  output_data_vec)
    {
        file 	<< "<DataArray type=\"Float64\" "
        		<< "Name=\"" << data->output_field_name <<"\" ";
        if (data->n_elem_ > 1)
        {
        	file
        		<< "NumberOfComponents=\"" << data->n_elem_ << "\" ";
        }
        file	<< "format=\"ascii\">"
                << endl;

        /* Set precision to max */
        file.precision(std::numeric_limits<double>::digits10);

        data->print_all(file);

        file << "\n</DataArray>" << endl;
    }
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
        this->write_vtk_data_ascii(output_data_vec_[NODE_DATA]);

        /* Write data in corners of elements */
        this->write_vtk_data_ascii(output_data_vec_[CORNER_DATA]);

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
    this->write_vtk_data_ascii(data_map);

    /* Write PointData end */
    file << "</CellData>" << endl;
}


void OutputVTK::write_vtk_vtu_tail(void)
{
    ofstream &file = this->_data_file;

    file << "</UnstructuredGrid>" << endl;
    file << "</VTKFile>" << endl;
}


void OutputVTK::write_vtk_vtu(void)
{
    ofstream &file = this->_data_file;
    Mesh *mesh = this->_mesh;

    /* Write header */
    this->write_vtk_vtu_head();

    /* When there is no discontinuous data, then write classical vtu */
    if ( this->output_data_vec_[CORNER_DATA].empty() )
    {
        /* Write Piece begin */
        file << "<Piece NumberOfPoints=\"" << mesh->n_nodes() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        this->write_vtk_geometry();

        /* Write VTK Topology */
        this->write_vtk_topology();

        /* Write VTK scalar and vector data on nodes to the file */
        this->write_vtk_node_data();

        /* Write VTK data on elements */
        this->write_vtk_element_data();

        /* Write Piece end */
        file << "</Piece>" << endl;

    } else {
        /* Write Piece begin */
        file << "<Piece NumberOfPoints=\"" << mesh->n_corners() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        this->write_vtk_discont_geometry();

        /* Write VTK Topology */
        this->write_vtk_discont_topology();

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





