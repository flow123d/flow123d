/*!
 *
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

#include <limits.h>
#include <mpi.h>
#include <boost/any.hpp>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <assert.h>

#include "system/xio.h"
#include "io/output.h"
#include "io/output_vtk.h"
#include "mesh/mesh.h"


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



void OutputVTK::write_vtk_vtu_head(void)
{
    ofstream &file = this->get_data_file();

    file << "<?xml version=\"1.0\"?>" << endl;
    // TODO: test endianess of platform (this would be important, when raw
    // data will be saved to the VTK file)
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    file << "<UnstructuredGrid>" << endl;
}

void OutputVTK::write_vtk_geometry(void)
{
    Mesh *mesh = this->get_mesh();
    ofstream &file = this->get_data_file();

    //NodeIter node;
    int tmp;

    /* Write Points begin*/
    file << "<Points>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    tmp = 0;
    /* Set floating point precision */
    this->get_data_file().precision(std::numeric_limits<double>::digits10);
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
    Mesh *mesh = this->get_mesh();
    ofstream &file = this->get_data_file();

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
    Mesh *mesh = this->get_mesh();
    ofstream &file = this->get_data_file();

    NodeIter node;
    unsigned int li;

    /* Write Points begin*/
    file << "<Points>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Set floating point precision */
    this->get_data_file().precision(std::numeric_limits<double>::digits10);
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
    Mesh *mesh = this->get_mesh();
    ofstream &file = this->get_data_file();

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




void OutputVTK::write_vtk_data_ascii(vector<OutputDataBase*> &vec_output_data)
{
    ofstream &file = this->get_data_file();

    for( OutputDataBase* data : vec_output_data)
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
        //file.precision(std::numeric_limits<float>::digits10);
        file.precision(std::numeric_limits<double>::digits10);

        for(unsigned int i=0; i < data->n_values; i ++) {
            data->print(file, i);
        }

        file << "\n</DataArray>" << endl;
    }
}




void OutputVTK::write_vtk_data_names(ofstream &file, vector<OutputDataBase*> &vec_output_data)
{
    file << "Scalars=\"";
	for( auto &data : vec_output_data)
		if (data->n_elem_ == OutputDataBase::scalar) file << data->output_field_name << ",";
	file << "\" ";

    file << "Vectors=\"";
	for( auto &data : vec_output_data)
		if (data->n_elem_ == OutputDataBase::vector) file << data->output_field_name << ",";
	file << "\" ";

    file << "Tensors=\"";
	for( auto &data : vec_output_data)
		if (data->n_elem_ == OutputDataBase::tensor) file << data->output_field_name << ",";
	file << "\"";
}


void OutputVTK::write_vtk_node_data(void)
{
    ofstream &file = this->get_data_file();

    // merge node and corner data
    vector<OutputDataBase*> node_corner_data(this->node_data);
    node_corner_data.insert(node_corner_data.end(), this->corner_data.begin(), this->corner_data.end());

    if( ! node_corner_data.empty() ) {
        /* Write <PointData begin */
        file << "<PointData ";
        write_vtk_data_names(file, node_corner_data);
        file << ">" << endl;

        /* Write data on nodes */
        if( ! this->node_data.empty() ) {
            this->write_vtk_data_ascii(this->node_data);
        }

        /* Write data in corners of elements */
        if( ! this->corner_data.empty() ) {
            this->write_vtk_data_ascii(this->corner_data);
        }

        /* Write PointData end */
        file << "</PointData>" << endl;
    }
}


void OutputVTK::write_vtk_element_data(void)
{
    ofstream &file = this->get_data_file();

    if(this->elem_data.empty() != true) {
        /* Write PointData begin */
        file << "<CellData ";
        write_vtk_data_names(file, this->elem_data);
        file << ">" << endl;

        /* Write own data */
        this->write_vtk_data_ascii(this->elem_data);

        /* Write PointData end */
        file << "</CellData>" << endl;
    }
}


void OutputVTK::write_vtk_vtu_tail(void)
{
    ofstream &file = this->get_data_file();

    file << "</UnstructuredGrid>" << endl;
    file << "</VTKFile>" << endl;
}


void OutputVTK::write_vtk_vtu(void)
{
    ofstream &file = this->get_data_file();
    Mesh *mesh = this->get_mesh();

    /* Write header */
    this->write_vtk_vtu_head();

    /* When there is no discontinuous data, then write classical vtu */
    if(this->corner_data.empty() == true)
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
        file << "<Piece NumberOfPoints=\"" << this->get_corner_count() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

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


int OutputVTK::write_data(void)
{
    //Mesh *mesh = this->output_time->get_mesh();
    char base_dir_name[PATH_MAX];
    char new_dir_name[PATH_MAX];
    char base_file_name[PATH_MAX];
    char base[PATH_MAX];
    char frame_file_name[PATH_MAX];
    ofstream *data_file = new ofstream;
    DIR *dir;
    int i, j, ret;
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    // Write header with mesh, when it hasn't been written to output file yet
    if(this->header_written == false) {
        this->write_head();
        this->header_written = true;
    }

    strncpy(base_dir_name, this->base_filename()->c_str(), PATH_MAX);

    /* Remove last file name from base_name and find position of last directory
     * delimiter: '/' */
    for(j=strlen(base_dir_name)-1; j>0; j--) {
        if(base_dir_name[j]=='/') {
            base_dir_name[j]='\0';
            break;
        }
    }

    strncpy(base_file_name, this->base_filename()->c_str(), PATH_MAX);

    /* Find, where is the '.' character of .pvd suffix of base_name */
    for(i=strlen(base_file_name)-1; i>=0; i--) {
        if(base_file_name[i]=='.') {
            break;
        }
    }

    /* Create base of pvd file. Example ./output/transport.pvd -> transport */
    strncpy(base, &base_file_name[j+1], i-j-1);
    base[i-j-1]='\0';

    /* New folder for output */
    sprintf(new_dir_name, "%s/%s", base_dir_name, base);

    /* Try to open directory */
    dir = opendir(new_dir_name);
    if(dir == NULL) {
        /* Directory doesn't exist. Create new one. */
        ret = mkdir(new_dir_name, 0777);

        if(ret != 0) {
            xprintf(Err, "Couldn't create directory: %s, error: %s\n", new_dir_name, strerror(errno));
        }
    } else {
        closedir(dir);
    }

    sprintf(frame_file_name, "%s/%s-%06d.vtu", new_dir_name, base, this->current_step);

    data_file->open(frame_file_name);
    if(data_file->is_open() == false) {
        xprintf(Err, "Could not write output to the file: %s\n", frame_file_name);
        return 0;
    } else {
        /* Set up data file */
        this->set_data_file(data_file);

        xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->base_filename()->c_str());

        /* Find first directory delimiter */
        for(i=strlen(frame_file_name); i>=0; i--) {
            if(frame_file_name[i]==DIR_DELIMITER) {
                break;
            }
        }

        /* Set floating point precision to max */
        this->get_base_file().precision(std::numeric_limits<double>::digits10);

        /* Strip out relative path and add "base/" string */
        this->get_base_file() << scientific << "<DataSet timestep=\"" << (isfinite(this->time)?this->time:0) << "\" group=\"\" part=\"0\" file=\"" << base << "/" << &frame_file_name[i+1] <<"\"/>" << endl;

        xprintf(MsgLog, "O.K.\n");

        xprintf(MsgLog, "%s: Writing output (frame %d) file %s ... ", __func__,
                this->current_step, frame_file_name);

        this->write_vtk_vtu();

        /* Close stream for file of current frame */
        data_file->close();
        delete data_file;
        this->set_data_file(NULL);

        xprintf(MsgLog, "O.K.\n");
    }

    return 1;
}


int OutputVTK::write_head(void)
{
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    xprintf(MsgLog, "%s: Writing output file (head) %s ... ", __func__,
            this->base_filename()->c_str() );

    this->get_base_file() << "<?xml version=\"1.0\"?>" << endl;
    this->get_base_file() << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    this->get_base_file() << "<Collection>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}


int OutputVTK::write_tail(void)
{
    int rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    xprintf(MsgLog, "%s: Writing output file (tail) %s ... ", __func__,
            this->base_filename()->c_str() );

    this->get_base_file() << "</Collection>" << endl;
    this->get_base_file() << "</VTKFile>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}


OutputVTK::OutputVTK(const Input::Record &in_rec) : OutputTime(in_rec)
{
	this->file_format = OutputTime::VTK;
    this->header_written = false;
}


OutputVTK::~OutputVTK()
{
    this->write_tail();
}

