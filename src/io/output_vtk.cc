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
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file    output_vtk.cc
 * @brief   The functions for outputs to VTK files.
 *
 */

#include <limits.h>

#include "system/xio.h"
#include "io/output.h"
#include "io/output_vtk.h"
#include "mesh/mesh.h"

// TODO: v tomto souboru se poflakuji vselijake funkce mimo tridy, to v objektovem navrhu nema co delat
// maji to byt privatni metody nejake tridy (OutputVTK)

/**
 * \brief Write header of VTK file (.vtu)
 * \param[in]	out	The output file
 */
static void write_vtk_vtu_head(Output *output)
{
    output->get_data_file() << "<?xml version=\"1.0\"?>" << endl;
    // TODO: test endianess of platform (this would be important, when raw
    // data will be saved to the VTK file)
    output->get_data_file() << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    output->get_data_file() << "<UnstructuredGrid>" << endl;
}

/**
 * \brief Write geometry (position of nodes) to the VTK file (.vtu)
 * \param[in] *out The output file
 */
static void write_vtk_geometry(Output *output)
{
    Mesh *mesh = output->get_mesh();

    NodeIter node;
    int tmp;

    /* Write Points begin*/
    output->get_data_file() << "<Points>" << endl;
    /* Write DataArray begin */
    output->get_data_file() << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    tmp = 0;
    /* Set floating point precision */
    output->get_data_file().precision(std::numeric_limits<double>::digits10);
    FOR_NODES(mesh, node ) {
        node->aux = tmp;   /* store index in the auxiliary variable */

        output->get_data_file() << scientific << node->getX() << " ";
        output->get_data_file() << scientific << node->getY() << " ";
        output->get_data_file() << scientific << node->getZ() << " ";

        tmp++;
    }
    /* Write DataArray end */
    output->get_data_file() << endl << "</DataArray>" << endl;
    /* Write Points end*/
    output->get_data_file() << "</Points>" << endl;
}

/**
 * \brief Write topology (connection of nodes) to the VTK file (.vtu)
 * \param[in] *out The output file
 */
static void write_vtk_topology(Output *output)
{
    //Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    Mesh *mesh = output->get_mesh();

    Node* node;
    ElementIter ele;
    int li, tmp;

    /* Write Cells begin*/
    output->get_data_file() << "<Cells>" << endl;
    /* Write DataArray begin */
    output->get_data_file() << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    FOR_ELEMENTS(mesh, ele) {
        FOR_ELEMENT_NODES(ele, li) {
            node = ele->node[li];
            output->get_data_file() << node->aux << " ";   /* Write connectivity */
        }
    }
    /* Write DataArray end */
    output->get_data_file() << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    output->get_data_file() << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    /* Write number of nodes for each element */
    tmp = 0;
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->type) {
        case LINE:
            tmp += VTK_LINE_SIZE;
            break;
        case TRIANGLE:
            tmp += VTK_TRIANGLE_SIZE;
            break;
        case TETRAHEDRON:
            tmp += VTK_TETRA_SIZE;
            break;
        }
        output->get_data_file() << tmp << " ";
    }
    /* Write DataArray end */
    output->get_data_file() << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    output->get_data_file() << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    /* Write type of nodes for each element */
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->type) {
        case LINE:
            output->get_data_file() << (int)VTK_LINE << " ";
            break;
        case TRIANGLE:
            output->get_data_file() << (int)VTK_TRIANGLE << " ";
            break;
        case TETRAHEDRON:
            output->get_data_file() << (int)VTK_TETRA << " ";
            break;
        }
    }
    /* Write DataArray end */
    output->get_data_file() << endl << "</DataArray>" << endl;

    /* Write Cells end*/
    output->get_data_file() << "</Cells>" << endl;
}

/**
 * \brief This function writes ascii data to VTK (.vtu) output file.
 * \param[in]   *output     The pointer at Output object
 * \param[in]   *out_data   The pointer at structure storing pointer at own data.
 */
static void write_vtk_ascii_data(Output *output, OutputData *out_data)
{
    ofstream &file = output->get_data_file();

    switch(out_data->type) {
    case OUT_VECTOR_INT_SCA:
        for( std::vector<int>::iterator item = ((std::vector<int>*)out_data->data)->begin();
                item != ((std::vector<int>*)out_data->data)->end();
                ++item) {
            file << *item << " ";
        }
        break;
    case OUT_VECTOR_INT_VEC:
        for( std::vector< vector<int> >::iterator vec = ((std::vector< vector<int> >*)out_data->data)->begin();
                vec != ((std::vector< vector<int> >*)out_data->data)->end();
                ++vec) {
            for (std::vector<int>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << "  ";
        }
        break;
    case OUT_VECTOR_FLOAT_SCA:
        file.precision(std::numeric_limits<float>::digits10);
        for( std::vector<float>::iterator item = ((std::vector<float>*)out_data->data)->begin();
                item != ((std::vector<float>*)out_data->data)->end();
                ++item) {
            file << scientific << *item << " ";
        }
        break;
    case OUT_VECTOR_FLOAT_VEC:
        file.precision(std::numeric_limits<float>::digits10);
        for( std::vector< vector<float> >::iterator vec = ((std::vector< vector<float> >*)out_data->data)->begin();
                vec != ((std::vector< vector<float> >*)out_data->data)->end();
                ++vec) {
            for (std::vector<float>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << scientific << *item << " ";
            }
            file << "  ";
        }
        break;
    case OUT_VECTOR_DOUBLE_SCA:
        file.precision(std::numeric_limits<double>::digits10);
        for( std::vector<double>::iterator item = ((std::vector<double>*)out_data->data)->begin();
                item != ((std::vector<double>*)out_data->data)->end();
                ++item) {
            file << scientific << *item << " ";
        }
        break;
    case OUT_VECTOR_DOUBLE_VEC:
        file.precision(std::numeric_limits<double>::digits10);
        for( std::vector< vector<double> >::iterator vec = ((std::vector< vector<double> >*)out_data->data)->begin();
                vec != ((std::vector< vector<double> >*)out_data->data)->end();
                ++vec) {
            for (std::vector<double>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << scientific << *item << " ";
            }
            file << "  ";
        }
        break;
    case OUT_ARRAY_INT_SCA:
        for(int i=0; i<out_data->num; i++) {
            file << ((int*)out_data->data)[i] << " ";
        }
        break;
    case OUT_ARRAY_FLOAT_SCA:
        file.precision(std::numeric_limits<float>::digits10);
        for(int i=0; i<out_data->num; i++) {
            file << scientific << ((float*)out_data->data)[i] << " ";
        }
        break;
    case OUT_ARRAY_DOUBLE_SCA:
        file.precision(std::numeric_limits<double>::digits10);
        for(int i=0; i<out_data->num; i++) {
            file << scientific << ((double*)out_data->data)[i] << " ";
        }
        break;
    default:
        xprintf(Err, "This type of data: %d is not supported by VTK file format\n", out_data->type);
        break;
    }
}

/**
 * \brief Write scalar data to the VTK file (.vtu)
 * \param[in]	out		The output file
 */
static void write_vtk_scalar_ascii(Output *output, OutputData *data)
{
    /* Write DataArray begin */
    output->get_data_file() << "<DataArray type=\"Float64\" Name=\"" << *data->getName() << "_[" << *data->getUnits() <<"]\" format=\"ascii\">" << endl;//, name);
    /* Write own data */
    write_vtk_ascii_data(output, data);

    /* Write DataArray end */
    output->get_data_file() << endl << "</DataArray>" << endl;
}

/**
 * \brief Write vector data to VTK file (.vtu)
 * \param[in]	out		The output file
 */
static void write_vtk_vector_ascii(Output *output, OutputData *data)
{
    /* Write DataArray begin */
    output->get_data_file() << "<DataArray type=\"Float64\" Name=\"" << *data->getName() << "_[" << *data->getUnits() << "]\" NumberOfComponents=\"" << data->getCompNum() << "\" format=\"ascii\">" << endl;

    /* Write own data */
    write_vtk_ascii_data(output, data);

    /* Write DataArray end */
    output->get_data_file() << endl << "</DataArray>" << endl;
}

/**
 * \brief Go through all vectors of scalars and vectors and call functions that
 * write these data to VTK file (.vtu)
 * \param[in]	*out		The output file
 */
static void write_vtk_data_ascii(Output *output, std::vector<OutputData> *data)
{
    /* Write data on nodes or elements */
    if(data != NULL) {
        for(OutputDataVec::iterator dta = data->begin();
                dta != data->end(); ++dta) {
            if((*dta).getCompNum()==1) {
                write_vtk_scalar_ascii(output, &(*dta));
            } else {
                write_vtk_vector_ascii(output, &(*dta));
            }
        }
    }
}

/**
 * \brief Write names of scalar and vector values to the VTK file (.vtu)
 * \param[in]	*out		The output file
 */
static void write_vtk_data_names(Output *output, vector<OutputData> *data)
{
    /* Write names of scalars */
    output->get_data_file() << "Scalars=\"";
    for(OutputDataVec::iterator dta = data->begin();
                dta != data->end(); ++dta) {
        if(dta->getCompNum() == 1) {
            output->get_data_file() << *dta->getName() << "_[" << *dta->getUnits() << "]";
            if((dta+1) != data->end()) {
                output->get_data_file() << ",";
            }
        }
    }
    output->get_data_file() << "\" ";

    /* Write names of vectors */
    output->get_data_file() << "Vectors=\"";
    for(OutputDataVec::iterator dta = data->begin();
                dta != data->end(); ++dta) {
        if(dta->getCompNum() == 3) {
            output->get_data_file() << *dta->getName() << "_[" << *dta->getUnits() << "]";
            if((dta+1) != data->end()) {
                output->get_data_file() << ",";
            }
        }
    }
    output->get_data_file() << "\"";
}

/**
 * \brief Write data on nodes to the VTK file (.vtu)
 * \param[in]	*out		The output file
 */
static void write_vtk_node_data(Output *output)
{
    std::vector<OutputData> *node_data = output->get_node_data();

    /* Write PointData begin */
    output->get_data_file() << "<PointData ";
    write_vtk_data_names(output, node_data);
    output->get_data_file() << ">" << endl;

    /* Write own data */
    write_vtk_data_ascii(output, node_data);

    /* Write PointData end */
    output->get_data_file() << "</PointData>" << endl;
}

/**
 * \brief Write data on elements to the VTK file (.vtu)
 * \param[in]	*output		The output object
 */
static void write_vtk_element_data(Output *output)
{
    std::vector<OutputData> *elem_data = output->get_elem_data();

    /* Write PointData begin */
    output->get_data_file() << "<CellData ";
    write_vtk_data_names(output, elem_data);
    output->get_data_file() << ">" << endl;

    /* Write own data */
    write_vtk_data_ascii(output, elem_data);

    /* Write PointData end */
    output->get_data_file() << "</CellData>" << endl;
}

/**
 * \brief Write tail of VTK file (.vtu)
 * \param[in]	*out	The output file
 */
static void write_vtk_vtu_tail(Output *output)
{
    output->get_data_file() << "</UnstructuredGrid>" << endl;
    output->get_data_file() << "</VTKFile>" << endl;
}

/**
 * \brief This function write all scalar and vector data on nodes and elements
 * to the VTK file (.vtu)
 * \param[in]	*out	The output file
 */
static void write_vtk_vtu(Output *output)
{
    Mesh *mesh = output->get_mesh();

    /* Write header */
    write_vtk_vtu_head(output);

    /* Write Piece begin */
    output->get_data_file() << "<Piece NumberOfPoints=\"" << mesh->node_vector.size() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

    /* Write VTK Geometry */
    write_vtk_geometry(output);

    /* Write VTK Topology */
    write_vtk_topology(output);

    /* Write VTK scalar and vector data on nodes to the file */
    write_vtk_node_data(output);

    /* Write VTK data on elements */
    write_vtk_element_data(output);

    /* Write Piece end */
    output->get_data_file() << "</Piece>" << endl;

    /* Write tail */
    write_vtk_vtu_tail(output);
}


/**
 * \brief This function output data to serial VTK file format
 * \param[in]   *output     The pointer at output object
 */
int write_vtk_vtu_data(Output *output)
{
    /* Serial VTK file format uses only one file */
    output->set_data_file( &output->get_base_file() );

    write_vtk_vtu(output);

    return 1;
}

/**
 * \brief This function write data to VTK (.pvd and .vtu) file format
 * for specific time
 * \param[in]   *output     The pointer at output object
 * \param[in]   time        The time from start
 * \param[in]   step        The number of steps from start
 */
int write_vtk_pvd_data(OutputTime *output, double time, int step)
{
    Mesh *mesh = output->get_mesh();
    char frame_file_name[PATH_MAX];
    ofstream *data_file = new ofstream;
    int i;

    sprintf(frame_file_name, "%s-%d.vtu", output->get_base_filename().c_str(), step);

    data_file->open(frame_file_name);
    if(data_file->is_open() == false) {
        xprintf(Err, "Could not write output to the file: %s\n", frame_file_name);
        return 0;
    } else {
        /* Set upd data file */
        output->set_data_file(data_file);

        xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, output->get_base_filename().c_str());

        /* Find first directory delimiter */
        for(i=strlen(frame_file_name); i>=0; i--) {
            if(frame_file_name[i]==DIR_DELIMITER) {
                break;
            }
        }

        /* Write reference to .vtu file of current frame to pvd file */
        if(i>0) {
            /* Strip out relative path, because vtu file is in the same directory as pvd file*/
            output->get_base_file() << "<DataSet timestep=\"" << step << "\" group=\"\" part=\"0\" file=\"" << &frame_file_name[i+1] <<"\"/>" << endl;
        } else {
            /* No path was found in string "frame_file_name" */
            output->get_base_file() << "<DataSet timestep=\"" << step << "\" group=\"\" part=\"0\" file=\"" << frame_file_name <<"\"/>" << endl;
        }

        xprintf(MsgLog, "O.K.\n");

        xprintf(MsgLog, "%s: Writing output (frame %d) file %s ... ", __func__, step, frame_file_name);

        /* Write header */
        write_vtk_vtu_head(output);

        /* Write Piece begin */
        output->get_data_file() << "<Piece NumberOfPoints=\"" << mesh->node_vector.size() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        write_vtk_geometry(output);
        /* Write VTK Topology */
        write_vtk_topology(output);

        /* Write VTK scalar and vector data on nodes */
        write_vtk_node_data(output);
        /* Write VTK scalar and vector data on elements */
        write_vtk_element_data(output);

        /* Write Piece end */
        output->get_data_file() << "</Piece>" << endl;

        /* Write tail */
        write_vtk_vtu_tail(output);

        /* Close stream for file of current frame */
        data_file->close();
        delete data_file;
        output->set_data_file(NULL);

        xprintf(MsgLog, "O.K.\n");
    }

    return 1;
}

/**
 * \brief This function writes header of VTK (.pvd) file format
 * \param[in]   *output     The pointer at output Object
 */
int write_vtk_pvd_head(OutputTime *output)
{
    xprintf(MsgLog, "%s: Writing output file (head) %s ... ", __func__, output->get_base_filename().c_str() );

    output->get_base_file() << "<?xml version=\"1.0\"?>" << endl;
    output->get_base_file() << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    output->get_base_file() << "<Collection>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

/**
 * \brief This function writes tail of VTK (.pvd) file format
 * \param[in]   *output     The pointer at output Object
 */
int write_vtk_pvd_tail(OutputTime *output)
{
    xprintf(MsgLog, "%s: Writing output file (tail) %s ... ", __func__, output->get_base_filename().c_str() );

    output->get_base_file() << "</Collection>" << endl;
    output->get_base_file() << "</VTKFile>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}
