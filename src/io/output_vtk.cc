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
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>


void OutputVTK::write_vtk_vtu_head(void)
{
    ofstream &file = this->output->get_data_file();

    file << "<?xml version=\"1.0\"?>" << endl;
    // TODO: test endianess of platform (this would be important, when raw
    // data will be saved to the VTK file)
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    file << "<UnstructuredGrid>" << endl;
}

void OutputVTK::write_vtk_geometry(void)
{
    Mesh *mesh = this->output->get_mesh();
    ofstream &file = this->output->get_data_file();

    NodeIter node;
    int tmp;

    /* Write Points begin*/
    file << "<Points>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    tmp = 0;
    /* Set floating point precision */
    this->output->get_data_file().precision(std::numeric_limits<double>::digits10);
    FOR_NODES(mesh, node ) {
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
    //Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    Mesh *mesh = this->output->get_mesh();
    ofstream &file = this->output->get_data_file();

    Node* node;
    ElementIter ele;
    int li, tmp;

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
        file << tmp << " ";
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    /* Write type of nodes for each element */
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->type) {
        case LINE:
            file << (int)VTK_LINE << " ";
            break;
        case TRIANGLE:
            file << (int)VTK_TRIANGLE << " ";
            break;
        case TETRAHEDRON:
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
    Mesh *mesh = this->output->get_mesh();
    ofstream &file = this->output->get_data_file();

    NodeIter node;
    int tmp, li;

    /* Write Points begin*/
    file << "<Points>" << endl;
    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    /* Write own coordinates */
    tmp = 0;
    /* Set floating point precision */
    this->output->get_data_file().precision(std::numeric_limits<double>::digits10);
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
    //Mesh* mesh = (Mesh*) ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
    Mesh *mesh = this->output->get_mesh();
    ofstream &file = this->output->get_data_file();

    Node* node;
    ElementIter ele;
    int li, tmp;

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
        file << tmp << " ";
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write DataArray begin */
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    /* Write type of nodes for each element */
    FOR_ELEMENTS(mesh, ele) {
        switch(ele->type) {
        case LINE:
            file << (int)VTK_LINE << " ";
            break;
        case TRIANGLE:
            file << (int)VTK_TRIANGLE << " ";
            break;
        case TETRAHEDRON:
            file << (int)VTK_TETRA << " ";
            break;
        }
    }
    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;

    /* Write Cells end*/
    file << "</Cells>" << endl;
}

void OutputVTK::write_vtk_ascii_cont_data(OutputData *data)
{
    ofstream &file = this->output->get_data_file();

    switch(data->type) {
    case OutputData::OUT_VECTOR_INT_SCA:
        for( std::vector<int>::iterator item = ((std::vector<int>*)data->data)->begin();
                item != ((std::vector<int>*)data->data)->end();
                ++item) {
            file << *item << " ";
        }
        break;
    case OutputData::OUT_VECTOR_INT_VEC:
        for( std::vector< vector<int> >::iterator vec = ((std::vector< vector<int> >*)data->data)->begin();
                vec != ((std::vector< vector<int> >*)data->data)->end();
                ++vec) {
            for (std::vector<int>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << "  ";
        }
        break;
    case OutputData::OUT_VECTOR_FLOAT_SCA:
        file.precision(std::numeric_limits<float>::digits10);
        for( std::vector<float>::iterator item = ((std::vector<float>*)data->data)->begin();
                item != ((std::vector<float>*)data->data)->end();
                ++item) {
            file << scientific << *item << " ";
        }
        break;
    case OutputData::OUT_VECTOR_FLOAT_VEC:
        file.precision(std::numeric_limits<float>::digits10);
        for( std::vector< vector<float> >::iterator vec = ((std::vector< vector<float> >*)data->data)->begin();
                vec != ((std::vector< vector<float> >*)data->data)->end();
                ++vec) {
            for (std::vector<float>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << scientific << *item << " ";
            }
            file << "  ";
        }
        break;
    case OutputData::OUT_VECTOR_DOUBLE_SCA:
        file.precision(std::numeric_limits<double>::digits10);
        for( std::vector<double>::iterator item = ((std::vector<double>*)data->data)->begin();
                item != ((std::vector<double>*)data->data)->end();
                ++item) {
            file << scientific << *item << " ";
        }
        break;
    case OutputData::OUT_VECTOR_DOUBLE_VEC:
        file.precision(std::numeric_limits<double>::digits10);
        for( std::vector< vector<double> >::iterator vec = ((std::vector< vector<double> >*)data->data)->begin();
                vec != ((std::vector< vector<double> >*)data->data)->end();
                ++vec) {
            for (std::vector<double>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << scientific << *item << " ";
            }
            file << "  ";
        }
        break;
    case OutputData::OUT_ARRAY_INT_SCA:
        for(int i=0; i<data->num; i++) {
            file << ((int*)data->data)[i] << " ";
        }
        break;
    case OutputData::OUT_ARRAY_FLOAT_SCA:
        file.precision(std::numeric_limits<float>::digits10);
        for(int i=0; i<data->num; i++) {
            file << scientific << ((float*)data->data)[i] << " ";
        }
        break;
    case OutputData::OUT_ARRAY_DOUBLE_SCA:
        file.precision(std::numeric_limits<double>::digits10);
        for(int i=0; i<data->num; i++) {
            file << scientific << ((double*)data->data)[i] << " ";
        }
        break;
    default:
        xprintf(Err, "This type of data: %d is not supported by VTK file format\n", data->type);
        break;
    }
}

void OutputVTK::write_vtk_ascii_discont_data(OutputData *data)
{
    ofstream &file = this->output->get_data_file();
    Mesh *mesh = this->output->get_mesh();
    Node* node;
    int li, tmp = 0;

    FOR_NODES(mesh, node) {
        node->aux = tmp;   /* store index in the auxiliary variable */
        tmp++;
    }

    switch(data->type) {
    case OutputData::OUT_VECTOR_INT_SCA:
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                file << scientific << ((std::vector<int>*)data->data)->at(node->aux) << " ";
            }
        }
        break;
    case OutputData::OUT_VECTOR_INT_VEC:
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                std::vector<int> &vec = ((std::vector< vector<int> >*)data->data)->at(node->aux);
                for (std::vector<int>::iterator item = vec.begin();
                        item != vec.end();
                        ++item) {
                    file << scientific << *item << " ";
                }
                file << "  ";
            }
        }
        break;
    case OutputData::OUT_VECTOR_FLOAT_SCA:
        file.precision(std::numeric_limits<float>::digits10);
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                file << scientific << ((std::vector<float>*)data->data)->at(node->aux) << " ";
            }
        }
        break;
    case OutputData::OUT_VECTOR_FLOAT_VEC:
        file.precision(std::numeric_limits<float>::digits10);
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                std::vector<float> &vec = ((std::vector< vector<float> >*)data->data)->at(node->aux);
                for (std::vector<float>::iterator item = vec.begin();
                        item != vec.end();
                        ++item) {
                    file << scientific << *item << " ";
                }
                file << "  ";
            }
        }
        break;
    case OutputData::OUT_VECTOR_DOUBLE_SCA:
        file.precision(std::numeric_limits<double>::digits10);
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                file << scientific << ((std::vector<double>*)data->data)->at(node->aux) << " ";
            }
        }
        break;
    case OutputData::OUT_VECTOR_DOUBLE_VEC:
        file.precision(std::numeric_limits<double>::digits10);
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                std::vector<double> &vec = ((std::vector< vector<double> >*)data->data)->at(node->aux);
                for (std::vector<double>::iterator item = vec.begin();
                        item != vec.end();
                        ++item) {
                    file << scientific << *item << " ";
                }
                file << "  ";
            }
        }
        break;
    case OutputData::OUT_ARRAY_INT_SCA:
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                file << ((int*)data->data)[node->aux] << " ";
            }
        }
        break;
    case OutputData::OUT_ARRAY_FLOAT_SCA:
        file.precision(std::numeric_limits<float>::digits10);
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                file << scientific << ((float*)data->data)[node->aux] << " ";
            }
        }
        break;
    case OutputData::OUT_ARRAY_DOUBLE_SCA:
        file.precision(std::numeric_limits<double>::digits10);
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, li) {
                node = ele->node[li];
                file << scientific << ((double*)data->data)[node->aux] << " ";
            }
        }
        break;
    default:
        xprintf(Err, "This type of data: %d is not supported by VTK file format\n", data->type);
        break;
    }
}

void OutputVTK::write_vtk_ascii_data(OutputData *data)
{
    /* Use write_vtk_ascii_discont_data only in situation, when there
     * are some discontinuous data and it is necessary to write continuous
     * data too */
    if(this->output->get_corner_data() != NULL &&
            this->output->get_corner_data()->empty()==false &&
            data->ref_type == OutputData::NODE_DATA)
    {
        this->write_vtk_ascii_discont_data(data);
    } else {
        this->write_vtk_ascii_cont_data(data);
    }
}

void OutputVTK::write_vtk_scalar_ascii(OutputData *data)
{
    ofstream &file = this->output->get_data_file();

    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" Name=\"" << *data->getName() << "_[" << *data->getUnits() <<"]\" format=\"ascii\">" << endl;//, name);
    /* Write own data */
    this->write_vtk_ascii_data(data);

    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;
}

void OutputVTK::write_vtk_vector_ascii(OutputData *data)
{
    ofstream &file = this->output->get_data_file();

    /* Write DataArray begin */
    file << "<DataArray type=\"Float64\" Name=\"" << *data->getName() << "_[" << *data->getUnits() << "]\" NumberOfComponents=\"" << data->getCompNum() << "\" format=\"ascii\">" << endl;

    /* Write own data */
    this->write_vtk_ascii_data(data);

    /* Write DataArray end */
    file << endl << "</DataArray>" << endl;
}

void OutputVTK::write_vtk_data_ascii(std::vector<OutputData> *data)
{
    /* Write data on nodes or elements */
    if(data != NULL) {
        for(OutputDataVec::iterator dta = data->begin();
                dta != data->end(); ++dta) {
            if((*dta).getCompNum()==1) {
                this->write_vtk_scalar_ascii(&(*dta));
            } else {
                this->write_vtk_vector_ascii(&(*dta));
            }
        }
    }
}

void OutputVTK::write_vtk_scalar_data_names(vector<OutputData> *data)
{
    ofstream &file = this->output->get_data_file();
    int tmp = 0;

    /* Write names of scalars */
    for(OutputDataVec::iterator dta = data->begin();
                dta != data->end(); ++dta) {
        if(dta->getCompNum() == 1) {
            file << *dta->getName() << "_[" << *dta->getUnits() << "]";
            file << ",";
        }
    }
}

void OutputVTK::write_vtk_vector_data_names(vector<OutputData> *data)
{
    ofstream &file = this->output->get_data_file();
    int tmp = 0;

    /* Write names of vectors */
    for(OutputDataVec::iterator dta = data->begin();
                dta != data->end(); ++dta) {
        if(dta->getCompNum() == 3) {
            file << *dta->getName() << "_[" << *dta->getUnits() << "]";
            file << ",";
        }
    }
}

void OutputVTK::write_vtk_node_data(void)
{
    ofstream &file = this->output->get_data_file();
    std::vector<OutputData> *node_data = this->output->get_node_data();
    std::vector<OutputData> *corner_data = this->output->get_corner_data();

    if((node_data != NULL && node_data->empty()==false) ||
            (corner_data != NULL && corner_data->empty()==false)) {
        /* Write <PointData begin */
        file << "<PointData ";

        /* Write names of scalars */
        file << "Scalars=\"";
        if(node_data != NULL && node_data->empty()==false) {
            this->write_vtk_scalar_data_names(node_data);
        }
        if(corner_data != NULL && corner_data->empty()==false) {
            this->write_vtk_scalar_data_names(corner_data);
        }
        file << "\" ";

        /* Write names of vectors */
        file << "Vectors=\"";
        if(node_data != NULL && node_data->empty()==false) {
            this->write_vtk_vector_data_names(node_data);
        }
        if(corner_data != NULL && corner_data->empty()==false) {
            this->write_vtk_vector_data_names(corner_data);
        }
        file << "\"";

        /* Write right bracket of <PointData */
        file << ">" << endl;

        /* Write own data */
        if(node_data != NULL) {
            this->write_vtk_data_ascii(node_data);
        }

        if(corner_data != NULL) {
            this->write_vtk_data_ascii(corner_data);
        }

        /* Write PointData end */
        file << "</PointData>" << endl;
    }
}

void OutputVTK::write_vtk_element_data(void)
{
    ofstream &file = this->output->get_data_file();
    std::vector<OutputData> *elem_data = this->output->get_elem_data();

    if(elem_data != NULL && elem_data->empty()==false) {
        /* Write PointData begin */
        file << "<CellData ";
        /* Write names of scalars */
        file << "Scalars=\"";
        this->write_vtk_scalar_data_names(elem_data);
        file << "\" ";

        file << "Vectors=\"";
        this->write_vtk_vector_data_names(elem_data);
        file << "\"";
        file << ">" << endl;

        /* Write own data */
        this->write_vtk_data_ascii(elem_data);

        /* Write PointData end */
        file << "</CellData>" << endl;
    }
}

void OutputVTK::write_vtk_vtu_tail(void)
{
    ofstream &file = this->output->get_data_file();

    file << "</UnstructuredGrid>" << endl;
    file << "</VTKFile>" << endl;
}

void OutputVTK::write_vtk_vtu(void)
{
    ofstream &file = this->output->get_data_file();
    Mesh *mesh = this->output->get_mesh();

    /* Write header */
    this->write_vtk_vtu_head();

    /* When there is no discontinuous data, then write classical vtu */
    if(this->output->get_corner_data() != NULL &&
            this->output->get_corner_data()->empty()==true)
    {
        /* Write Piece begin */
        file << "<Piece NumberOfPoints=\"" << mesh->node_vector.size() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

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
        file << "<Piece NumberOfPoints=\"" << this->output->get_corner_count() << "\" NumberOfCells=\"" << mesh->n_elements() <<"\">" << endl;

        /* Write VTK Geometry */
        this->write_vtk_discont_geometry();

        /* Write VTK Topology */
        this->write_vtk_discont_topology();

        /* Write VTK scalar and vector data on nodes to the file */
        this->write_vtk_node_data();

        /* Write VTK scalar and vector data in corners to the file */
        //this->write_vtk_corner_data();

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
    /* Serial VTK file format uses only one file */
    this->output->set_data_file( &this->output->get_base_file() );

    write_vtk_vtu();

    return 1;
}

int OutputVTK::write_data(double time)
{
    Mesh *mesh = this->output->get_mesh();
    char base_dir_name[PATH_MAX];
    char new_dir_name[PATH_MAX];
    char base_file_name[PATH_MAX];
    char base[PATH_MAX];
    char frame_file_name[PATH_MAX];
    ofstream *data_file = new ofstream;
    DIR *dir;
    int i, j, ret;
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    strncpy(base_dir_name, this->output_time->get_base_filename().c_str(), PATH_MAX);

    /* Remove last file name from base_name and find position of last directory
     * delimiter: '/' */
    for(j=strlen(base_dir_name)-1; j>0; j--) {
        if(base_dir_name[j]=='/') {
            base_dir_name[j]='\0';
            break;
        }
    }

    strncpy(base_file_name, this->output_time->get_base_filename().c_str(), PATH_MAX);

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

    sprintf(frame_file_name, "%s/%s-%06d.vtu", new_dir_name, base, this->output_time->current_step);

    data_file->open(frame_file_name);
    if(data_file->is_open() == false) {
        xprintf(Err, "Could not write output to the file: %s\n", frame_file_name);
        return 0;
    } else {
        /* Set up data file */
        this->output_time->set_data_file(data_file);

        xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, this->output_time->get_base_filename().c_str());

        /* Find first directory delimiter */
        for(i=strlen(frame_file_name); i>=0; i--) {
            if(frame_file_name[i]==DIR_DELIMITER) {
                break;
            }
        }

        /* Set floating point precision to max */
        this->output_time->get_base_file().precision(std::numeric_limits<double>::digits10);

        /* Strip out relative path and add "base/" string */
        this->output_time->get_base_file() << scientific << "<DataSet timestep=\"" << time << "\" group=\"\" part=\"0\" file=\"" << base << "/" << &frame_file_name[i+1] <<"\"/>" << endl;

        xprintf(MsgLog, "O.K.\n");

        xprintf(MsgLog, "%s: Writing output (frame %d) file %s ... ", __func__,
                this->output_time->current_step, frame_file_name);

        this->write_vtk_vtu();

        /* Close stream for file of current frame */
        data_file->close();
        delete data_file;
        this->output->set_data_file(NULL);

        xprintf(MsgLog, "O.K.\n");
    }

    return 1;
}

int OutputVTK::write_head(void)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    xprintf(MsgLog, "%s: Writing output file (head) %s ... ", __func__,
            this->output->get_base_filename().c_str() );

    this->output->get_base_file() << "<?xml version=\"1.0\"?>" << endl;
    this->output->get_base_file() << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    this->output->get_base_file() << "<Collection>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

int OutputVTK::write_tail(void)
{
    int rank=0;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    /* It's possible now to do output to the file only in the first process */
    if(rank!=0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return 0;
    }

    xprintf(MsgLog, "%s: Writing output file (tail) %s ... ", __func__,
            this->output->get_base_filename().c_str() );

    this->output->get_base_file() << "</Collection>" << endl;
    this->output->get_base_file() << "</VTKFile>" << endl;

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

OutputVTK::OutputVTK(Output *_output)
{
    this->output = _output;
    Mesh *mesh = this->output->get_mesh();
    unsigned int li;

    this->write_head();
}

OutputVTK::OutputVTK(OutputTime *_output_time, const Input::Record &in_rec)
{
    this->output = _output_time;
    this->output_time = _output_time;

    this->write_head();
}

OutputVTK::OutputVTK(OutputTime *_output_time)
{
    this->output = _output_time;
    this->output_time = _output_time;

    this->write_head();
}

OutputVTK::~OutputVTK()
{
    this->write_tail();
}

Input::Type::Record & OutputVTK::get_input_type()
{
	using namespace Input::Type;
	static Record rec("vtk", "Parameters of vtk output format.");

	if (!rec.is_finished()) {

		// It is derived from abstract class
		rec.derive_from(OutputFormat::get_input_type());

		// The variant
		static Selection variant_sel("VTK variant (ascii or binary)");
	    variant_sel.add_value(OutputVTK::VARIANT_ASCII, "ascii",
	    		"ASCII variant of VTK file format");
	    variant_sel.add_value(OutputVTK::VARIANT_BINARY, "binary",
	    		"Binary variant of VTK file format (not supported yet)");
	    variant_sel.finish();

		rec.declare_key("variant", variant_sel, Default("ascii"),
				"Variant of output stream file format.");

		// The parallel or serial variant
		rec.declare_key("parallel", Bool(), Default("false"),
				"Parallel or serial version of file format.");

		// The compression
		static Selection compression_sel("Type of compression of VTK file format");
		compression_sel.add_value(OutputVTK::COMPRESSION_NONE, "none",
				"Data in VTK file format are not compressed");
		compression_sel.add_value(OutputVTK::COMPRESSION_GZIP, "zlib",
				"Data in VTK file format are compressed using zlib (not supported yet)");
		compression_sel.finish();

		rec.declare_key("compression", compression_sel, Default("none"),
				"Compression used in output stream file format.");

		rec.finish();
	}

	return rec;
}
