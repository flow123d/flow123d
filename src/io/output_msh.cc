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
 * @file    output_msh.cc
 * @brief   The functions for outputs to GMSH files.
 *
 */

#include "io/output.h"
#include "io/output_msh.h"
#include "xio.h"
#include "mesh/mesh.h"

// TODO: v tomto souboru se poflakuji vselijake funkce mimo tridy, to v objektovem navrhu nema co delat
// maji to byt privatni metody nejake tridy (patrne Output)

/**
 * \brief This function write header of GMSH (.msh) file format
 * \param[in]   *output     The pointer at output object
 */
void write_msh_header(Output *output)
{
    // Write simple header
    output->get_base_file() << "$MeshFormat" << endl;
    output->get_base_file() << "2" << " 0 " << sizeof(double) << endl;
    output->get_base_file() << "$EndMeshFormat" << endl;
}

/**
 * \brief This function writes geometry (position of nodes) to GMSH (.msh) file
 * format
* \param[in]   *output     The pointer at output object
 */
void write_msh_geometry(Output *output)
{
    Mesh* mesh = output->get_mesh();
    NodeIter nod;
    long int id = 1;

    // Write information about nodes
    output->get_base_file() << "$Nodes" << endl;
    output->get_base_file() <<  mesh->node_vector.size() << endl;
    FOR_NODES(mesh,  nod ) {
        output->get_base_file() << id << " " << nod->getX() << " " << nod->getY() << " " << nod->getZ() << endl;
        id++;
    }
    output->get_base_file() << "$EndNodes" << endl;
}

/**
 * \brief This function writes topology (connection of nodes) to the GMSH (.msh)
 * file format
 * \param[in]   *output     The pointer at output object
 */
void write_msh_topology(Output *output)
{
    Mesh* mesh = output->get_mesh();
    ElementIter elm;
    long int id = 1;
    int i;

    // Write information about elements
    output->get_base_file() << "$Elements" << endl;
    output->get_base_file() << mesh->n_elements() << endl;
    FOR_ELEMENTS(mesh, elm) {
        // element_id element_type 3_other_tags material region partition
        output->get_base_file() << id << " " << elm->type << " 3 " << elm->mid << " " << elm->rid << " " << elm->pid;
        FOR_ELEMENT_NODES(elm, i)
            output->get_base_file() << " " << elm->node[i]->id;
        output->get_base_file() << endl;
        id++;
    }
    output->get_base_file() << "$EndElements" << endl;
}

/**
 * \brief This function writes ascii data to GMSH (.msh) output file.
 * \param[in]   *output     The pointer at Output object
 * \param[in]   *out_data   The pointer at structure storing pointer at own data.
 */
void write_msh_ascii_data(Output *output, OutputData *out_data)
{
    ofstream &file = output->get_base_file();
    long int id = 1;

    switch(out_data->type) {
    case OUT_VECTOR_INT_SCA:
        for( std::vector<int>::iterator item = ((std::vector<int>*)out_data->data)->begin();
                item != ((std::vector<int>*)out_data->data)->end();
                ++item, ++id) {
            file << id << " " << *item << endl;
        }
        break;
    case OUT_VECTOR_INT_VEC:
        for( std::vector< vector<int> >::iterator vec = ((std::vector< vector<int> >*)out_data->data)->begin();
                vec != ((std::vector< vector<int> >*)out_data->data)->end();
                ++vec, ++id)
        {
            file << id << " ";
            for (std::vector<int>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << endl;
        }
        break;
    case OUT_VECTOR_FLOAT_SCA:
        for( std::vector<float>::iterator item = ((std::vector<float>*)out_data->data)->begin();
                item != ((std::vector<float>*)out_data->data)->end();
                ++item, ++id) {
            file << id << " " << *item << endl;
        }
        break;
    case OUT_VECTOR_FLOAT_VEC:
        for( std::vector< vector<float> >::iterator vec = ((std::vector< vector<float> >*)out_data->data)->begin();
                vec != ((std::vector< vector<float> >*)out_data->data)->end();
                ++vec)
        {
            file << id << " ";
            for (std::vector<float>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << endl;
        }
        break;
    case OUT_VECTOR_DOUBLE_SCA:
        for( std::vector<double>::iterator item = ((std::vector<double>*)out_data->data)->begin();
                item != ((std::vector<double>*)out_data->data)->end();
                ++item, ++id) {
            file << id << " " << *item << endl;
        }
        break;
    case OUT_VECTOR_DOUBLE_VEC:
        for( std::vector< vector<double> >::iterator vec = ((std::vector< vector<double> >*)out_data->data)->begin();
                vec != ((std::vector< vector<double> >*)out_data->data)->end();
                ++vec, ++id)
        {
            file << id << " ";
            for (std::vector<double>::iterator item = vec->begin();
                    item != vec->end();
                    ++item) {
                file << *item << " ";
            }
            file << endl;
        }
        break;
    case OUT_ARRAY_INT_SCA:
        for(int i=0; i<out_data->num; i++, id++) {
            file << id << " " << ((int*)out_data->data)[i] << endl;
        }
        break;
    case OUT_ARRAY_FLOAT_SCA:
        for(int i=0; i<out_data->num; i++, id++) {
            file << id << " " << ((float*)out_data->data)[i] << endl;
        }
        break;
    case OUT_ARRAY_DOUBLE_SCA:
        for(int i=0; i<out_data->num; i++, id++) {
            file << id << " " << ((double*)out_data->data)[i] << endl;
        }
        break;
    default:
        xprintf(Err, "This type of data: %d is not supported by MSH file format\n", out_data->type);
        break;
    }
}

/**
 * \brief This function write all data on nodes to output file. This function
 * is used for static and dynamic data
 * \param[in]   *output     The pointer at Output object
 * \param[in]   time        The time from start
 * \param[in]   step        The number of steps from start
 */
void write_msh_node_data(Output *output, double time, int step)
{
    std::vector<OutputData> *node_data = output->get_node_data();

    if(node_data != NULL) {
        for(OutputDataVec::iterator dta = node_data->begin();
                    dta != node_data->end();
                    ++dta)
        {
            output->get_base_file() << "$NodeData" << endl;

            output->get_base_file() << "1" << endl;     // one string tag
            output->get_base_file() << "\"" << *dta->getName() << "_[" << *dta->getUnits() <<"]\"" << endl;

            output->get_base_file() << "1" << endl;     // one real tag
            output->get_base_file() << time << endl;    // first real tag = time

            output->get_base_file() << "3" << endl;     // 3 integer tags
            output->get_base_file() << step << endl;     // step number (start = 0)
            output->get_base_file() << dta->getCompNum() << endl;   // number of components
            output->get_base_file() << dta->getValueNum() << endl;  // number of values

            write_msh_ascii_data(output, &(*dta));

            output->get_base_file() << "$EndNodeData" << endl;
        }
    }
}

/**
 * \brief This function write all data on elements to output file. This
 * function is used for static and dynamic data
 * \param[in]   *output     The pointer at output object
 * \param[in]   time        The time from start
 * \param[in]   step        The number of steps from start
 */
void write_msh_elem_data(Output *output, double time, int step)
{
    std::vector<OutputData> *elem_data = output->get_elem_data();

    if(elem_data != NULL) {
        for(OutputDataVec::iterator dta = elem_data->begin();
                    dta != elem_data->end();
                    ++dta)
        {
            output->get_base_file() << "$ElementData" << endl;

            output->get_base_file() << "1" << endl;     // one string tag
            output->get_base_file() << "\"" << *dta->getName() << "_[" << *dta->getUnits() <<"]\"" << endl;

            output->get_base_file() << "1" << endl;     // one real tag
            output->get_base_file() << time << endl;    // first real tag = time

            output->get_base_file() << "3" << endl;     // 3 integer tags
            output->get_base_file() << step << endl;     // step number (start = 0)
            output->get_base_file() << dta->getCompNum() << endl;   // number of components
            output->get_base_file() << dta->getValueNum() << endl;  // number of values

            write_msh_ascii_data(output, &(*dta));

            output->get_base_file() << "$EndElementData" << endl;
        }
    }
}

/**
 * \brief Write head of .msh file format
 * \param[in]   *output     The pointer at Output object
 * \return      This function returns 1
 */
int write_msh_data(Output *output)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, output->get_base_filename().c_str());

    write_msh_header(output);

    write_msh_geometry(output);

    write_msh_topology(output);

    write_msh_node_data(output, 0.0, 0);

    write_msh_elem_data(output, 0.0, 0);

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

/**
 * \brief Write head of .msh file format
 * \param[in]   *output     The pointer at OutputTime object
 * \return      This function returns 1
 */
int write_msh_head(OutputTime *output)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, output->get_base_filename().c_str());

    write_msh_header(output);

    write_msh_geometry(output);

    write_msh_topology(output);

    xprintf(MsgLog, "O.K.\n");

    return 1;
}


/**
 * \brief Write data of .msh file format for current step
 * \param[in]   *output     The pointer at OutputTime object
 * \param[in]   time        The time from start
 * \param[in]   step        The number of steps from start
 * \return      This function returns 1
 */
int write_msh_time_data(OutputTime *output, double time, int step)
{
    xprintf(MsgLog, "%s: Writing output file %s ... ", __func__, output->get_base_filename().c_str());

    write_msh_node_data(output, time, step);

    write_msh_elem_data(output, time, step);

    // Flush stream to be sure everything is in the file now
    output->get_base_file().flush();

    xprintf(MsgLog, "O.K.\n");

    return 1;
}

/**
 * \brief Write tail of .msh file format
 *
 * It is stupid file format. It doesn't write anything special at the end of
 * the file
 *
 * \param[in]   *output     The pointer at OutputTime object
 * \return      This function returns 1
 */
int write_msh_tail(OutputTime *output)
{
    return 1;
}
