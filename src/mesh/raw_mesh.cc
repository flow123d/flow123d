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
 * @file    raw_mesh.cc
 * @ingroup mesh
 * @brief   Mesh construction
 */


#include "mesh/raw_mesh.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/partitioning.hh"
#include "la/sparse_graph.hh"
#include "input/reader_to_storage.hh"
#include "input/type_record.hh"

namespace IT = Input::Type;


RawMesh::RawMesh()
: comm_(MPI_COMM_WORLD)
{
    istringstream is("{mesh_file=\"\", optimize_mesh=false}");
    Input::ReaderToStorage reader;
    IT::Record &in_rec = const_cast<IT::Record &>(Mesh::get_input_type());
    in_rec.finish();
    reader.read_stream(is, in_rec, Input::FileFormat::format_JSON);
    in_record_ = reader.get_root_interface<Input::Record>();
}

RawMesh::RawMesh(Mesh *mesh)
: in_record_(mesh->in_record_),
  comm_(mesh->get_comm())
{
    this->init(mesh);
}

RawMesh::RawMesh(RawMesh &other)
: in_record_(other.in_record_),
  comm_(other.comm_),
  nodes_(other.nodes_)
{}

RawMesh::~RawMesh()
{}

void RawMesh::init(Mesh *mesh)
{
    this->nodes_ = mesh->nodes_;

    unsigned int size = mesh->n_elements() + mesh->bc_mesh()->n_elements();
    element_vec_.clear();
    element_vec_.reserve(size);
    std::vector<unsigned int> node_ids(4);
    for (auto ele : mesh->elements_range()) {
        for (unsigned int n=0; n<ele->n_nodes(); n++)
            node_ids[n] = ele.node(n).idx();
        this->add_element(mesh, ele->dim(), ele->region_idx(), ele->pid(), node_ids);
    }
    for (auto ele : mesh->bc_mesh()->elements_range()) {
        for (unsigned int n=0; n<ele->n_nodes(); n++)
            node_ids[n] = ele.node(n).idx();
        this->add_element(mesh, ele->dim(), ele->region_idx(), ele->pid(), node_ids);
    }
    delete graph;

}

void RawMesh::add_element(Mesh *mesh, unsigned int dim, RegionIdx region_idx, unsigned int partition_id, std::vector<unsigned int> node_ids)
{
    element_vec_.push_back( Element() );
    Element * elem = &element_vec_.back(); //[element_vec_.size()-1];
    elem->init(dim, region_idx);
    elem->pid_ = partition_id;

	unsigned int ni=0;
	for (; ni<elem->n_nodes(); ni++) {
		elem->nodes_[ni] = mesh->node_index(node_ids[ni]);
	}
	for( ;ni < 4; ni++) elem->nodes_[ni] = undef_idx;
}

