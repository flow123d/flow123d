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
	static std::vector<std::vector<std::vector<unsigned int>>> side_nodes =
	    {
		    {{0}, {1}},
		    {{0,1}, {0,2}, {1,2}},
		    {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}
	    };

    this->nodes_ = mesh->nodes_;

    unsigned int size = mesh->n_elements();
    element_vec_.clear();
    element_vec_.reserve(size);
    std::vector<unsigned int> node_ids(4);
    for (auto ele : mesh->elements_range()) {
        for (unsigned int n=0; n<ele->n_nodes(); n++)
            node_ids[n] = ele.node(n).idx();
        this->add_element(mesh, ele->dim(), ele->region_idx(), ele->pid(), node_ids);
    }
//    for (auto ele : mesh->bc_mesh()->elements_range()) {
//        for (unsigned int n=0; n<ele->n_nodes(); n++)
//            node_ids[n] = ele.node(n).idx();
//        this->add_element(mesh, ele->dim(), ele->region_idx(), ele->pid(), node_ids);
//    }

    int mesh_size = mesh->n_elements();
    init_el_ds_ = new Distribution(DistributionLocalized(), mesh_size, comm_ );  // initial distr.
    SparseGraph * graph = new SparseGraphMETIS(*init_el_ds_);

    int num_of_procs = init_el_ds_->np();
    if (mesh_size < num_of_procs) { // check if decomposing is possible
        THROW( Partitioning::ExcDecomposeMesh() << Partitioning::EI_NElems( mesh_size ) << Partitioning::EI_NProcs( num_of_procs ) );
    }

    unsigned int i_elem = 0;
    for (auto ele : element_vec_) {
        if ( !init_el_ds_->is_local( i_elem ) )
            continue;
        std::vector<unsigned int> nodes_list(ele.n_nodes()-1);
		std::vector<unsigned int> isec_elm_list;
		auto &elem_side_nodes = side_nodes[ele.dim()-1];
		std::cout << "Element: " << i_elem << ", dim: " << ele.dim() << std::endl;
		for (unsigned int i=0; i<ele.n_nodes(); ++i) std::cout << " " << ele.node_idx(i);
		    for (unsigned int j=0; j<nodes_list.size(); j++) {
		        nodes_list[j] = ele.node_idx(elem_side_nodes[i][j]);
            mesh->intersect_element_lists(nodes_list, isec_elm_list);
            for (auto insec_elem : isec_elm_list)
            	if (insec_elem != i_elem ) graph->set_edge( i_elem, insec_elem );
		}
		i_elem++;
    }

    // compute partitioning
    loc_part_ = new LongIdx[init_el_ds_->lsize()];
    graph->partition(loc_part_);
    delete graph;

}

void RawMesh::add_element(FMT_UNUSED Mesh *mesh, unsigned int dim, RegionIdx region_idx, unsigned int partition_id, std::vector<unsigned int> node_ids)
{
    element_vec_.push_back( Element() );
    Element * elem = &element_vec_.back(); //[element_vec_.size()-1];
    elem->init(dim, region_idx);
    elem->pid_ = partition_id;

	unsigned int ni=0;
	for (; ni<elem->n_nodes(); ni++) {
		elem->nodes_[ni] = node_ids[ni]; //mesh->node_index(node_ids[ni]);
	}
	for( ;ni < 4; ni++) elem->nodes_[ni] = undef_idx;
}

