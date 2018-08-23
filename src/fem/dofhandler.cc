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
 * @file    dofhandler.cc
 * @brief   Declaration of class which handles the ordering of degrees of freedom (dof) and mappings between local and global dofs.
 * @author  Jan Stebel
 */

#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "mesh/mesh.h"
#include "mesh/duplicate_nodes.h"
#include "mesh/partitioning.hh"
#include "mesh/long_idx.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/neighbours.h"
#include "la/distribution.hh"




DOFHandlerBase::~DOFHandlerBase()
{
  if (dof_ds_ != nullptr) delete dof_ds_;
}






DOFHandlerMultiDim::DOFHandlerMultiDim(Mesh& _mesh)
	: DOFHandlerBase(_mesh),
	  ds_(nullptr)
{
	make_elem_partitioning();
}


unsigned int DOFHandlerMultiDim::n_dofs(ElementAccessor<3> cell) const
{
    switch (cell->dim()) {
        case 1:
            return fe<1>(cell)->n_dofs();
            break;
        case 2:
            return fe<2>(cell)->n_dofs();
            break;
        case 3:
            return fe<3>(cell)->n_dofs();
            break;
    }
}


const Dof &DOFHandlerMultiDim::cell_dof(ElementAccessor<3> cell, unsigned int idof) const
{
    switch (cell.dim())
    {
        case 1:
            return fe<1>(cell)->dof(idof);
            break;
        case 2:
            return fe<2>(cell)->dof(idof);
            break;
        case 3:
            return fe<3>(cell)->dof(idof);
            break;
    }
}


void DOFHandlerMultiDim::init_cell_starts()
{
    // get number of dofs per element and then set up cell_starts
    cell_starts = std::vector<LongIdx>(mesh_->n_elements()+1, 0);
    for (unsigned int loc_el=0; loc_el<el_ds_->lsize(); ++loc_el)
    {
        auto cell = mesh_->element_accessor(el_index(loc_el));
        cell_starts[row_4_el[cell.idx()]+1] = n_dofs(cell);
        max_elem_dofs_ = max((int)max_elem_dofs_, (int)n_dofs(cell));
    }
    for (unsigned int gid=0; gid<ghost_4_loc.size(); ++gid)
    {
        auto cell = mesh_->element_accessor(ghost_4_loc[gid]);
        cell_starts[row_4_el[cell.idx()]+1] = n_dofs(cell);
        max_elem_dofs_ = max((int)max_elem_dofs_, (int)n_dofs(cell));
    }
    for (unsigned int i=0; i<mesh_->n_elements(); ++i)
        cell_starts[i+1] += cell_starts[i];
}


void DOFHandlerMultiDim::init_node_dof_starts(std::vector<LongIdx> &node_dof_starts)
{
    // initialize dofs on nodes
    // We must separate dofs for dimensions because FE functions
    // may be discontinuous on nodes shared by different
    // dimensions.
    unsigned int n_node_dofs = 0;
    
    for (unsigned int nid=0; nid<mesh_->tree->n_nodes(); nid++)
    {
      node_dof_starts.push_back(n_node_dofs);
      n_node_dofs += ds_->n_node_dofs(nid);
    }
    node_dof_starts.push_back(n_node_dofs);
}


void DOFHandlerMultiDim::init_node_status(std::vector<short int> &node_status)
{
    // mark local dofs
    for (unsigned int loc_el=0; loc_el<el_ds_->lsize(); loc_el++)
    {
      CellIterator cell = mesh_->element_accessor(el_index(loc_el));

      for (unsigned int n=0; n<cell->dim()+1; n++)
      {
        unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[n];
        node_status[nid] = VALID_NODE;
      }
    }
    
    // unmark dofs on ghost cells from lower procs
    for (unsigned int gid=0; gid<ghost_4_loc.size(); gid++)
    {
      CellIterator cell = mesh_->element_accessor(ghost_4_loc[gid]);
      if (cell.proc() < el_ds_->myp())
      {
        for (unsigned int n=0; n<cell->dim()+1; n++)
        {
          unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[n];
          node_status[nid] = INVALID_NODE;
        }
      }
    }
}


void DOFHandlerMultiDim::receive_ghost_dofs(unsigned int proc, vector<LongIdx> &dofs)
{
    // send number of elements required from the other processor
    unsigned int n_elems = ghost_proc_el[proc].size();
    MPI_Send(&n_elems, 1, MPI_UNSIGNED, proc, 0, MPI_COMM_WORLD);
    
    // send indices of elements required
    MPI_Send(&(ghost_proc_el[proc][0]), n_elems, MPI_LONG_IDX, proc, 1, MPI_COMM_WORLD);
    
    // receive numbers of dofs on required elements
    vector<unsigned int> n_dofs(n_elems);
    MPI_Recv(&(n_dofs[0]), n_elems, MPI_UNSIGNED, proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // receive dofs on required elements
    unsigned int n_dofs_sum = 0;
    for (auto nd : n_dofs) n_dofs_sum += nd;
    dofs.resize(n_dofs_sum);
    MPI_Recv(&(dofs[0]), n_dofs_sum, MPI_LONG_IDX, proc, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void DOFHandlerMultiDim::send_ghost_dofs(unsigned int proc)
{
    // receive number of elements required by the other processor
    unsigned int n_elems;
    MPI_Recv(&n_elems, 1, MPI_UNSIGNED, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // receive indices of elements required
    vector<LongIdx> elems(n_elems);
    MPI_Recv(&(elems[0]), n_elems, MPI_LONG_IDX, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    // send numbers of dofs on required elements
    vector<unsigned int> n_dofs;
    for (LongIdx el : elems) n_dofs.push_back(cell_starts[row_4_el[el]+1] - cell_starts[row_4_el[el]]);
    MPI_Send(&(n_dofs[0]), n_elems, MPI_UNSIGNED, proc, 2, MPI_COMM_WORLD);
    
    // send dofs on the required elements
    vector<LongIdx> dofs;
    for (LongIdx el : elems)
        dofs.insert(dofs.end(), dof_indices.begin()+cell_starts[row_4_el[el]], dof_indices.begin()+cell_starts[row_4_el[el]+1]);
    MPI_Send(&(dofs[0]), dofs.size(), MPI_LONG_IDX, proc, 3, MPI_COMM_WORLD);
}


void DOFHandlerMultiDim::update_local_dofs(unsigned int proc,
                                           const std::vector<bool> &update_cells,
                                           const std::vector<LongIdx> &dofs, 
                                           const std::vector<LongIdx> &node_dof_starts,
                                           std::vector<LongIdx> &node_dofs
                                          )
{
    // update dof_indices on ghost cells
    unsigned int dof_offset=0;
    for (unsigned int gid=0; gid<ghost_proc_el[proc].size(); gid++)
    {
        auto cell = mesh_->element_accessor(ghost_proc_el[proc][gid]);
        
        for (unsigned dof=0; dof<n_dofs(cell); dof++)
            dof_indices[cell_starts[row_4_el[cell.idx()]]+dof] = dofs[dof_offset+dof];
        
        vector<unsigned int> loc_node_dof_count(cell->n_nodes(), 0);
        for (unsigned int idof = 0; idof<n_dofs(cell); ++idof)
        {
            if (cell_dof(cell, idof).dim == 0)
            {   // update nodal dof
                unsigned int dof_nface_idx = cell_dof(cell, idof).n_face_idx;
                unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[dof_nface_idx];
                    
                if (node_dofs[node_dof_starts[nid]+loc_node_dof_count[dof_nface_idx]] == INVALID_DOF)
                    node_dofs[node_dof_starts[nid]+loc_node_dof_count[dof_nface_idx]] = dofs[dof_offset+idof];
                
                loc_node_dof_count[dof_nface_idx]++;
            }
        }
        
        dof_offset += n_dofs(cell);
    }
    
    // update dof_indices on local elements
    for (unsigned int loc_el=0; loc_el < el_ds_->lsize(); loc_el++)
    {
        if (!update_cells[loc_el]) continue;
        
        CellIterator cell = mesh_->element_accessor(el_index(loc_el));
        
        // loop over element dofs
        vector<unsigned int> loc_node_dof_count(cell->n_nodes(), 0);
        for (unsigned int idof = 0; idof<n_dofs(cell); ++idof)
        {
            if (cell_dof(cell,idof).dim == 0)
            {
                unsigned int dof_nface_idx = cell_dof(cell, idof).n_face_idx;
                if (dof_indices[cell_starts[row_4_el[cell.idx()]]+idof] == INVALID_DOF)
                {   // update nodal dof
                    unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[dof_nface_idx];
                    dof_indices[cell_starts[row_4_el[cell.idx()]]+idof] = node_dofs[node_dof_starts[nid]+loc_node_dof_count[dof_nface_idx]];
                }
                loc_node_dof_count[dof_nface_idx]++;
            }
        }
    }
}


void DOFHandlerMultiDim::distribute_dofs(std::shared_ptr<DiscreteSpace> ds,
        bool sequential,
        const unsigned int offset)
{
	// First check if dofs are already distributed.
	OLD_ASSERT(ds_ == nullptr, "Attempt to distribute DOFs multiple times!");
    
    ds_ = ds;
    global_dof_offset = offset;

    std::vector<LongIdx> node_dofs, node_dof_starts;
    std::vector<short int> node_status(mesh_->tree->n_nodes(), INVALID_NODE);
    std::vector<bool> update_cells(el_ds_->lsize(), false);
    unsigned int next_free_dof = offset;

    init_cell_starts();
    init_node_dof_starts(node_dof_starts);
    node_dofs.resize(node_dof_starts[node_dof_starts.size()-1]);
    init_node_status(node_status);
    
    // Distribute dofs on local elements.
    dof_indices.resize(cell_starts[cell_starts.size()-1]);
    for (unsigned int loc_el=0; loc_el < el_ds_->lsize(); loc_el++)
    {
      CellIterator cell = mesh_->element_accessor(el_index(loc_el));
      
      // loop over element dofs
      vector<unsigned int> loc_node_dof_count(cell->n_nodes(), 0);
      for (unsigned int idof = 0; idof<n_dofs(cell); ++idof)
      {
        unsigned int dof_dim = cell_dof(cell, idof).dim;
        unsigned int dof_nface_idx = cell_dof(cell, idof).n_face_idx;
        
        if (dof_dim == 0)
        {   // add dofs shared by nodes
            unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[dof_nface_idx];
            unsigned int node_dof_idx = node_dof_starts[nid]+loc_node_dof_count[dof_nface_idx];
                
            switch (node_status[nid])
            {
            case VALID_NODE:
                for (int i=0; i<node_dof_starts[nid+1] - node_dof_starts[nid]; i++)
                    node_dofs[node_dof_starts[nid]+i] = next_free_dof++;
                node_status[nid] = ASSIGNED_NODE;
                break;
            case INVALID_NODE:
                node_dofs[node_dof_idx] = INVALID_DOF;
                update_cells[loc_el] = true;
                break;
            }
            dof_indices[cell_starts[row_4_el[cell.idx()]]+idof] = node_dofs[node_dof_idx];
            loc_node_dof_count[dof_nface_idx]++;
        }
        else if (dof_dim == cell.dim())
            // add dofs owned only by the element
            dof_indices[cell_starts[row_4_el[cell.idx()]]+idof] = next_free_dof++;
        else
            ASSERT(false).error("Unsupported dof n_face.");
      }
    }
    node_status.clear();
    
    lsize_ = next_free_dof - offset;

    // communicate n_dofs across all processes
    dof_ds_ = new Distribution(lsize_, PETSC_COMM_WORLD);
    n_global_dofs_ = dof_ds_->size();
    
    // shift dof indices
    loffset_ = dof_ds_->get_starts_array()[dof_ds_->myp()];
    if (loffset_ > 0)
      for (unsigned int i=0; i<dof_indices.size(); i++)
        if (dof_indices[i] != INVALID_DOF)
          dof_indices[i] += loffset_;
    
    // communicate dofs from ghost cells
    // first propagate from lower procs to higher procs and then vice versa
    for (unsigned int from_higher = 0; from_higher < 2; from_higher++)
    {
        for (unsigned int proc : ghost_proc)
        {
            if ((proc > el_ds_->myp()) == from_higher)
            { // receive dofs from master processor
                vector<LongIdx> dofs;
                receive_ghost_dofs(proc, dofs);
                
                // update dof_indices and node_dofs on ghost elements
                update_local_dofs(proc, update_cells, dofs, node_dof_starts, node_dofs);
            }
            else
                send_ghost_dofs(proc);
        }
    }
    update_cells.clear();
    node_dofs.clear();
    node_dof_starts.clear();
    
    if (sequential) create_sequential();
}


void DOFHandlerMultiDim::create_sequential()
{
  ASSERT(dof_indices.size() > 0).error("Attempt to create sequential dof handler from uninitialized or already sequential object!");
  
  // Auxiliary vectors cell_starts_loc and dof_indices_loc contain
  // only local element data (without ghost elements).
  // Then it is possible to create sequential vectors by simple reduce/gather operation.
  vector<LongIdx> cell_starts_loc(mesh_->n_elements()+1, 0);
  vector<LongIdx> dof_indices_loc;
  
  // construct cell_starts_loc
  for (unsigned int loc_el=0; loc_el<el_ds_->lsize(); ++loc_el)
  {
    auto cell = mesh_->element_accessor(el_index(loc_el));
    cell_starts_loc[row_4_el[cell.idx()]+1] = n_dofs(cell);
  }
  for (unsigned int i=0; i<mesh_->n_elements(); ++i)
    cell_starts_loc[i+1] += cell_starts_loc[i];

  // construct dof_indices_loc
  dof_indices_loc.resize(cell_starts_loc[mesh_->n_elements()]);
  for (unsigned int loc_el=0; loc_el<el_ds_->lsize(); ++loc_el)
  {
    auto cell = mesh_->element_accessor(el_index(loc_el));
    for (unsigned int idof=0; idof<n_dofs(cell); idof++)
        dof_indices_loc[cell_starts_loc[row_4_el[cell.idx()]]+idof] = dof_indices[cell_starts[row_4_el[cell.idx()]]+idof];
  }
  
  // Parallel data are no more used.
  cell_starts.clear();
  dof_indices.clear();
  
  Distribution distr(dof_indices_loc.size(), PETSC_COMM_WORLD);
  cell_starts_seq.resize(mesh_->n_elements()+1);
  dof_indices_seq.resize(distr.size());
  
  MPI_Allreduce( cell_starts_loc.data(),
                 cell_starts_seq.data(),
                 cell_starts_loc.size(),
                 MPI_LONG_IDX,
                 MPI_SUM,
                 MPI_COMM_WORLD );
  
  MPI_Allgatherv( dof_indices_loc.data(),
                  dof_indices_loc.size(),
                  MPI_LONG_IDX,
                  dof_indices_seq.data(),
                  (const int *)distr.get_lsizes_array(),
                  (const int *)distr.get_starts_array(),
                  MPI_LONG_IDX,
                  MPI_COMM_WORLD );
}



unsigned int DOFHandlerMultiDim::get_dof_indices(const CellIterator &cell, std::vector<int> &indices) const
{
  unsigned int ndofs = 0;
  if ( cell_starts_seq.size() > 0 && dof_indices_seq.size() > 0)
  {
    ndofs = cell_starts_seq[row_4_el[cell.idx()]+1]-cell_starts_seq[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = dof_indices_seq[cell_starts_seq[row_4_el[cell.idx()]]+k];
  }
  else
  {
    ndofs = cell_starts[row_4_el[cell.idx()]+1]-cell_starts[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = dof_indices[cell_starts[row_4_el[cell.idx()]]+k];
  }
  
  return ndofs;
}



unsigned int DOFHandlerMultiDim::get_loc_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const
{
  unsigned int ndofs = 0;
  if ( cell_starts_seq.size() > 0 && dof_indices_seq.size() > 0)
  {
    ndofs = cell_starts_seq[row_4_el[cell.idx()]+1]-cell_starts_seq[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = dof_indices_seq[cell_starts_seq[row_4_el[cell.idx()]]+k] - loffset_;
  }
  else
  {
    ndofs = cell_starts[row_4_el[cell.idx()]+1]-cell_starts[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = dof_indices[cell_starts[row_4_el[cell.idx()]]+k] - loffset_;
  }

  return ndofs;
}


DOFHandlerMultiDim::~DOFHandlerMultiDim()
{}



void DOFHandlerMultiDim::make_elem_partitioning()
{
	// create local arrays of elements
    el_ds_ = mesh_->get_el_ds();
    el_4_loc = mesh_->get_el_4_loc();
    row_4_el = mesh_->get_row_4_el();

    // create local array of edges
    for (unsigned int iedg=0; iedg<mesh_->n_edges(); iedg++)
    {
        bool is_edge_local = false;
        Edge *edg = &mesh_->edges[iedg];
        for (int sid=0; sid<edg->n_sides; sid++)
        	if ( el_is_local(edg->side(sid)->element().idx()) )
        	{
        		is_edge_local = true;
        		break;
        	}
        if (is_edge_local)
        	edg_4_loc.push_back(iedg);
    }

    // create local array of neighbours
	for (unsigned int inb=0; inb<mesh_->n_vb_neighbours(); inb++)
	{
		Neighbour *nb = &mesh_->vb_neighbours_[inb];
		if ( el_is_local(nb->element().idx())
				|| el_is_local(nb->side()->element().idx()) )
			nb_4_loc.push_back(inb);
	}
	
	// create array of local nodes
	std::vector<bool> node_is_local(mesh_->tree->n_nodes(), false);
    for (unsigned int loc_el=0; loc_el<el_ds_->lsize(); loc_el++)
    {
      CellIterator cell = mesh_->element_accessor(el_index(loc_el));
      unsigned int obj_idx = mesh_->tree->obj_4_el()[cell.idx()];
      for (unsigned int nid=0; nid<cell->n_nodes(); nid++)
        node_is_local[mesh_->tree->objects(cell->dim())[obj_idx].nodes[nid]] = true;
    }
    for (unsigned int nid=0; nid<mesh_->tree->n_nodes(); nid++)
      if (node_is_local[nid])
        node_4_loc.push_back(nid);
    
    // create array of local ghost cells
    for ( auto cell : mesh_->bulk_elements_range() )
    {
      if (cell.proc() != el_ds_->myp())
      {
        bool has_local_node = false;
        unsigned int obj_idx = mesh_->tree->obj_4_el()[cell.idx()];
        for (unsigned int nid=0; nid<cell->n_nodes(); nid++)
          if (node_is_local[mesh_->tree->objects(cell->dim())[obj_idx].nodes[nid]])
          {
            has_local_node = true;
            break;
          }
        if (has_local_node)
        {
            ghost_4_loc.push_back(cell.idx());
            ghost_proc.insert(cell.proc());
            ghost_proc_el[cell.proc()].push_back(cell.idx());
        }
      }
    }
}


bool DOFHandlerMultiDim::el_is_local(int index) const
{
	return el_ds_->is_local(row_4_el[index]);
}


std::size_t DOFHandlerMultiDim::hash() const {
	return this->n_global_dofs_;
}




