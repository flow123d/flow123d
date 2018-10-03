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
#include "fem/dof_cell_accessor.hh"
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



//template<unsigned int dim, unsigned int spacedim> inline
//DOFHandler<dim,spacedim>::DOFHandler(Mesh & _mesh)
//: DOFHandlerBase(_mesh),
//  finite_element(0)
//{
//	object_dofs = new int**[mesh->n_elements()];
//	for (int i=0; i<mesh->n_elements(); i++)
//		object_dofs[i] = NULL;
//}
//
//
//template<unsigned int dim, unsigned int spacedim> inline
//DOFHandler<dim,spacedim>::~DOFHandler()
//{
//	for (int i=0; i<mesh->n_elements(); i++)
//		if (object_dofs[i] != NULL)
//		{
//			for (int j=0; j<mesh->element[i].dim(); j++)
//				if (object_dofs[i][j] != NULL)
//					delete[] object_dofs[i][j];
//
//			delete[] object_dofs[i];
//		}
//	delete[] object_dofs;
//}
//
//
//template<unsigned int dim, unsigned int spacedim> inline
//void DOFHandler<dim,spacedim>::distribute_dofs(FiniteElement<dim,spacedim> & fe, const unsigned int offset)
//{
//	// First check if dofs are already distributed.
//	OLD_ASSERT(finite_element == 0, "Attempt to distribute DOFs multiple times!");
//
//    unsigned int next_free_dof = offset;
//    unsigned int n_obj_dofs[dim+1];
//
//    finite_element = &fe;
//    global_dof_offset = offset;
//
//    for (unsigned int dm=0; dm <= dim; dm++)
//    {
//        n_obj_dofs[dm] = 0;
//        for (unsigned int m=0; m<dof_multiplicities.size(); m++)
//            n_obj_dofs[dm] += fe.n_object_dofs(dm, dof_multiplicities[m])*dof_multiplicities[m];
//    }
//
//    // Broadcast partition of elements to all processes.
//    LongIdx *loc_part;
//    int myp = mesh->get_part()->get_init_distr()->myp();
//    if (myp == 0)
//    {
//    	loc_part = (int*)mesh->get_part()->get_loc_part();
//    }
//    else
//    {
//    	loc_part = new LongIdx[mesh->n_elements()];
//    }
//    MPI_Bcast(loc_part, mesh->n_elements(), MPI_INT, 0, mesh->get_part()->get_init_distr()->get_comm());
//
//    // Distribute element dofs.
//    // First we distribute dofs on elements associated to process 0 and so on.
//    for (int proc=0; proc<mesh->get_part()->get_init_distr()->np(); proc++)
//    {
//    	if (proc == myp)
//    		loffset_ = next_free_dof;
//
//    	for (auto cell : mesh->elements_range()) {
//    		if (cell->dim() != dim) continue;
//
//    		if (loc_part[cell.index()] != proc) continue;
//
//    		// distribute dofs
//			// TODO: For the moment we distribute only dofs associated to the cell
//			//       In the future we want to distribute dofs on vertices, lines,
//			//       and triangles as well.
//			object_dofs[cell.index()] = new int*[dim+1];
//			for (int i=0; i<dim+1; i++)
//				object_dofs[cell.index()][i] = NULL;
//			object_dofs[cell.index()][dim] = new int[n_obj_dofs[dim]];
//
//			for (unsigned int i=0; i<n_obj_dofs[dim]; i++)
//			   object_dofs[cell.index()][dim][i] = next_free_dof++;
//    	}
//
//    	if (proc == myp)
//    		lsize_ = next_free_dof - loffset_;
//    }
//
//    // Finally we free the unused array loc_part.
//    if (mesh->get_part()->get_init_distr()->myp() != 0)
//    	delete[] loc_part;
//
////    for (auto cell : mesh->elements_range()) {
////        // skip cells of different dimension
////        if (cell->dim() != dim) continue;
////
////        // distribute dofs
////        // TODO: For the moment we distribute only dofs associated to the cell
////        //       In the future we want to distribute dofs on vertices, lines,
////        //       and triangles as well.
////        object_dofs[dim][cell] = new int[n_obj_dofs[dim]];
////        for (unsigned int i=0; i<n_obj_dofs[dim]; i++)
////           object_dofs[dim][cell][i] = next_free_dof++;
////    }
//
//
////    for (auto cell : mesh->elements_range()) {
////        // skip cells of different dimension
////        if (cell->dim != dim) continue;
////
////        // distribute dofs
////        for (int dm=0; dm<=dim; dm++)
////        {
////            for (int i=0; i<cell->n_sides_by_dim(dm); i++)
////            {
////                void *side = cell->side_by_dim(dm, i);
////                // check if node has already assigned dofs, otherwise
////                // distribute
////                if (object_dofs[dm].find(side) == object_dofs[dm].end())
////                {
////                    object_dofs[dm][side] = new int[n_obj_dofs[dm]];
////                    for (int i=0; i<n_obj_dofs[dm]; i++)
////                        object_dofs[dm][side][i] = next_free_dof++;
////                }
////            }
////        }
////    }
//
//    n_dofs = next_free_dof - offset;
//}
//
//
//
//template<unsigned int dim, unsigned int spacedim> inline
//const unsigned int DOFHandler<dim,spacedim>::n_local_dofs()
//{
//    return finite_element->n_dofs();
//}
//
//
//
//template<unsigned int dim, unsigned int spacedim>
//void DOFHandler<dim,spacedim>::get_dof_indices(const CellIterator &cell, unsigned int indices[])
//{
////    void *side;
////    unsigned int offset, pid;
//
//    for (unsigned int k=0; k<finite_element->n_object_dofs(dim,DOF_SINGLE); k++)
//    	indices[k] = object_dofs[cell.index()][dim][k];
//
////    indices.clear();
////
////    get_object_dof_indices<0>(cell, indices);
////    get_object_dof_indices<1>(cell, indices);
////    get_object_dof_indices<2>(cell, indices);
////    get_object_dof_indices<3>(cell, indices);
//}
//
////template<unsigned int dim> template<unsigned int obj_dim> inline void DOFHandler<dim>::get_object_dof_indices(const CellIterator &cell, unsigned int indices[])
////{
//    // TODO: implement for lower dimensional objects
//
////    void *side;
////    unsigned int offset, pid;
////
////    // loop over cell points/lines/triangles/tetrahedra
////    for (int i=0; i<n_simplex_objects<dim>(obj_dim); i++)
////    {
////        side   = cell->side_by_dim(obj_dim,i);
////        pid    = permutation_id<dim,obj_dim>(cell,i);
////        offset = 0;
////        // loop over dof multiplicities (single dofs, pairs, triples, sextuples)
////        for (vector<unsigned int>::iterator m=dof_multiplicities.begin(); m!=dof_multiplicities.end(); m++)
////        {
////            // loop over particular single dofs/dof pairs/triples/sextuples
////            for (int j=0; j<finite_element.n_object_dofs(obj_dim,*m); j++)
////            {
////                // loop over particular dofs (the single dof/2 dofs in the pair etc.)
////                for (int k=0; k<*m; k++)
////                    indices.push_back(object_dofs[obj_dim][side][offset+Simplex<obj_dim>::pair_permutations[pid][k]]);
////
////                offset += *m;
////            }
////        }
////    }
////}
//
//template<unsigned int dim, unsigned int spacedim> inline
//void DOFHandler<dim,spacedim>::get_dof_values(const CellIterator &cell, const Vec &values, double local_values[])
//{
//    unsigned int indices[finite_element->n_dofs()];
//
//    get_dof_indices(cell, indices);
//    VecGetValues(values, finite_element->n_dofs(), (PetscInt *)indices, local_values);
//}



DOFHandlerMultiDim::DOFHandlerMultiDim(Mesh& _mesh, bool make_elem_part)
	: DOFHandlerBase(_mesh),
	  ds_(nullptr)
{
	if (make_elem_part) make_elem_partitioning();
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
    for (auto cell : this->local_range())
    {
        cell_starts[row_4_el[cell.element_idx()]+1] = n_dofs(cell.element_accessor());
        max_elem_dofs_ = max((int)max_elem_dofs_, (int)n_dofs(cell.element_accessor()));
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
	for (auto cell : this->own_range())
    {
      for (unsigned int n=0; n<cell->dim()+1; n++)
      {
        unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.element_idx()]].nodes[n];
        node_status[nid] = VALID_NODE;
      }
    }
    
    // unmark dofs on ghost cells from lower procs
	for (auto cell : this->ghost_range())
    {
      if (cell.element_accessor().proc() < el_ds_->myp())
      {
        for (unsigned int n=0; n<cell->dim()+1; n++)
        {
          unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.element_idx()]].nodes[n];
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
    for (auto cell : this->own_range())
    {
        if (!update_cells[cell.local_idx()]) continue;
        
        // loop over element dofs
        vector<unsigned int> loc_node_dof_count(cell->n_nodes(), 0);
        for (unsigned int idof = 0; idof<n_dofs(cell.element_accessor()); ++idof)
        {
            if (cell_dof(cell.element_accessor(),idof).dim == 0)
            {
                unsigned int dof_nface_idx = cell_dof(cell.element_accessor(), idof).n_face_idx;
                if (dof_indices[cell_starts[row_4_el[cell.element_idx()]]+idof] == INVALID_DOF)
                {   // update nodal dof
                    unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.element_idx()]].nodes[dof_nface_idx];
                    dof_indices[cell_starts[row_4_el[cell.element_idx()]]+idof] = node_dofs[node_dof_starts[nid]+loc_node_dof_count[dof_nface_idx]];
                }
                loc_node_dof_count[dof_nface_idx]++;
            }
        }
    }
}


void DOFHandlerMultiDim::distribute_dofs(std::shared_ptr<DiscreteSpace> ds,
        bool sequential)
{
	// First check if dofs are already distributed.
	OLD_ASSERT(ds_ == nullptr, "Attempt to distribute DOFs multiple times!");
    
    ds_ = ds;

    std::vector<LongIdx> node_dofs, node_dof_starts;
    std::vector<short int> node_status(mesh_->tree->n_nodes(), INVALID_NODE);
    std::vector<bool> update_cells(el_ds_->lsize(), false);
    unsigned int next_free_dof = 0;

    init_cell_starts();
    init_node_dof_starts(node_dof_starts);
    node_dofs.resize(node_dof_starts[node_dof_starts.size()-1]);
    init_node_status(node_status);
    
    // Distribute dofs on local elements.
    dof_indices.resize(cell_starts[cell_starts.size()-1]);
    for (auto cell : this->own_range())
    {
      
      // loop over element dofs
      vector<unsigned int> loc_node_dof_count(cell->n_nodes(), 0);
      for (unsigned int idof = 0; idof<n_dofs(cell.element_accessor()); ++idof)
      {
        unsigned int dof_dim = cell_dof(cell.element_accessor(), idof).dim;
        unsigned int dof_nface_idx = cell_dof(cell.element_accessor(), idof).n_face_idx;
        
        if (dof_dim == 0)
        {   // add dofs shared by nodes
            unsigned int nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.element_idx()]].nodes[dof_nface_idx];
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
                update_cells[cell.local_idx()] = true;
                break;
            }
            dof_indices[cell_starts[row_4_el[cell.element_idx()]]+idof] = node_dofs[node_dof_idx];
            loc_node_dof_count[dof_nface_idx]++;
        }
        else if (dof_dim == cell->dim())
            // add dofs owned only by the element
            dof_indices[cell_starts[row_4_el[cell.element_idx()]]+idof] = next_free_dof++;
        else
            ASSERT(false).error("Unsupported dof n_face.");
      }
    }
    node_status.clear();
    
    lsize_ = next_free_dof;

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
  for (auto cell : this->own_range())
  {
    cell_starts_loc[row_4_el[cell.element_idx()]+1] = n_dofs(cell.element_accessor());
  }
  for (unsigned int i=0; i<mesh_->n_elements(); ++i)
    cell_starts_loc[i+1] += cell_starts_loc[i];

  // construct dof_indices_loc
  dof_indices_loc.resize(cell_starts_loc[mesh_->n_elements()]);
  for (auto cell : this->own_range())
  {
    for (unsigned int idof=0; idof<n_dofs(cell.element_accessor()); idof++)
        dof_indices_loc[cell_starts_loc[row_4_el[cell.element_idx()]]+idof] = dof_indices[cell_starts[row_4_el[cell.element_idx()]]+idof];
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



unsigned int DOFHandlerMultiDim::get_dof_indices(const ElementAccessor<3> &cell, std::vector<int> &indices) const
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



unsigned int DOFHandlerMultiDim::get_loc_dof_indices(const ElementAccessor<3> &cell, std::vector<LongIdx> &indices) const
{
  unsigned int ndofs = 0;
  if ( cell_starts_seq.size() > 0 && dof_indices_seq.size() > 0)
  {
    ndofs = cell_starts_seq[row_4_el[cell.idx()]+1]-cell_starts_seq[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = cell_starts_seq[row_4_el[cell.idx()]]+k;
  }
  else
  {
    ndofs = cell_starts[row_4_el[cell.idx()]+1]-cell_starts[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = cell_starts[row_4_el[cell.idx()]]+k;
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
	for (auto cell : this->own_range())
    {
      unsigned int obj_idx = mesh_->tree->obj_4_el()[cell.element_idx()];
      for (unsigned int nid=0; nid<cell->n_nodes(); nid++)
        node_is_local[mesh_->tree->objects(cell->dim())[obj_idx].nodes[nid]] = true;
    }
    for (unsigned int nid=0; nid<mesh_->tree->n_nodes(); nid++)
      if (node_is_local[nid])
        node_4_loc.push_back(nid);
    
    // create array of local ghost cells
    for ( auto cell : mesh_->elements_range() )
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


Range<DofCellAccessor, DOFHandlerMultiDim> DOFHandlerMultiDim::own_range() const {
    return Range<DofCellAccessor, DOFHandlerMultiDim>(this, 0, el_ds_->lsize());
}


Range<DofCellAccessor, DOFHandlerMultiDim> DOFHandlerMultiDim::local_range() const {
    return Range<DofCellAccessor, DOFHandlerMultiDim>(this, 0, el_ds_->lsize()+ghost_4_loc.size());
}


Range<DofCellAccessor, DOFHandlerMultiDim> DOFHandlerMultiDim::ghost_range() const {
    return Range<DofCellAccessor, DOFHandlerMultiDim>(this, el_ds_->lsize(), el_ds_->lsize()+ghost_4_loc.size());
}

