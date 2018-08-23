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
//    	for (auto cell : mesh->bulk_elements_range()) {
//    		if (cell->dim() != dim) continue;
//
//    		if (loc_part[cell.idx()] != proc) continue;
//
//    		// distribute dofs
//			// TODO: For the moment we distribute only dofs associated to the cell
//			//       In the future we want to distribute dofs on vertices, lines,
//			//       and triangles as well.
//			object_dofs[cell.idx()] = new int*[dim+1];
//			for (int i=0; i<dim+1; i++)
//				object_dofs[cell.idx()][i] = NULL;
//			object_dofs[cell.idx()][dim] = new int[n_obj_dofs[dim]];
//
//			for (unsigned int i=0; i<n_obj_dofs[dim]; i++)
//			   object_dofs[cell.idx()][dim][i] = next_free_dof++;
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
////    for (auto cell : mesh->bulk_elements_range()) {
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
////    for (auto cell : mesh->bulk_elements_range()) {
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
//    	indices[k] = object_dofs[cell.idx()][dim][k];
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



DOFHandlerMultiDim::DOFHandlerMultiDim(Mesh& _mesh)
	: DOFHandlerBase(_mesh),
	  ds_(nullptr)
{
//	object_dofs = new LongIdx**[mesh_->n_elements()];
//	for (unsigned int i=0; i<mesh_->n_elements(); i++)
//		object_dofs[i] = NULL;

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


void DOFHandlerMultiDim::distribute_dofs(std::shared_ptr<DiscreteSpace> ds,
		bool sequential,
		const unsigned int offset)
{
	// First check if dofs are already distributed.
	OLD_ASSERT(ds_ == nullptr, "Attempt to distribute DOFs multiple times!");
    
    ds_ = ds;
    global_dof_offset = offset;

    const int INVALID_NODE = -1;
    const int VALID_NODE = -2;
    const int ASSIGNED_NODE = -3;
    unsigned int next_free_dof = offset;
    unsigned int n_node_dofs = 0;
    unsigned int n_dof_indices = 0;

    // initialize dofs on nodes
    // We must separate dofs for dimensions because FE functions
    // may be discontinuous on nodes shared by different
    // dimensions.
    std::vector<int> node_dofs, node_dof_starts, node_status;
    for (unsigned int nid=0; nid<mesh_->tree->n_nodes(); nid++)
    {
      node_dof_starts.push_back(n_node_dofs);
      n_node_dofs += ds->n_node_dofs(nid);
    }
    node_dof_starts.push_back(n_node_dofs);
    node_dofs.resize(n_node_dofs);
    node_status.resize(mesh_->tree->n_nodes(), INVALID_NODE);
    
    // mark local dofs
    for (unsigned int loc_el=0; loc_el<el_ds_->lsize(); loc_el++)
    {
      CellIterator cell = mesh_->element_accessor(el_index(loc_el));

      n_dof_indices += ds->n_elem_dofs(cell);
      
      unsigned int nid;
      for (unsigned int n=0; n<cell->dim()+1; n++)
      {
        nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[n];
        node_status[nid] = VALID_NODE;
        n_dof_indices += ds->n_node_dofs(nid);
      }
    }
    
    // unmark dofs on ghost cells from lower procs
    for (unsigned int gid=0; gid<ghost_4_loc.size(); gid++)
    {
      CellIterator cell = mesh_->element_accessor(ghost_4_loc[gid]);
      if (mesh_->get_el_ds()->get_proc(mesh_->get_row_4_el()[cell.idx()]) < el_ds_->myp())
      {
        unsigned int nid;
        for (unsigned int n=0; n<cell->dim()+1; n++)
        {
          nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[n];
          node_status[nid] = INVALID_NODE;
        }
      }
    }
    
    // Distribute element dofs for local process.
    cell_starts = std::vector<LongIdx>(mesh_->n_elements()+1, 0);
    dof_indices.reserve(n_dof_indices);
    for (unsigned int loc_el=0; loc_el < el_ds_->lsize(); loc_el++)
    {
      CellIterator cell = mesh_->element_accessor(el_index(loc_el));

      unsigned int ncell_dofs = 0;
      switch (cell->dim()) {
        case 1:
            ncell_dofs = fe<1>(cell)->n_dofs();
            break;
        case 2:
            ncell_dofs = fe<2>(cell)->n_dofs();
            break;
        case 3:
            ncell_dofs = fe<3>(cell)->n_dofs();
            break;
      }
      
      vector<unsigned int> loc_node_dof_count(cell->n_nodes(), 0);
      for (unsigned int idof = 0; idof<ncell_dofs; ++idof)
      {
        unsigned int dof_dim = 0, dof_nface_idx = 0;
        switch (cell.dim()) {
        case 1:
            dof_dim = fe<1>(cell)->dof(idof).dim;
            dof_nface_idx = fe<1>(cell)->dof(idof).n_face_idx;
            break;
        case 2:
            dof_dim = fe<2>(cell)->dof(idof).dim;
            dof_nface_idx = fe<2>(cell)->dof(idof).n_face_idx;
            break;
        case 3:
            dof_dim = fe<3>(cell)->dof(idof).dim;
            dof_nface_idx = fe<3>(cell)->dof(idof).n_face_idx;
            break;
        }
        
        if (dof_dim == 0) {
            // add dofs shared by nodes
            unsigned int nid;
            nid = mesh_->tree->objects(cell->dim())[mesh_->tree->obj_4_el()[cell.idx()]].nodes[dof_nface_idx];
                
            switch (node_status[nid]) {
            case VALID_NODE:
                for (int i=0; i<node_dof_starts[nid+1] - node_dof_starts[nid]; i++)
                    node_dofs[node_dof_starts[nid]+i] = next_free_dof++;
                node_status[nid] = ASSIGNED_NODE;
                break;
            case INVALID_NODE:
                for (int i=0; i<node_dof_starts[nid+1] - node_dof_starts[nid]; i++)
                    node_dofs[node_dof_starts[nid]+i] = 0;
                break;
            }
            dof_indices.push_back(node_dofs[node_dof_starts[nid]+(loc_node_dof_count[dof_nface_idx]++)]);
        } else if (dof_dim == cell.dim()) {
            // add dofs owned only by the element
            dof_indices.push_back(next_free_dof++);
        } else {
            ASSERT(false).error("Unsupported dof n_face.");
        }
      }
      cell_starts[row_4_el[el_4_loc[loc_el]]+1] = dof_indices.size();
      
      max_elem_dofs_ = max((int)max_elem_dofs_, (int)dof_indices.size() - cell_starts[row_4_el[el_4_loc[loc_el]]]);
    }
    for (unsigned int i=1; i<cell_starts.size(); i++)
      if (cell_starts[i] < cell_starts[i-1])
        cell_starts[i] = cell_starts[i-1];
    
    lsize_ = next_free_dof - offset;

    // communicate n_dofs across all processes
    dof_ds_ = new Distribution(lsize_, PETSC_COMM_WORLD);
    n_global_dofs_ = dof_ds_->size();
    
    // shift dof indices
    loffset_ = dof_ds_->get_starts_array()[dof_ds_->myp()];
    if (dof_ds_->myp() > 0)
      for (unsigned int i=0; i<dof_indices.size(); i++)
        if (dof_indices[i] != INVALID_NODE)
          dof_indices[i] += loffset_;
    
    // communicate dofs from ghost cells
    // - determine ids of neighbouring processors
    // - send number of elements required from the other processor
    // - receive number of elements required by the other processor
    // - send indices of elements required
    // - receive indices of elements required
    // - receive dofs on the required elements
    // - send dofs on the required elements
    if (sequential) create_sequential();
}


void DOFHandlerMultiDim::create_sequential()
{
  // on each processor create sequential vectors of all dofs on the whole mesh.
  Distribution distr(dof_indices.size(), PETSC_COMM_WORLD);
  cell_starts_seq.resize(cell_starts.size());
  dof_indices_seq.resize(distr.size());
  
  MPI_Allreduce( &cell_starts[0], &cell_starts_seq[0], cell_starts.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allgatherv( &dof_indices[0], dof_indices.size(), MPI_INT, &dof_indices_seq[0], (const int *)distr.get_lsizes_array(), (const int *)distr.get_starts_array(), MPI_INT, MPI_COMM_WORLD);
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
    OLD_ASSERT(el_is_local(cell.idx()), "Attempt to use DOF handler on nonlocal element!");
    ndofs = cell_starts_seq[row_4_el[cell.idx()]+1]-cell_starts_seq[row_4_el[cell.idx()]];
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
    OLD_ASSERT(el_is_local(cell.idx()), "Attempt to use DOF handler on nonlocal element!");
    ndofs = cell_starts_seq[row_4_el[cell.idx()]+1]-cell_starts_seq[row_4_el[cell.idx()]];
    for (unsigned int k=0; k<ndofs; k++)
      indices[k] = dof_indices[cell_starts[row_4_el[cell.idx()]]+k] - loffset_;
  }

  return ndofs;
}

// void DOFHandlerMultiDim::get_dof_values(const CellIterator &cell, const Vec &values, double local_values[]) const
// {
// 	int ndofs=0;
// 
// 	switch (cell->dim())
// 	{
// 	case 1:
// 		ndofs = fe1d_->n_dofs();
// 		break;
// 	case 2:
// 		ndofs = fe2d_->n_dofs();
// 		break;
// 	case 3:
// 		ndofs = fe3d_->n_dofs();
// 		break;
// 	}
// 
//     std::vector<LongIdx> indices(ndofs);
// 
//     get_dof_indices(cell, indices);
//     VecGetValues(values, ndofs, (PetscInt *)indices, local_values);
// }

DOFHandlerMultiDim::~DOFHandlerMultiDim()
{
// 	for ( auto ele : mesh_->bulk_elements_range() )
// 		if (object_dofs[elem.idx()] != NULL)
// 		{
// 			for (unsigned int j=0; j<elem->dim(); j++)
// 				if (object_dofs[elem.idx()][j] != NULL)
// 					delete[] object_dofs[elem.idx()][j];
// 
// 			delete[] object_dofs[elem.idx()];
// 		}
// 	delete[] object_dofs;
}



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
      if (mesh_->get_el_ds()->get_proc(mesh_->get_row_4_el()[cell.idx()]) != el_ds_->myp())
      {
        bool has_local_node = false;
        unsigned int obj_idx = mesh_->tree->obj_4_el()[cell.idx()];
        for (unsigned int nid=0; nid<cell->n_nodes(); nid++)
          if (node_is_local[mesh_->tree->objects(cell->dim())[obj_idx].nodes[nid]])
          {
            has_local_node = true;
            break;
          }
        if (has_local_node) ghost_4_loc.push_back(cell.idx());
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





// template<> FiniteElement<1> *DOFHandlerMultiDim::fe<1>() const { return fe1d_; }
// template<> FiniteElement<2> *DOFHandlerMultiDim::fe<2>() const { return fe2d_; }
// template<> FiniteElement<3> *DOFHandlerMultiDim::fe<3>() const { return fe3d_; }


//template class DOFHandler<0,3>;
//template class DOFHandler<1,3>;
//template class DOFHandler<2,3>;
//template class DOFHandler<3,3>;
