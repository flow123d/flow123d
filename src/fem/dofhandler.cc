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
#include "mesh/partitioning.hh"
#include "mesh/long_idx.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "mesh/neighbours.h"
#include "la/distribution.hh"







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
	  fe1d_(0),
	  fe2d_(0),
	  fe3d_(0)
{
	object_dofs = new LongIdx**[mesh_->n_elements()];
	for (unsigned int i=0; i<mesh_->n_elements(); i++)
		object_dofs[i] = NULL;

	if (make_elem_part) make_elem_partitioning();
}

void DOFHandlerMultiDim::distribute_dofs(FiniteElement<1>& fe1d,
		FiniteElement<2>& fe2d, FiniteElement<3>& fe3d,
		const unsigned int offset)
{
	// First check if dofs are already distributed.
	OLD_ASSERT((fe1d_ == 0) && (fe2d_ == 0) && (fe3d_ == 0), "Attempt to distribute DOFs multiple times!");

    unsigned int next_free_dof = offset;
    // n_obj_dofs[element_dim][general_face_dim]
    // for every dim of ref. element, number of DoFs on its generalized face (face, edge, vertex)
    unsigned int n_obj_dofs[4][4];


    fe1d_ = &fe1d;
   	fe2d_ = &fe2d;
   	fe3d_ = &fe3d;
    global_dof_offset = offset;
    max_elem_dofs_ = max(fe1d_->n_dofs(), max(fe2d_->n_dofs(), fe3d_->n_dofs()));

    for (unsigned int dm=0; dm <= 3; dm++)
        for (unsigned int dd=1; dd <= 3; dd++)
            n_obj_dofs[dd][dm] = 0;
    for (unsigned int d=0; d<fe1d.n_dofs(); d++)
        n_obj_dofs[1][fe1d.dof(d).dim]++;
    for (unsigned int d=0; d<fe2d.n_dofs(); d++)
        n_obj_dofs[2][fe2d.dof(d).dim]++;
    for (unsigned int d=0; d<fe3d.n_dofs(); d++)
        n_obj_dofs[3][fe3d.dof(d).dim]++;

    // Broadcast partition of elements to all processes.
    LongIdx *loc_part;
    unsigned int myp = mesh_->get_part()->get_init_distr()->myp();
    if (myp == 0)
    {
    	loc_part = (LongIdx*)mesh_->get_local_part();
    }
    else
    {
    	loc_part = new int[mesh_->n_elements()];
    }
    MPI_Bcast(loc_part, mesh_->n_elements(), MPI_INT, 0, mesh_->get_part()->get_init_distr()->get_comm());

    // Distribute element dofs.
    // First we distribute dofs on elements associated to process 0 and so on.
    for (unsigned int proc=0; proc<mesh_->get_part()->get_init_distr()->np(); proc++)
    {
    	if (proc == myp)
    		loffset_ = next_free_dof;

    	for (auto cell : mesh_->elements_range()) {
    		if (loc_part[ cell.idx() ] != (int)proc) continue;

    		unsigned int dim = cell->dim();
    		unsigned int dof_index = cell.idx();

    		// distribute dofs
			// TODO: For the moment we distribute only dofs associated to the cell
			//       In the future we want to distribute dofs on vertices, lines,
			//       and triangles as well.
			object_dofs[ dof_index ] = new LongIdx*[dim+1];
			for (unsigned int i=0; i<dim+1; i++)
				object_dofs[ dof_index ][i] = NULL;
			object_dofs[ dof_index ][dim] = new int[n_obj_dofs[dim][dim]];

			for (unsigned int i=0; i<n_obj_dofs[dim][dim]; i++)
			   object_dofs[ dof_index ][dim][i] = next_free_dof++;
    	}

    	if (proc == myp) {
    		lsize_ = next_free_dof - loffset_;
    		ds_ = new Distribution(lsize_, PETSC_COMM_WORLD);
    	}
    }

    // Finally we free the unused array loc_part.
    if (mesh_->get_part()->get_init_distr()->myp() != 0)
    	delete[] loc_part;

    n_dofs = next_free_dof - offset;
}

unsigned int DOFHandlerMultiDim::get_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const
{
	unsigned int dim = cell->dim();
    unsigned int n_objects_dofs;
    
	switch (dim)
	{
	case 1:
		n_objects_dofs = fe1d_->n_dofs();
		break;
	case 2:
		n_objects_dofs = fe2d_->n_dofs();
		break;
	case 3:
		n_objects_dofs = fe3d_->n_dofs();
		break;
	}
	
	ASSERT_LE(n_objects_dofs, indices.size()).error();

	for (unsigned int k = 0; k < n_objects_dofs; k++)
        indices[k] = object_dofs[cell.idx()][dim][k];

	return n_objects_dofs;
}

unsigned int DOFHandlerMultiDim::get_loc_dof_indices(const CellIterator &cell, std::vector<LongIdx> &indices) const
{
	/*
	 * TODO: This method is currently wrong (the shift by loc_offset_ is not enough for
	 * local <-> global mapping. Only local to global mapping should be used in future.
	 */
    unsigned int dim = cell->dim();
    unsigned int n_objects_dofs;
    
    switch (dim)
    {
    case 1:
        n_objects_dofs = fe1d_->n_dofs();
        break;
    case 2:
        n_objects_dofs = fe2d_->n_dofs();
        break;
    case 3:
        n_objects_dofs = fe3d_->n_dofs();
        break;
    }

    ASSERT_LE(n_objects_dofs, indices.size()).error();

    for (unsigned int k = 0; k < n_objects_dofs; k++)
        indices[k] = object_dofs[cell.idx()][dim][k] - loffset_;

    return n_objects_dofs;
}

void DOFHandlerMultiDim::get_dof_values(const CellIterator &cell, const Vec &values, double local_values[]) const
{
	int ndofs=0;

	switch (cell->dim())
	{
	case 1:
		ndofs = fe1d_->n_dofs();
		break;
	case 2:
		ndofs = fe2d_->n_dofs();
		break;
	case 3:
		ndofs = fe3d_->n_dofs();
		break;
	}

    std::vector<LongIdx> indices(ndofs);

    get_dof_indices(cell, indices);
    VecGetValues(values, ndofs, (PetscInt *)indices[0], local_values);
}

DOFHandlerMultiDim::~DOFHandlerMultiDim()
{
	for ( auto ele : mesh_->elements_range() ) {
		if (object_dofs[ele.idx()] != NULL)
		{
			for (unsigned int j=0; j<=ele->dim(); j++)
				if (object_dofs[ele.idx()][j] != NULL)
					delete[] object_dofs[ele.idx()][j];

			delete[] object_dofs[ele.idx()];
		}
	}
	delete[] object_dofs;
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
}


bool DOFHandlerMultiDim::el_is_local(int index) const
{
	return el_ds_->is_local(row_4_el[index]);
}


std::size_t DOFHandlerMultiDim::hash() const {
	return this->n_dofs;
}





template<> FiniteElement<1> *DOFHandlerMultiDim::fe<1>() const { return fe1d_; }
template<> FiniteElement<2> *DOFHandlerMultiDim::fe<2>() const { return fe2d_; }
template<> FiniteElement<3> *DOFHandlerMultiDim::fe<3>() const { return fe3d_; }


//template class DOFHandler<0,3>;
//template class DOFHandler<1,3>;
//template class DOFHandler<2,3>;
//template class DOFHandler<3,3>;
