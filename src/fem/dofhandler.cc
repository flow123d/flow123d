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
 * @file
 * @brief Declaration of class which handles the ordering of degrees of freedom (dof) and mappings between local and global dofs.
 * @author Jan Stebel
 */


#include "fem/dofhandler.hh"
#include "fem/finite_element.hh"
#include "mesh/mesh.h"
#include "mesh/partitioning.hh"
#include "la/distribution.hh"







DOFHandlerMultiDim::DOFHandlerMultiDim(Mesh& _mesh)
	: DOFHandlerBase(_mesh),
	  fe1d_(0),
	  fe2d_(0),
	  fe3d_(0)
{
	object_dofs = new int**[mesh_->n_elements()];
	for (unsigned int i=0; i<mesh_->n_elements(); i++)
		object_dofs[i] = NULL;

	make_elem_partitioning();
}

void DOFHandlerMultiDim::distribute_dofs(FiniteElement<1, 3>& fe1d,
		FiniteElement<2, 3>& fe2d, FiniteElement<3, 3>& fe3d,
		const unsigned int offset)
{
	// First check if dofs are already distributed.
	ASSERT((fe1d_ == 0) && (fe2d_ == 0) && (fe3d_ == 0), "Attempt to distribute DOFs multiple times!");

    unsigned int next_free_dof = offset;
    // n_obj_dofs[element_dim][general_face_dim]
    // for every dim of ref. element, number of DoFs on its generalized face (face, edge, vertex)
    unsigned int n_obj_dofs[4][4];


    fe1d_ = &fe1d;
   	fe2d_ = &fe2d;
   	fe3d_ = &fe3d;
    global_dof_offset = offset;

    for (unsigned int dm=0; dm <= 1; dm++)
    {
        n_obj_dofs[1][dm] = 0;
        for (unsigned int m=0; m<dof_multiplicities.size(); m++)
            n_obj_dofs[1][dm] += fe1d.n_object_dofs(dm, dof_multiplicities[m])*dof_multiplicities[m];
    }
    for (unsigned int dm=0; dm <= 2; dm++)
	{
		n_obj_dofs[2][dm] = 0;
		for (unsigned int m=0; m<dof_multiplicities.size(); m++)
			n_obj_dofs[2][dm] += fe2d.n_object_dofs(dm, dof_multiplicities[m])*dof_multiplicities[m];
	}
    for (unsigned int dm=0; dm <= 3; dm++)
	{
		n_obj_dofs[3][dm] = 0;
		for (unsigned int m=0; m<dof_multiplicities.size(); m++)
			n_obj_dofs[3][dm] += fe3d.n_object_dofs(dm, dof_multiplicities[m])*dof_multiplicities[m];
	}

    // Broadcast partition of elements to all processes.
    int *loc_part;
    unsigned int myp = mesh_->get_part()->get_init_distr()->myp();
    if (myp == 0)
    {
    	loc_part = (int*)mesh_->get_part()->get_loc_part();
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

    	FOR_ELEMENTS(mesh_, cell)
		{
    		if (loc_part[cell.index()] != (int)proc) continue;

    		unsigned int dim = cell->dim();

    		// distribute dofs
			// TODO: For the moment we distribute only dofs associated to the cell
			//       In the future we want to distribute dofs on vertices, lines,
			//       and triangles as well.
			object_dofs[cell.index()] = new int*[dim+1];
			for (unsigned int i=0; i<dim+1; i++)
				object_dofs[cell.index()][i] = NULL;
			object_dofs[cell.index()][dim] = new int[n_obj_dofs[dim][dim]];

			for (unsigned int i=0; i<n_obj_dofs[dim][dim]; i++)
			   object_dofs[cell.index()][dim][i] = next_free_dof++;
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

void DOFHandlerMultiDim::get_dof_indices(const CellIterator &cell, unsigned int indices[]) const
{
	unsigned int dim = cell->dim();
	switch (dim)
	{
	case 1:
		for (unsigned int k=0; k<fe1d_->n_object_dofs(dim,DOF_SINGLE); k++)
			indices[k] = object_dofs[cell.index()][dim][k];
		break;
	case 2:
		for (unsigned int k=0; k<fe2d_->n_object_dofs(dim,DOF_SINGLE); k++)
			indices[k] = object_dofs[cell.index()][dim][k];
		break;
	case 3:
		for (unsigned int k=0; k<fe3d_->n_object_dofs(dim,DOF_SINGLE); k++)
			indices[k] = object_dofs[cell.index()][dim][k];
		break;
	}
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

    unsigned int indices[ndofs];

    get_dof_indices(cell, indices);
    VecGetValues(values, ndofs, (PetscInt *)indices, local_values);
}

DOFHandlerMultiDim::~DOFHandlerMultiDim()
{
	for (ElementFullIter elem=mesh_->element.begin(); elem!=mesh_->element.end(); ++elem)
		if (object_dofs[elem.index()] != NULL)
		{
			for (unsigned int j=0; j<elem->dim(); j++)
				if (object_dofs[elem.index()][j] != NULL)
					delete[] object_dofs[elem.index()][j];

			delete[] object_dofs[elem.index()];
		}
	delete[] object_dofs;
}



void DOFHandlerMultiDim::make_elem_partitioning()
{
	// create local arrays of elements
    int *id_4_old = new int[mesh_->n_elements()];
    int i = 0;
    FOR_ELEMENTS(mesh_, ele) id_4_old[i++] = ele.index();
    mesh_->get_part()->id_maps(mesh_->n_elements(), id_4_old, el_ds_, el_4_loc, row_4_el);
    delete[] id_4_old;

    // create local array of edges
    for (unsigned int iedg=0; iedg<mesh_->edges.size(); iedg++)
    {
        bool is_edge_local = false;
        Edge *edg = &mesh_->edges[iedg];
        for (int sid=0; sid<edg->n_sides; sid++)
        	if (el_is_local(edg->side(sid)->element().index()))
        	{
        		is_edge_local = true;
        		break;
        	}
        if (is_edge_local)
        	edg_4_loc.push_back(iedg);
    }

    // create local array of neighbours
	for (unsigned int inb=0; inb<mesh_->vb_neighbours_.size(); inb++)
	{
		Neighbour *nb = &mesh_->vb_neighbours_[inb];
		if (el_is_local(mesh_->element.index(nb->element()))
				|| el_is_local(nb->side()->element().index()))
			nb_4_loc.push_back(inb);
	}
}


bool DOFHandlerMultiDim::el_is_local(int index) const
{
	return el_ds_->is_local(row_4_el[index]);
}





template<> FiniteElement<1,3> *DOFHandlerMultiDim::fe<1>() const { return fe1d_; }
template<> FiniteElement<2,3> *DOFHandlerMultiDim::fe<2>() const { return fe2d_; }
template<> FiniteElement<3,3> *DOFHandlerMultiDim::fe<3>() const { return fe3d_; }

