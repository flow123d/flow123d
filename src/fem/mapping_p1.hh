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
 * @brief Class MappingP1 implements the affine transformation of
 *        the unit cell onto the actual cell.
 * @author Jan Stebel
 */

#ifndef MAPPING_P1_HH_
#define MAPPING_P1_HH_

#include "fem/mapping.hh"




/**
 * Class MappingP1 implements the affine transformation of
 * the reference cell onto the actual cell.
 */
template<unsigned int dim, unsigned int spacedim>
class MappingP1 : public Mapping<dim,spacedim>
{
public:
    /**
     * Constructor.
     */
    MappingP1();

    MappingInternalData *initialize(const Quadrature<dim> &q, UpdateFlags flags);

    UpdateFlags update_each(UpdateFlags flags);

    void fill_fe_values(const typename DOFHandler<dim>::CellIterator &cell,
                            const Quadrature<dim> &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data);

    void fill_fe_side_values(const typename DOFHandler<dim>::CellIterator &cell,
                            const Side &side,
                            const Quadrature<dim> &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data);


private:
    /**
     * Auxiliary matrix of gradients of shape functions (used for
     * computation of the Jacobian).
     */
    mat::fixed<dim,dim+1> grad;

};








#endif /* MAPPING_P1_HH_ */
