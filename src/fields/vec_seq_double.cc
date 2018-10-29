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
 * @file    vec_seq_double.cc
 * @brief
 */


#include "fields/vec_seq_double.hh"
#include "fields/fe_value_handler.hh"
#include "fields/field_fe.hh"
#include "fem/fe_p.hh"
#include "fem/mapping_p1.hh"
#include "fem/fe_system.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/range_wrapper.hh"

template <int spacedim, class Value>
std::shared_ptr<FieldFE<spacedim, Value> > VectorSeqDouble::create_field(Mesh & mesh, unsigned int n_comp)
{
	static MappingP1<1,3> map1;
	static MappingP1<2,3> map2;
	static MappingP1<3,3> map3;

	if ( (dh_ == nullptr) || (dh_->mesh()->n_elements() != mesh.n_elements()) ) {
		switch (n_comp) { // by number of components
			case 1: { // scalar
				fe0_ = new FE_P_disc<0>(0);
				fe1_ = new FE_P_disc<1>(0);
				fe2_ = new FE_P_disc<2>(0);
				fe3_ = new FE_P_disc<3>(0);
				break;
			}
			case 3: { // vector
				std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
				std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
				std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
				std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
				fe0_ = new FESystem<0>(fe0_ptr, FEType::FEVector, 3);
				fe1_ = new FESystem<1>(fe1_ptr, FEType::FEVector, 3);
				fe2_ = new FESystem<2>(fe2_ptr, FEType::FEVector, 3);
				fe3_ = new FESystem<3>(fe3_ptr, FEType::FEVector, 3);
				break;
			}
			case 9: { // tensor
				std::shared_ptr< FiniteElement<0> > fe0_ptr = std::make_shared< FE_P_disc<0> >(0);
				std::shared_ptr< FiniteElement<1> > fe1_ptr = std::make_shared< FE_P_disc<1> >(0);
				std::shared_ptr< FiniteElement<2> > fe2_ptr = std::make_shared< FE_P_disc<2> >(0);
				std::shared_ptr< FiniteElement<3> > fe3_ptr = std::make_shared< FE_P_disc<3> >(0);
				fe0_ = new FESystem<0>(fe0_ptr, FEType::FETensor, 9);
				fe1_ = new FESystem<1>(fe1_ptr, FEType::FETensor, 9);
				fe2_ = new FESystem<2>(fe2_ptr, FEType::FETensor, 9);
				fe3_ = new FESystem<3>(fe3_ptr, FEType::FETensor, 9);
				break;
			}
			default:
				ASSERT(false).error("Should not happen!\n");
		}

		DOFHandlerMultiDim dh_par(mesh);
		std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( &mesh, fe0_, fe1_, fe2_, fe3_);
		dh_par.distribute_dofs(ds);
        dh_ = dh_par.sequential();
	}

	VectorSeqDouble *data_vec = new VectorSeqDouble();
	data_vec->resize( this->size() );
	std::shared_ptr< FieldFE<spacedim, Value> > field_ptr = std::make_shared< FieldFE<spacedim, Value> >();
	field_ptr->set_fe_data(dh_, &map1, &map2, &map3, data_vec);
	return field_ptr;
}

template <int spacedim, class Value>
void VectorSeqDouble::fill_output_data(std::shared_ptr<FieldFE<spacedim, Value> > field_ptr)
{
	ASSERT_EQ(field_ptr->dh_->hash(), dh_->hash());
	unsigned int ndofs = dh_->max_elem_dofs();
	unsigned int idof; // iterate over indices
	std::vector<LongIdx> indices(ndofs);

	/*for (auto cell : dh_->own_range()) {
		cell.get_dof_indices(indices);
		for(idof=0; idof<ndofs; idof++) (*field_ptr->data_vec_)[ indices[idof] ] = (*data_ptr_)[ ndofs*cell.elm_idx()+idof ];
	}*/

	for (auto ele : dh_->mesh()->elements_range()) {
		dh_->get_dof_indices(ele, indices);
		for(idof=0; idof<ndofs; idof++) (*field_ptr->data_vec_)[ indices[idof] ] = (*data_ptr_)[ ndofs*ele.idx()+idof ];
	}
}


// Instantiations of templated methods
template std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> VectorSeqDouble::create_field(Mesh & mesh, unsigned int n_comp);
template std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> VectorSeqDouble::create_field(Mesh & mesh, unsigned int n_comp);
template void VectorSeqDouble::fill_output_data(std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> field_ptr);
template void VectorSeqDouble::fill_output_data(std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> field_ptr);
