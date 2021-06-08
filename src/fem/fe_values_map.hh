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
 * @file    fe_values_map.hh
 * @author  Jan Stebel, David Flanderka
 */

#ifndef FE_VALUES_MAP_HH_
#define FE_VALUES_MAP_HH_

#include <vector>                             // for vector
#include "fem/fe_values.hh"                   // for FEValues
#include "fem/element_values.hh"              // for ElementValues


/**
 * @brief Helper class allows update values and gradients of FEValues of FEScalar type.
 */
template<unsigned int spacedim = 3>
class MapScalar {
public:
    /// Empty method.
    inline void fill_values_vec(FMT_UNUSED FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            FMT_UNUSED const typename FEValues<spacedim>::FEInternalData &fe_data) {}

    /// Update shape_values of given FEValues object.
    inline void update_values(FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEScalar);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
                fe_values.shape_values[i][j] = fe_data.ref_shape_values[i][j][0];
    }

    /// Update shape_gradients of given FEValues object.
    inline void update_gradients(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEScalar);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
                fe_values.shape_gradients[i][j] = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j];
    }
};


/**
 * @brief Helper class allows update values and gradients of FEValues of FEVectorPiola type
 */
template<unsigned int spacedim = 3>
class MapPiola {
public:
    /// Empty method.
    inline void fill_values_vec(FMT_UNUSED FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            FMT_UNUSED const typename FEValues<spacedim>::FEInternalData &fe_data) {}

    /// Update shape_values of given FEValues object.
    inline void update_values(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEVectorPiola);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = elm_values.jacobian(i)*fe_data.ref_shape_values[i][j]/elm_values.determinant(i);
                for (unsigned int c=0; c<spacedim; c++)
                    fe_values.shape_values[i][j*spacedim+c] = fv_vec(c);
            }
    }

    /// Update shape_gradients of given FEValues object.
    inline void update_gradients(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEVectorPiola);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j] * trans(elm_values.jacobian(i))
                        / elm_values.determinant(i);
                for (unsigned int c=0; c<spacedim; c++)
                    fe_values.shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
};


/**
 * @brief Helper class allows update values and gradients of FEValues of FEVectorContravariant type
 */
template<unsigned int spacedim = 3>
class MapContravariant {
public:
    /// Empty method.
    inline void fill_values_vec(FMT_UNUSED FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            FMT_UNUSED const typename FEValues<spacedim>::FEInternalData &fe_data) {}

    /// Update shape_values of given FEValues object.
    inline void update_values(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEVectorContravariant);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = elm_values.jacobian(i) * fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    fe_values.shape_values[i][j*spacedim+c] = fv_vec[c];
            }
    }

    /// Update shape_gradients of given FEValues object.
    inline void update_gradients(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEVectorContravariant);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j] * trans(elm_values.jacobian(i));
                for (unsigned int c=0; c<spacedim; c++)
                    fe_values.shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
};


/**
 * @brief Helper class allows update values and gradients of FEValues of FEVector type
 */
template<unsigned int spacedim = 3>
class MapVector {
public:
    /// Empty method.
    inline void fill_values_vec(FMT_UNUSED FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            FMT_UNUSED const typename FEValues<spacedim>::FEInternalData &fe_data) {}

    /// Update shape_values of given FEValues object.
    inline void update_values(FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEVector);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    fe_values.shape_values[i][j*spacedim+c] = fv_vec[c];
            }
    }

    /// Update shape_gradients of given FEValues object.
    inline void update_gradients(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEVector);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j];
                for (unsigned int c=0; c<spacedim; c++)
                    fe_values.shape_gradients[i][j*spacedim+c] = grads.col(c);
            }
    }
};


/**
 * @brief Helper class allows update values and gradients of FEValues of FETensor type
 */
template<unsigned int spacedim = 3>
class MapTensor {
public:
    /// Empty method.
    inline void fill_values_vec(FMT_UNUSED FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            FMT_UNUSED const typename FEValues<spacedim>::FEInternalData &fe_data) {}

    /// Update shape_values of given FEValues object.
    inline void update_values(FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FETensor);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::vec fv_vec = fe_data.ref_shape_values[i][j];
                for (unsigned int c=0; c<spacedim*spacedim; c++)
                    fe_values.shape_values[i][j*spacedim*spacedim+c] = fv_vec[c];
            }
    }

    /// Update shape_gradients of given FEValues object.
    inline void update_gradients(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FETensor);

        for (unsigned int i = 0; i < fe_data.n_points; i++)
            for (unsigned int j = 0; j < fe_data.n_dofs; j++)
            {
                arma::mat grads = trans(elm_values.inverse_jacobian(i)) * fe_data.ref_shape_grads[i][j];
                for (unsigned int c=0; c<spacedim*spacedim; c++)
                    fe_values.shape_gradients[i][j*spacedim*spacedim+c] = grads.col(c);
            }
    }
};


/**
 * @brief Helper class allows update values and gradients of FEValues of FEMixedSystem type
 */
template<unsigned int spacedim = 3>
class MapSystem {
public:
    /// Fill fe_values_vec of components of mixed system FEValues object.
    inline void fill_values_vec(FEValues<spacedim> &fe_values, const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEMixedSystem);

        unsigned int comp_offset = 0;
        for (unsigned int f=0; f<fe_values.fe_sys_dofs_.size(); f++)
        {
            // fill fe_values for base FE
            typename FEValues<spacedim>::FEInternalData vec_fe_data(fe_data, fe_values.fe_sys_dofs_[f], comp_offset, fe_values.fe_sys_n_components_[f]);
            fe_values.fe_values_vec[f].fill_data(elm_values, vec_fe_data);

            comp_offset += fe_values.fe_sys_n_components_[f];
        }
    }

    /// Update shape_values of given FEValues object.
    inline void update_values(FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEMixedSystem);

        arma::vec fv_vec;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_values.fe_sys_dofs_.size(); f++)
        {
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_values.fe_sys_dofs_[f].size(); n++)
                    for (unsigned int c=0; c<fe_values.fe_sys_n_space_components_[f]; c++)
                        fe_values.shape_values[i][shape_offset+fe_values.n_components_*n+comp_offset+c] =
                                fe_values.fe_values_vec[f].shape_values[i][n*fe_values.fe_sys_n_space_components_[f]+c];

            comp_offset += fe_values.fe_sys_n_space_components_[f];
            shape_offset += fe_values.fe_sys_dofs_[f].size()*fe_values.n_components_;
        }
    }

    /// Update shape_gradients of given FEValues object.
    inline void update_gradients(FEValues<spacedim> &fe_values, FMT_UNUSED const ElementValues<spacedim> &elm_values,
            const typename FEValues<spacedim>::FEInternalData &fe_data) {
        ASSERT_DBG(fe_values.fe_type_ == FEMixedSystem);

        arma::mat grads;
        unsigned int comp_offset = 0;
        unsigned int shape_offset = 0;
        for (unsigned int f=0; f<fe_values.fe_sys_dofs_.size(); f++)
        {
            // gather fe_values in vectors for FESystem
            for (unsigned int i=0; i<fe_data.n_points; i++)
                for (unsigned int n=0; n<fe_values.fe_sys_dofs_[f].size(); n++)
                    for (unsigned int c=0; c<fe_values.fe_sys_n_space_components_[f]; c++)
                        fe_values.shape_gradients[i][shape_offset+fe_values.n_components_*n+comp_offset+c] =
                                fe_values.fe_values_vec[f].shape_gradients[i][n*fe_values.fe_sys_n_space_components_[f]+c];

            comp_offset += fe_values.fe_sys_n_space_components_[f];
            shape_offset += fe_values.fe_sys_dofs_[f].size()*fe_values.n_components_;
        }
    }
};


#endif // FE_VALUES_MAP_HH_
