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
 * @file    assembly_flow_output.hh
 * @brief
 */

#ifndef ASSEMBLY_FLOW_OUTPUT_HH_
#define ASSEMBLY_FLOW_OUTPUT_HH_

#include <fstream>

#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"
#include "flow/darcy_flow_mh_output.hh"
#include "flow/darcy_flow_lmh.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_fe.hh"

#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "fem/fe_p.hh"
#include "fem/fe_values.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_values_views.hh"
#include "fem/fe_system.hh"
#include "quadrature/quadrature_lib.hh"


/**
 * Calculate approximation of L2 norm for:
 * 1) difference between regularized pressure and analytical solution (using FunctionPython)
 * 2) difference between RT velocities and analytical solution
 * 3) difference of divergence
 *
 * TODO:
 * 1) implement field objects
 * 2) implement DG_P2 finite elements
 * 3) implement pressure postprocessing (result is DG_P2 field)
 * 4) implement calculation of L2 norm for two field (compute the norm and values on individual elements as P0 field)
 *
 */
template <unsigned int dim>
class L2DifferenceAssembly : public AssemblyBase<dim>
{
public:
    typedef typename DarcyLMH::EqFields EqFields;
    typedef typename DarcyFlowMHOutput::DiffEqData EqData;

    static constexpr const char * name() { return "L2DifferenceAssembly"; }

    /// Constructor.
    L2DifferenceAssembly(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(2), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->conductivity;
        this->used_fields_ += eq_fields_->cross_section;
        this->used_fields_ += eq_fields_->ref_pressure;
        this->used_fields_ += eq_fields_->ref_velocity;
        this->used_fields_ += eq_fields_->ref_divergence;
    }

    /// Destructor.
    virtual ~L2DifferenceAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;

        UpdateFlags flags = update_values | update_JxW_values | update_quadrature_points;
        fe_p0_ = std::make_shared< FE_P_disc<dim> >(0);
        fe_values_.initialize(*this->quad_, *fe_p0_, flags);
        fv_rt_.initialize(*this->quad_, fe_rt_, flags);

        fluxes_.resize(dim+1);
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> ele = cell.elm();
        fv_rt_.reinit(ele);
        fe_values_.reinit(ele);

        // get coefficients on the current element
        for (unsigned int li = 0; li < ele->n_sides(); li++) {
            fluxes_[li] = eq_data_->flow_data_->full_solution.get( cell.get_loc_dof_indices()[li] );
            //pressure_traces[li] = diff_eq_data_->dh->side_scalar( *(ele->side( li ) ) );
        }
        // TODO: replace with DHCell getter when available for FESystem component
        double pressure_mean = eq_data_->flow_data_->full_solution.get( cell.get_loc_dof_indices()[ cell.n_dofs()/2 ] );

        double velocity_diff=0, divergence_diff=0, pressure_diff=0, diff;

        // 1d:  mean_x_squared = 1/6 (v0^2 + v1^2 + v0*v1)
        // 2d:  mean_x_squared = 1/12 (v0^2 + v1^2 +v2^2 + v0*v1 + v0*v2 + v1*v2)
        double mean_x_squared=0;
        for(unsigned int i_node=0; i_node < ele->n_nodes(); i_node++ )
            for(unsigned int j_node=0; j_node < ele->n_nodes(); j_node++ )
            {
                mean_x_squared += (i_node == j_node ? 2.0 : 1.0) / ( 6 * dim )   // multiply by 2 on diagonal
                        * arma::dot( *ele.node(i_node), *ele.node(j_node));
            }

        unsigned int i_point=0, oposite_node;
        for (auto p : this->bulk_points(element_patch_idx) )
        {
            q_point_ = fe_values_.point(i_point);

            conductivity_ = eq_fields_->conductivity(p);
            cross_ = eq_fields_->cross_section(p);
            ref_pressure_ = eq_fields_->ref_pressure(p);
            ref_flux_ = eq_fields_->ref_velocity(p);
            ref_divergence_ = eq_fields_->ref_divergence(p);

            // compute postprocesed pressure
            diff = 0;
            for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
                oposite_node = RefElement<dim>::oposite_node(i_shape);

                diff += fluxes_[ i_shape ] *
                                   (  arma::dot( q_point_, q_point_ )/ 2
                                    - mean_x_squared / 2
                                    - arma::dot( q_point_, *ele.node(oposite_node) )
                                    + arma::dot( ele.centre(), *ele.node(oposite_node) )
                                   );
            }

            diff = - (1.0 / conductivity_) * diff / dim / ele.measure() / cross_ + pressure_mean ;
            diff = ( diff - ref_pressure_);
            pressure_diff += diff * diff * fe_values_.JxW(i_point);

            // velocity difference
            flux_in_q_point_.zeros();
            for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) {
                flux_in_q_point_ += fluxes_[ i_shape ]
                                  * fv_rt_.vector_view(0).value(i_shape, i_point)
                                  / cross_;
            }

            flux_in_q_point_ -= ref_flux_;
            velocity_diff += dot(flux_in_q_point_, flux_in_q_point_) * fe_values_.JxW(i_point);

            // divergence diff
            diff = 0;
            for(unsigned int i_shape=0; i_shape < ele->n_sides(); i_shape++) diff += fluxes_[ i_shape ];
            diff = ( diff / ele.measure() / cross_ - ref_divergence_);
            divergence_diff += diff * diff * fe_values_.JxW(i_point);

            i_point++;
        }

        // DHCell constructed with diff fields DH, get DOF indices of actual element
        DHCellAccessor sub_dh_cell = cell.cell_with_other_dh(eq_data_->dh_.get());
        IntIdx idx = sub_dh_cell.get_loc_dof_indices()[0];

        auto velocity_data = eq_data_->vel_diff_ptr->vec();
        velocity_data.set( idx, sqrt(velocity_diff) );
        eq_data_->velocity_error[dim-1] += velocity_diff;
        if (dim == 2 && eq_data_->velocity_mask.size() != 0 ) {
            eq_data_->mask_vel_error += (eq_data_->velocity_mask[ ele.idx() ])? 0 : velocity_diff;
        }

        auto pressure_data = eq_data_->pressure_diff_ptr->vec();
        pressure_data.set( idx, sqrt(pressure_diff) );
        eq_data_->pressure_error[dim-1] += pressure_diff;

        auto div_data = eq_data_->div_diff_ptr->vec();
        div_data.set( idx, sqrt(divergence_diff) );
        eq_data_->div_error[dim-1] += divergence_diff;
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
    	DebugOut() << "l2 norm output\n";

    	eq_data_->mask_vel_error=0;
        for(unsigned int j=0; j<3; j++){
        	eq_data_->pressure_error[j] = 0;
        	eq_data_->velocity_error[j] = 0;
        	eq_data_->div_error[j] = 0;
        }
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        ofstream os( FilePath("solution_error", FilePath::output_file) );

        // square root for L2 norm
        for(unsigned int j=0; j<3; j++){
            eq_data_->pressure_error[j] = sqrt(eq_data_->pressure_error[j]);
            eq_data_->velocity_error[j] = sqrt(eq_data_->velocity_error[j]);
            eq_data_->div_error[j] = sqrt(eq_data_->div_error[j]);
        }
        eq_data_->mask_vel_error = sqrt(eq_data_->mask_vel_error);

        os 	<< "l2 norm output\n\n"
        	<< "pressure error 1d: " << eq_data_->pressure_error[0] << endl
        	<< "pressure error 2d: " << eq_data_->pressure_error[1] << endl
        	<< "pressure error 3d: " << eq_data_->pressure_error[2] << endl
        	<< "velocity error 1d: " << eq_data_->velocity_error[0] << endl
        	<< "velocity error 2d: " << eq_data_->velocity_error[1] << endl
        	<< "velocity error 3d: " << eq_data_->velocity_error[2] << endl
        	<< "masked velocity error 2d: " << eq_data_->mask_vel_error <<endl
        	<< "div error 1d: " << eq_data_->div_error[0] << endl
        	<< "div error 2d: " << eq_data_->div_error[1] << endl
            << "div error 3d: " << eq_data_->div_error[2];
    }


protected:
    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    /// Data objects shared with Flow equation
    EqFields *eq_fields_;
    EqData *eq_data_;

    /// following is used for calculation of postprocessed pressure difference
    /// and comparison to analytical solution
    std::shared_ptr<FE_P_disc<dim>> fe_p0_;
    FEValues<3> fe_values_;

    /// FEValues for velocity.
    FE_RT0<dim> fe_rt_;
    FEValues<3> fv_rt_;

    std::vector<double> fluxes_;                                   ///< Precomputed fluxes on element sides.
    arma::vec3 flux_in_q_point_;                                   ///< Precomputed flux in quadrature point.
    arma::vec3 q_point_;                                           ///< Local coords of quadrature point.
    arma::vec3 ref_flux_;                                          ///< Field result.
    double ref_pressure_, ref_divergence_, conductivity_, cross_;  ///< Field results.

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


/**
 * Compute output of internal flow data.
 */
template <unsigned int dim>
class OutputInternalFlowAssembly : public AssemblyBase<dim>
{
public:
    typedef typename DarcyLMH::EqFields EqFields;
    typedef typename DarcyFlowMHOutput::RawOutputEqData EqData;

    static constexpr const char * name() { return "OutputInternalFlowAssembly"; }

    /// Constructor.
    OutputInternalFlowAssembly(EqFields *eq_fields, EqData *eq_data)
    : AssemblyBase<dim>(0), eq_fields_(eq_fields), eq_data_(eq_data) {
        this->active_integrals_ = ActiveIntegrals::bulk;
        this->used_fields_ += eq_fields_->field_ele_velocity;
    }

    /// Destructor.
    virtual ~OutputInternalFlowAssembly() {}

    /// Initialize auxiliary vectors and other data members
    void initialize(ElementCacheMap *element_cache_map) {
        this->element_cache_map_ = element_cache_map;
    }


    /// Assemble integral over element
    inline void cell_integral(DHCellAccessor cell, unsigned int element_patch_idx)
    {
        ASSERT_EQ(cell.dim(), dim).error("Dimension of element mismatch!");

        ElementAccessor<3> ele = cell.elm();
        LocDofVec indices = cell.get_loc_dof_indices();

        std::stringstream ss;
        // pressure
        ss << fmt::format("{} {} ", cell.elm().input_id(), eq_data_->flow_data_->full_solution.get(indices[ele->n_sides()]));

        // velocity at element center
        auto p = *( this->bulk_points(element_patch_idx).begin() );
        flux_in_center_ = eq_fields_->field_ele_velocity(p);
        for (unsigned int i = 0; i < 3; i++)
        	ss << flux_in_center_[i] << " ";

        // number of sides
        ss << ele->n_sides() << " ";

        // use node permutation to permute sides
        auto &new_to_old_node = ele.orig_nodes_order();
        std::vector<uint> old_to_new_side(ele->n_sides());
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            // According to RefElement<dim>::opposite_node()
            uint new_opp_node = ele->n_sides() - i - 1;
            uint old_opp_node = new_to_old_node[new_opp_node];
            uint old_iside = ele->n_sides() - old_opp_node - 1;
            old_to_new_side[old_iside] = i;
        }

        // pressure on edges
        // unsigned int lid = ele->n_sides() + 1;
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            uint new_lid = ele->n_sides() + 1 + old_to_new_side[i];
            ss << eq_data_->flow_data_->full_solution.get(indices[new_lid]) << " ";
        }
        // fluxes on sides
        for (unsigned int i = 0; i < ele->n_sides(); i++) {
            uint new_iside = old_to_new_side[i];
            ss << eq_data_->flow_data_->full_solution.get(indices[new_iside]) << " ";
        }

        // remove last white space
        string line = ss.str();
        eq_data_->raw_output_strings_[cell.elm_idx()] = line.substr(0, line.size()-1);
    }

    /// Implements @p AssemblyBase::begin.
    void begin() override
    {
    	DebugOut() << "Internal flow data output\n";

    	eq_data_->raw_output_strings_.resize( eq_data_->flow_data_->dh_->n_own_cells() );
    }

    /// Implements @p AssemblyBase::end.
    void end() override
    {
        //char dbl_fmt[ 16 ]= "%.8g ";
        // header
        eq_data_->raw_output_file <<  "// fields:\n//ele_id    ele_presure    flux_in_barycenter[3]    n_sides   side_pressures[n]    side_fluxes[n]\n";
        eq_data_->raw_output_file <<  fmt::format("$FlowField\nT={}\n", eq_data_->time_->t());
        eq_data_->raw_output_file <<  fmt::format("{}\n" , eq_data_->flow_data_->dh_->mesh()->n_elements() );

        auto permutation_vec = eq_data_->flow_data_->dh_->mesh()->element_permutations();
        for (unsigned int i_elem=0; i_elem<eq_data_->flow_data_->dh_->n_own_cells(); ++i_elem) {
            eq_data_->raw_output_file << eq_data_->raw_output_strings_[ permutation_vec[i_elem] ] << endl;
        }

        eq_data_->raw_output_file << "$EndFlowField\n" << endl;
    }


protected:
    /// Sub field set contains fields used in calculation.
    FieldSet used_fields_;

    /// Data objects shared with Flow equation
    EqFields *eq_fields_;
    EqData *eq_data_;

    arma::vec3 flux_in_center_;

    template < template<IntDim...> class DimAssembly>
    friend class GenericAssembly;

};


#endif /* ASSEMBLY_FLOW_OUTPUT_HH_ */
