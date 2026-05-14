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
 * @file    transport_dg.cc
 * @brief   Discontinuous Galerkin method for equation of transport with dispersion.
 * @author  Jan Stebel
 */

#include <boost/lexical_cast.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "system/sys_profiler.hh"
#include "mechanics/elasticity.hh"
#include "mechanics/assembly_elasticity.hh"

#include "io/output_time.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_values.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/fe_system.hh"
#include "fields/field_fe.hh"
#include "la/linsys_PETSC.hh"
#include "la/linsys_PERMON.hh"
#include "coupling/balance.hh"
#include "mesh/neighbours.h"
#include "coupling/generic_assembly.hh"

#include "fields/multi_field.hh"
#include "fields/generic_field.hh"
#include "fields/field_model.hh"
#include "input/factory.hh"




using namespace Input::Type;

namespace {

struct ResidualSplit {
    double norm = 0.0;
    double free_norm = 0.0;
    double support_norm = 0.0;
    double filtered_norm = 0.0;
    double filtered_free_norm = 0.0;
    double filtered_support_norm = 0.0;
    double free_interface_norm = 0.0;
    double free_interior_norm = 0.0;
    double filtered_free_interface_norm = 0.0;
    double filtered_free_interior_norm = 0.0;
    double solution_interface_jump_norm = 0.0;
    double solution_interface_jump_max = 0.0;
    double free_max = 0.0;
    double support_max = 0.0;
    int free_nnz = 0;
    int support_nnz = 0;
    int free_count = 0;
    int support_count = 0;
    int solution_interface_jump_nnz = 0;
};

std::vector<int> make_interface_mask(Mat matrix, PetscInt global_size, MPI_Comm comm)
{
    PetscBool is_matis = PETSC_FALSE;
    chkerr(PetscObjectTypeCompare((PetscObject)matrix, MATIS, &is_matis));
    if (!is_matis) return std::vector<int>(global_size, 0);

    ISLocalToGlobalMapping l2g_mapping = NULL;
    const PetscInt *l2g_arr = NULL;
    PetscInt nmap = 0;
    std::vector<int> local_count(global_size, 0), global_count(global_size, 0);

    chkerr(MatISGetLocalToGlobalMapping(matrix, &l2g_mapping, NULL));
    chkerr(ISLocalToGlobalMappingGetIndices(l2g_mapping, &l2g_arr));
    chkerr(ISLocalToGlobalMappingGetSize(l2g_mapping, &nmap));
    for (PetscInt i=0; i<nmap; i++) {
        ASSERT_LT(static_cast<unsigned int>(l2g_arr[i]), static_cast<unsigned int>(global_size))
            .error("MatIS local-to-global index out of bounds of residual vector.\n");
        local_count[l2g_arr[i]]++;
    }
    chkerr(ISLocalToGlobalMappingRestoreIndices(l2g_mapping, &l2g_arr));
    chkerr(MPI_Allreduce(local_count.data(), global_count.data(), global_size, MPI_INT, MPI_SUM, comm));

    std::vector<int> interface_mask(global_size, 0);
    for (PetscInt i=0; i<global_size; i++) {
        interface_mask[i] = (global_count[i] > 1);
    }
    return interface_mask;
}

ResidualSplit compute_residual_split_by_dirichlet_dofs(Elasticity::EqData *eq_data)
{
    Mat matrix = *eq_data->ls->get_matrix();
    Vec rhs = *eq_data->ls->get_rhs();
    const Vec &solution = eq_data->ls->get_solution();
    MPI_Comm comm = PetscObjectComm((PetscObject)rhs);
    Vec residual = NULL;
    PetscInt global_size = 0, own_begin = 0, own_end = 0;
    const PetscScalar *residual_array = NULL;
    ResidualSplit split;

    chkerr(VecGetSize(rhs, &global_size));
    std::vector<int> interface_mask = make_interface_mask(matrix, global_size, comm);
    std::vector<int> local_support(global_size, 0), global_support(global_size, 0);
    const auto &l2g = eq_data->dh_->get_local_to_global_map();
    for (unsigned int local_dof : eq_data->dirichlet_dofs) {
        ASSERT_LT(local_dof, static_cast<unsigned int>(l2g.size()))
            .error("Dirichlet DOF index out of bounds of local-to-global map.\n");
        ASSERT_LT(static_cast<unsigned int>(l2g[local_dof]), static_cast<unsigned int>(global_size))
            .error("Dirichlet DOF index out of bounds of residual vector.\n");
        local_support[l2g[local_dof]] = 1;
    }
    chkerr(MPI_Allreduce(local_support.data(), global_support.data(), global_size, MPI_INT, MPI_MAX, comm));

    PetscBool is_matis = PETSC_FALSE;
    chkerr(PetscObjectTypeCompare((PetscObject)matrix, MATIS, &is_matis));
    if (is_matis) {
        ISLocalToGlobalMapping l2g_mapping = NULL;
        const PetscInt *matis_l2g_arr = NULL;
        Vec solution_local = NULL;
        IS from = NULL, to = NULL;
        VecScatter solution_scatter = NULL;
        const PetscScalar *solution_local_arr = NULL;
        PetscInt nmap = 0;
        std::vector<int> local_solution_count(global_size, 0), global_solution_count(global_size, 0);
        std::vector<double> local_solution_min(global_size, std::numeric_limits<double>::infinity());
        std::vector<double> local_solution_max(global_size, -std::numeric_limits<double>::infinity());
        std::vector<double> global_solution_min(global_size, 0.0), global_solution_max(global_size, 0.0);

        chkerr(MatISGetLocalToGlobalMapping(matrix, &l2g_mapping, NULL));
        chkerr(ISLocalToGlobalMappingGetIndices(l2g_mapping, &matis_l2g_arr));
        chkerr(ISLocalToGlobalMappingGetSize(l2g_mapping, &nmap));
        chkerr(VecCreateSeq(PETSC_COMM_SELF, nmap, &solution_local));
        chkerr(ISCreateGeneral(comm, nmap, matis_l2g_arr, PETSC_COPY_VALUES, &from));
        chkerr(ISCreateStride(PETSC_COMM_SELF, nmap, 0, 1, &to));
        chkerr(VecScatterCreate(solution, from, solution_local, to, &solution_scatter));
        chkerr(VecScatterBegin(solution_scatter, solution, solution_local, INSERT_VALUES, SCATTER_FORWARD));
        chkerr(VecScatterEnd(solution_scatter, solution, solution_local, INSERT_VALUES, SCATTER_FORWARD));
        chkerr(VecGetArrayRead(solution_local, &solution_local_arr));

        for (PetscInt i=0; i<nmap; i++) {
            PetscInt global_idx = matis_l2g_arr[i];
            ASSERT_LT(static_cast<unsigned int>(global_idx), static_cast<unsigned int>(global_size))
                .error("MatIS local-to-global index out of bounds of solution vector.\n");
            double value = static_cast<double>(PetscRealPart(solution_local_arr[i]));
            local_solution_count[global_idx]++;
            local_solution_min[global_idx] = std::min(local_solution_min[global_idx], value);
            local_solution_max[global_idx] = std::max(local_solution_max[global_idx], value);
        }

        chkerr(MPI_Allreduce(local_solution_count.data(), global_solution_count.data(), global_size, MPI_INT, MPI_SUM, comm));
        chkerr(MPI_Allreduce(local_solution_min.data(), global_solution_min.data(), global_size, MPI_DOUBLE, MPI_MIN, comm));
        chkerr(MPI_Allreduce(local_solution_max.data(), global_solution_max.data(), global_size, MPI_DOUBLE, MPI_MAX, comm));

        double jump_norm_sq = 0.0;
        for (PetscInt i=0; i<global_size; i++) {
            if (global_solution_count[i] > 1) {
                double jump = global_solution_max[i] - global_solution_min[i];
                jump_norm_sq += jump * jump;
                split.solution_interface_jump_max = std::max(split.solution_interface_jump_max, jump);
                if (jump > 1e-12) split.solution_interface_jump_nnz++;
            }
        }
        split.solution_interface_jump_norm = std::sqrt(jump_norm_sq);

        chkerr(VecRestoreArrayRead(solution_local, &solution_local_arr));
        chkerr(VecScatterDestroy(&solution_scatter));
        chkerr(ISDestroy(&to));
        chkerr(ISDestroy(&from));
        chkerr(VecDestroy(&solution_local));
        chkerr(ISLocalToGlobalMappingRestoreIndices(l2g_mapping, &matis_l2g_arr));
    }

    chkerr(VecDuplicate(rhs, &residual));
    chkerr(MatMult(matrix, solution, residual));
    chkerr(VecAXPY(residual, -1.0, rhs));
    chkerr(VecGetOwnershipRange(residual, &own_begin, &own_end));
    chkerr(VecGetArrayRead(residual, &residual_array));

    double local_free_norm_sq = 0.0, local_support_norm_sq = 0.0;
    double local_free_interface_norm_sq = 0.0, local_free_interior_norm_sq = 0.0;
    double local_free_max = 0.0, local_support_max = 0.0;
    int local_free_nnz = 0, local_support_nnz = 0;
    int local_free_count = 0, local_support_count = 0;
    for (PetscInt i=0; i<own_end-own_begin; i++) {
        PetscInt global_idx = own_begin + i;
        double val_abs = PetscAbsScalar(residual_array[i]);
        if (global_support[global_idx]) {
            local_support_norm_sq += val_abs * val_abs;
            local_support_max = std::max(local_support_max, val_abs);
            local_support_count++;
            if (val_abs > 1e-12) local_support_nnz++;
        } else {
            local_free_norm_sq += val_abs * val_abs;
            if (interface_mask[global_idx]) {
                local_free_interface_norm_sq += val_abs * val_abs;
            } else {
                local_free_interior_norm_sq += val_abs * val_abs;
            }
            local_free_max = std::max(local_free_max, val_abs);
            local_free_count++;
            if (val_abs > 1e-12) local_free_nnz++;
        }
    }

    chkerr(VecRestoreArrayRead(residual, &residual_array));
    chkerr(VecDestroy(&residual));

    double free_norm_sq = 0.0, support_norm_sq = 0.0;
    double free_interface_norm_sq = 0.0, free_interior_norm_sq = 0.0;
    chkerr(MPI_Allreduce(&local_free_norm_sq, &free_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_support_norm_sq, &support_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_free_interface_norm_sq, &free_interface_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_free_interior_norm_sq, &free_interior_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_free_max, &split.free_max, 1, MPI_DOUBLE, MPI_MAX, comm));
    chkerr(MPI_Allreduce(&local_support_max, &split.support_max, 1, MPI_DOUBLE, MPI_MAX, comm));
    chkerr(MPI_Allreduce(&local_free_nnz, &split.free_nnz, 1, MPI_INT, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_support_nnz, &split.support_nnz, 1, MPI_INT, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_free_count, &split.free_count, 1, MPI_INT, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_support_count, &split.support_count, 1, MPI_INT, MPI_SUM, comm));

    split.free_norm = std::sqrt(free_norm_sq);
    split.support_norm = std::sqrt(support_norm_sq);
    split.norm = std::sqrt(free_norm_sq + support_norm_sq);
    split.free_interface_norm = std::sqrt(free_interface_norm_sq);
    split.free_interior_norm = std::sqrt(free_interior_norm_sq);

    if (!is_matis) {
        split.filtered_norm = split.norm;
        split.filtered_free_norm = split.free_norm;
        split.filtered_support_norm = split.support_norm;
        split.filtered_free_interface_norm = split.free_interface_norm;
        split.filtered_free_interior_norm = split.free_interior_norm;
        return split;
    }

    Mat Aloc = NULL, Aloc_filtered = NULL;
    Mat matrix_filtered = NULL;
    IS isnz = NULL, isnz_tmp = NULL;
    ISLocalToGlobalMapping l2g_mapping = NULL;
    ISLocalToGlobalMapping filtered_mapping = NULL;
    const PetscInt *l2g_arr = NULL, *nz_arr = NULL;
    PetscInt *l2g_filtered_arr = NULL;
    Vec residual_filtered = NULL;
    PetscInt nmap = 0, nrows = 0;
    double local_filtered_free_norm_sq = 0.0, local_filtered_support_norm_sq = 0.0;
    double local_filtered_free_interface_norm_sq = 0.0, local_filtered_free_interior_norm_sq = 0.0;
    double filtered_free_norm_sq = 0.0, filtered_support_norm_sq = 0.0;
    double filtered_free_interface_norm_sq = 0.0, filtered_free_interior_norm_sq = 0.0;

    chkerr(MatISGetLocalMat(matrix, &Aloc));
    chkerr(MatFindNonzeroRows(Aloc, &isnz_tmp));
    if (isnz_tmp != NULL) {
        const PetscInt *nz_idx = NULL;
        chkerr(ISGetLocalSize(isnz_tmp, &nrows));
        chkerr(ISGetIndices(isnz_tmp, &nz_idx));
        chkerr(ISCreateGeneral(PETSC_COMM_SELF, nrows, nz_idx, PETSC_COPY_VALUES, &isnz));
        chkerr(ISRestoreIndices(isnz_tmp, &nz_idx));
    } else {
        chkerr(MatGetLocalSize(Aloc, &nrows, NULL));
        chkerr(ISCreateStride(PETSC_COMM_SELF, nrows, 0, 1, &isnz));
    }
    chkerr(ISDestroy(&isnz_tmp));
    chkerr(MatCreateSubMatrix(Aloc, isnz, isnz, MAT_INITIAL_MATRIX, &Aloc_filtered));
    chkerr(ISGetIndices(isnz, &nz_arr));
    chkerr(MatISGetLocalToGlobalMapping(matrix, &l2g_mapping, NULL));
    chkerr(ISLocalToGlobalMappingGetIndices(l2g_mapping, &l2g_arr));
    chkerr(ISLocalToGlobalMappingGetSize(l2g_mapping, &nmap));
    chkerr(PetscMalloc1(nrows, &l2g_filtered_arr));
    for (PetscInt i=0; i<nrows; i++) {
        ASSERT_LT(static_cast<unsigned int>(nz_arr[i]), static_cast<unsigned int>(nmap))
            .error("Filtered residual local row index out of bounds of local-to-global mapping.\n");
        l2g_filtered_arr[i] = l2g_arr[nz_arr[i]];
    }
    chkerr(ISRestoreIndices(isnz, &nz_arr));
    chkerr(ISLocalToGlobalMappingRestoreIndices(l2g_mapping, &l2g_arr));
    chkerr(ISLocalToGlobalMappingCreate(comm, 1, nrows, l2g_filtered_arr, PETSC_OWN_POINTER, &filtered_mapping));
    chkerr(MatCreateIS(comm, 1, eq_data->dh_->lsize(), eq_data->dh_->lsize(),
                       PETSC_DETERMINE, PETSC_DETERMINE, filtered_mapping, filtered_mapping, &matrix_filtered));
    chkerr(MatISSetLocalMat(matrix_filtered, Aloc_filtered));
    chkerr(MatAssemblyBegin(matrix_filtered, MAT_FINAL_ASSEMBLY));
    chkerr(MatAssemblyEnd(matrix_filtered, MAT_FINAL_ASSEMBLY));

    chkerr(VecDuplicate(rhs, &residual_filtered));
    chkerr(MatMult(matrix_filtered, solution, residual_filtered));
    chkerr(VecAXPY(residual_filtered, -1.0, rhs));
    chkerr(VecGetOwnershipRange(residual_filtered, &own_begin, &own_end));
    chkerr(VecGetArrayRead(residual_filtered, &residual_array));
    for (PetscInt i=0; i<own_end-own_begin; i++) {
        PetscInt global_idx = own_begin + i;
        double val_abs = PetscAbsScalar(residual_array[i]);
        if (global_support[global_idx]) {
            local_filtered_support_norm_sq += val_abs * val_abs;
        } else {
            local_filtered_free_norm_sq += val_abs * val_abs;
            if (interface_mask[global_idx]) {
                local_filtered_free_interface_norm_sq += val_abs * val_abs;
            } else {
                local_filtered_free_interior_norm_sq += val_abs * val_abs;
            }
        }
    }
    chkerr(VecRestoreArrayRead(residual_filtered, &residual_array));

    chkerr(MPI_Allreduce(&local_filtered_free_norm_sq, &filtered_free_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_filtered_support_norm_sq, &filtered_support_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_filtered_free_interface_norm_sq, &filtered_free_interface_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    chkerr(MPI_Allreduce(&local_filtered_free_interior_norm_sq, &filtered_free_interior_norm_sq, 1, MPI_DOUBLE, MPI_SUM, comm));
    split.filtered_free_norm = std::sqrt(filtered_free_norm_sq);
    split.filtered_support_norm = std::sqrt(filtered_support_norm_sq);
    split.filtered_norm = std::sqrt(filtered_free_norm_sq + filtered_support_norm_sq);
    split.filtered_free_interface_norm = std::sqrt(filtered_free_interface_norm_sq);
    split.filtered_free_interior_norm = std::sqrt(filtered_free_interior_norm_sq);

    chkerr(VecDestroy(&residual_filtered));
    chkerr(MatDestroy(&matrix_filtered));
    chkerr(ISLocalToGlobalMappingDestroy(&filtered_mapping));
    chkerr(MatDestroy(&Aloc_filtered));
    chkerr(ISDestroy(&isnz));
    chkerr(MatISRestoreLocalMat(matrix, &Aloc));
    return split;
}

} // namespace



const Record & Elasticity::get_input_type() {
    std::string equation_name = std::string(name_) + "_FE";
	return IT::Record(
                std::string(equation_name),
                "FEM for linear elasticity.")
           .copy_keys(EquationBase::record_template())
		   .copy_keys(EquationBase::user_fields_template(equation_name))
           .declare_key("balance", Balance::get_input_type(), Default("{}"),
                    "Settings for computing balance.")
           .declare_key("output_stream", OutputTime::get_input_type(), Default::obligatory(),
                    "Parameters of output stream.")
           .declare_key("solver", LinSys_PERMON::get_input_type(), Default::obligatory(),
				"Linear solver for elasticity.")
		   .declare_key("input_fields", Array(
		        Elasticity::EqFields()
		            .make_field_descriptor_type(equation_name)),
		        IT::Default::obligatory(),
		        "Input fields of the equation.")
           .declare_key("output",
                EqFields().output_fields.make_output_type(equation_name, ""),
                IT::Default("{ \"fields\": [ \"displacement\" ] }"),
                "Setting of the field output.")
           .declare_key("contact", Bool(), IT::Default("false"), "Indicates the use of contact conditions on fractures.")
           .declare_key("feti", Bool(), IT::Default("false"), "Use FETI domain decomposition.")
           .declare_key("view_matrices", Bool(), IT::Default("false"), "Print system matrices and vectors.")
           .declare_key("fix_nullspace", Bool(), IT::Default("false"), "Correct formulation to preserve rotational DOF in nullspace of stiffness matrix.")
           .declare_key("dirichlet_by_eq", Bool(), IT::Default("false"), "Enforce Dirichlet b.c. by equality constraints instead of penalization.")
		   .close();
}

const int Elasticity::registrar =
		Input::register_class< Elasticity, Mesh &, const Input::Record>(std::string(name_) + "_FE") +
		Elasticity::get_input_type().size();



double lame_mu(double young, double poisson)
{
    return young*0.5/(poisson+1.);
}


double lame_lambda(double young, double poisson)
{
    return young*poisson/((poisson+1.)*(1.-2.*poisson));
}

// Functor computing lame_mu
struct fn_lame_mu {
	inline double operator() (double young, double poisson) {
        return young * 0.5 / (poisson+1.);
    }
};

// Functor computing lame_lambda
struct fn_lame_lambda {
	inline double operator() (double young, double poisson) {
        return young * poisson / ((poisson+1.)*(1.-2.*poisson));
    }
};

// Functor computing base of dirichlet_penalty (without dividing by side meassure)
struct fn_dirichlet_penalty {
	inline double operator() (double lame_mu, double lame_lambda) {
        return 1e0 * (2 * lame_mu + lame_lambda);
    }
};






const Selection & Elasticity::EqFields::get_bc_type_selection() {
    return Selection("Elasticity_BC_Type", "Types of boundary conditions for mechanics.")
            .add_value(bc_type_displacement, "displacement",
                  "Prescribed displacement.")
            .add_value(bc_type_displacement_normal, "displacement_n",
                  "Prescribed displacement in the normal direction to the boundary.")
            .add_value(bc_type_traction, "traction",
                  "Prescribed traction.")
            .add_value(bc_type_stress, "stress",
                  "Prescribed stress tensor.")
            .close();
}


Elasticity::EqFields::EqFields()
{
    *this+=bc_type
        .name("bc_type")
        .description(
        "Type of boundary condition.")
        .units( UnitSI::dimensionless() )
        .input_default("\"traction\"")
        .input_selection( get_bc_type_selection() )
        .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
        
    *this+=bc_displacement
        .name("bc_displacement")
        .description("Prescribed displacement on boundary.")
        .units( UnitSI().m() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=bc_traction
        .name("bc_traction")
        .description("Prescribed traction on boundary.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=bc_stress
        .name("bc_stress")
        .description("Prescribed stress on boundary.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);
        
    *this+=load
        .name("load")
        .description("Prescribed bulk load.")
        .units( UnitSI().N().m(-3) )
        .input_default("0.0")
        .flags_add(in_rhs);
    
    *this+=young_modulus
        .name("young_modulus")
        .description("Young's modulus.")
        .units( UnitSI().Pa() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
    
    *this+=poisson_ratio
        .name("poisson_ratio")
        .description("Poisson's ratio.")
        .units( UnitSI().dimensionless() )
         .input_default("0.0")
        .flags_add(in_main_matrix & in_rhs);
        
    *this+=fracture_sigma
            .name("fracture_sigma")
            .description(
            "Coefficient of transfer of forces through fractures.")
            .units( UnitSI::dimensionless() )
            .input_default("1.0")
            .flags_add(in_main_matrix & in_rhs);

    *this+=initial_stress
        .name("initial_stress")
        .description("Initial stress tensor.")
        .units( UnitSI().Pa() )
        .input_default("0.0")
        .flags_add(in_rhs);

    *this += region_id.name("region_id")
    	        .units( UnitSI::dimensionless())
    	        .flags(FieldFlag::equation_external_output);
                
    *this += subdomain.name("subdomain")
      .units( UnitSI::dimensionless() )
      .flags(FieldFlag::equation_external_output);
      
    *this+=cross_section
      .name("cross_section")
      .units( UnitSI().m(3).md() )
      .flags(input_copy & in_time_term & in_main_matrix & in_rhs);
    
    *this+=cross_section_min
      .name("cross_section_min")
      .description("Minimal cross-section of fractures.")
      .units( UnitSI().m(3).md() )
      .input_default("0.0");
      
    *this+=potential_load
      .name("potential_load")
      .units( UnitSI().m() )
      .flags(input_copy & in_rhs);

    *this+=output_field
            .name("displacement")
            .description("Displacement vector field output.")
            .units( UnitSI().m() )
            .flags(equation_result);
    
    *this += output_stress
            .name("stress")
            .description("Stress tensor output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_von_mises_stress
            .name("von_mises_stress")
            .description("von Mises stress output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);
    
    *this += output_mean_stress
            .name("mean_stress")
            .description("mean stress output.")
            .units( UnitSI().Pa() )
            .flags(equation_result);

    *this += output_cross_section
            .name("cross_section_updated")
            .description("Cross-section after deformation - output.")
            .units( UnitSI().m() )
            .flags(equation_result);
            
    *this += output_divergence
            .name("displacement_divergence")
            .description("Displacement divergence output.")
            .units( UnitSI().dimensionless() )
            .flags(equation_result);

    *this += lame_mu.name("lame_mu")
            .description("Field lame_mu.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    *this += lame_lambda.name("lame_lambda")
            .description("Field lame_lambda.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    *this += dirichlet_penalty.name("dirichlet_penalty")
            .description("Field dirichlet_penalty.")
            .input_default("0.0")
            .units( UnitSI().Pa() );

    // add all input fields to the output list
    output_fields += *this;

    this->add_coords_field();
    this->set_default_fieldset();

}

Elasticity::EqData::~EqData()
{
    if (ls!=nullptr) delete ls;
    if (constraint_matrix!=nullptr) MatDestroy(&constraint_matrix);
    if (constraint_vec!=nullptr) VecDestroy(&constraint_vec);
    if (dirichlet_matrix!=nullptr) MatDestroy(&dirichlet_matrix);
    if (dirichlet_vec!=nullptr) VecDestroy(&dirichlet_vec);
}

void Elasticity::EqData::create_dh(Mesh * mesh, unsigned int fe_order)
{
	ASSERT_EQ(fe_order, 1)(fe_order).error("Unsupported polynomial order for finite elements in Elasticity");
    MixedPtr<FE_P> fe_p(1);
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);

    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
	dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh);

	dh_->distribute_dofs(ds);


    MixedPtr<FE_P_disc> fe_p_disc(0);
    dh_scalar_ = make_shared<DOFHandlerMultiDim>(*mesh);
	std::shared_ptr<DiscreteSpace> ds_scalar = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_p_disc);
	dh_scalar_->distribute_dofs(ds_scalar);


    MixedPtr<FiniteElement> fe_t = mixed_fe_system(MixedPtr<FE_P_disc>(0), FEType::FETensor, 9);
    dh_tensor_ = make_shared<DOFHandlerMultiDim>(*mesh);
	std::shared_ptr<DiscreteSpace> dst = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe_t);
	dh_tensor_->distribute_dofs(dst);
}


Elasticity::Elasticity(Mesh & init_mesh, const Input::Record in_rec, TimeGovernor *tm)
        : EquationBase(init_mesh, in_rec),
		  input_rec(in_rec),
		  stiffness_assembly_(nullptr),
		  rhs_assembly_(nullptr),
          constraint_assembly_(nullptr),
		  output_fields_assembly_(nullptr),
          dirichlet_assembly_(nullptr)
{
	// Can not use name() + "constructor" here, since START_TIMER only accepts const char *
	// due to constexpr optimization.
	START_TIMER("Mechanics constructor");

    eq_data_ = std::make_shared<EqData>();
    eq_fields_ = std::make_shared<EqFields>();
    this->eq_fieldset_ = eq_fields_;
    
    auto time_rec = in_rec.val<Input::Record>("time");
    if (tm == nullptr)
    {
        time_ = new TimeGovernor(time_rec);
    }
    else
    {
        TimeGovernor time_from_rec(time_rec);
        ASSERT( time_from_rec.is_default() ).error("Duplicate key 'time', time in elasticity is already initialized from parent class!");
        time_ = tm;
    }


    // Set up physical parameters.
    eq_fields_->set_mesh(init_mesh);
    eq_fields_->region_id = GenericField<3>::region_id(*mesh_);
    eq_fields_->subdomain = GenericField<3>::subdomain(*mesh_);
    eq_data_->balance_ = this->balance();
    
    // create finite element structures and distribute DOFs
    eq_data_->create_dh(mesh_, 1);
    DebugOut().fmt("Mechanics: solution size {}\n", eq_data_->dh_->n_global_dofs());
    END_TIMER("Mechanics constructor");
}


void Elasticity::initialize()
{
    output_stream_ = OutputTime::create_output_stream("mechanics", input_rec.val<Input::Record>("output_stream"), time().get_unit_conversion());
    
    eq_fields_->set_components({"displacement"});
    eq_fields_->set_input_list( input_rec.val<Input::Array>("input_fields"), time() );
    
//     balance_ = std::make_shared<Balance>("mechanics", mesh_);
//     balance_->init_from_input(input_rec.val<Input::Record>("balance"), *time_);
    // initialization of balance object
//     eq_data_->balance_idx_ = {
//         balance_->add_quantity("force_x"),
//         balance_->add_quantity("force_y"),
//         balance_->add_quantity("force_z")
//     };
//     balance_->units(UnitSI().kg().m().s(-2));

    // create shared pointer to a FieldFE, pass FE data and push this FieldFE to output_field on all regions
    eq_fields_->output_field_ptr = create_field_fe<3, FieldValue<3>::VectorFixed>(eq_data_->dh_);
    eq_fields_->output_field.set(eq_fields_->output_field_ptr, 0.);
    
    // setup output stress
    eq_fields_->output_stress_ptr = create_field_fe<3, FieldValue<3>::TensorFixed>(eq_data_->dh_tensor_);
    eq_fields_->output_stress.set(eq_fields_->output_stress_ptr, 0.);
    
    // setup output von Mises stress
    eq_fields_->output_von_mises_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_von_mises_stress.set(eq_fields_->output_von_mises_stress_ptr, 0.);

    // setup output mean stress
    eq_fields_->output_mean_stress_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_mean_stress.set(eq_fields_->output_mean_stress_ptr, 0.);
    
    // setup output cross-section
    eq_fields_->output_cross_section_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_cross_section.set(eq_fields_->output_cross_section_ptr, 0.);
    
    // setup output divergence
    eq_fields_->output_div_ptr = create_field_fe<3, FieldValue<3>::Scalar>(eq_data_->dh_scalar_);
    eq_fields_->output_divergence.set(eq_fields_->output_div_ptr, 0.);

    // read optional user fields
    Input::Array user_fields_arr;
    if (input_rec.opt_val("user_fields", user_fields_arr)) {
       	this->init_user_fields(user_fields_arr, eq_fields_->output_fields);
    }
    
    eq_fields_->output_fields.set_mesh(*mesh_);
    eq_fields_->output_field.output_type(OutputTime::CORNER_DATA);

    // set time marks for writing the output
    eq_fields_->output_fields.initialize(output_stream_, mesh_, input_rec.val<Input::Record>("output"), this->time());

    // set instances of FieldModel
    eq_fields_->lame_mu.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_mu(), eq_fields_->young_modulus, eq_fields_->poisson_ratio), 0.0);
    eq_fields_->lame_lambda.set(Model<3, FieldValue<3>::Scalar>::create(fn_lame_lambda(), eq_fields_->young_modulus, eq_fields_->poisson_ratio), 0.0);
    eq_fields_->dirichlet_penalty.set(Model<3, FieldValue<3>::Scalar>::create(fn_dirichlet_penalty(), eq_fields_->lame_mu, eq_fields_->lame_lambda), 0.0);

    // equation default PETSc solver options
    std::string petsc_default_opts;
    petsc_default_opts = "-permon_pc_type hypre";
    
    // allocate matrix and vector structures
#ifndef FLOW123D_HAVE_PERMON
    ASSERT(false).error("Flow123d was not built with PERMON library, therefore contact conditions are unsupported.");
#endif //FLOW123D_HAVE_PERMON
    LinSys_PERMON *ls = new LinSys_PERMON(*eq_data_->dh_, petsc_default_opts);
    has_contact_ = input_rec.val<bool>("contact");
    if (has_contact_) {
        // allocate constraint matrix and vector
        eq_data_->constraint_idx.clear();
        eq_data_->constraint_idx_local.clear();
        unsigned int n_own_constraints = 0; // count locally owned cells with neighbours
        for (auto cell : eq_data_->dh_->own_range())
            if (cell.elm()->n_neighs_vb() > 0)
                n_own_constraints++;
        unsigned int n_constraints = 0; // count all cells with neighbours
        for (auto elm : mesh_->elements_range())
            if (elm->n_neighs_vb() > 0)
                eq_data_->constraint_idx[elm.idx()] = n_constraints++;
        std::vector<PetscInt> constraint_rows_l2g;
        constraint_rows_l2g.reserve(n_own_constraints);
        unsigned int local_constraint_idx = 0;
        for (auto cell : eq_data_->dh_->own_range()) {
            if (cell.elm()->n_neighs_vb() > 0) {
                eq_data_->constraint_idx_local[cell.elm_idx()] = local_constraint_idx++;
                constraint_rows_l2g.push_back(eq_data_->constraint_idx[cell.elm_idx()]);
            }
        }
        unsigned int nnz = eq_data_->dh_->ds()->fe()[1_d]->n_dofs()*mesh_->max_edge_sides(1) +
                        eq_data_->dh_->ds()->fe()[2_d]->n_dofs()*mesh_->max_edge_sides(2) +
                        eq_data_->dh_->ds()->fe()[3_d]->n_dofs()*mesh_->max_edge_sides(3);
        ISLocalToGlobalMapping row_l2g = NULL, col_l2g = NULL;
        Mat constraint_matrix_local = NULL;
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, n_own_constraints, constraint_rows_l2g.data(), PETSC_COPY_VALUES, &row_l2g);
        ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, eq_data_->dh_->get_local_to_global_map().size(),
                                     eq_data_->dh_->get_local_to_global_map().data(), PETSC_COPY_VALUES, &col_l2g);
        MatCreateIS(PETSC_COMM_WORLD, 1, n_own_constraints, eq_data_->dh_->lsize(), PETSC_DECIDE, PETSC_DECIDE,
                    row_l2g, col_l2g, &eq_data_->constraint_matrix);
        MatISGetLocalMat(eq_data_->constraint_matrix, &constraint_matrix_local);
        MatSeqAIJSetPreallocation(constraint_matrix_local, nnz, NULL);
        MatISRestoreLocalMat(eq_data_->constraint_matrix, &constraint_matrix_local);
        MatCreateVecs(eq_data_->constraint_matrix, NULL, &eq_data_->constraint_vec);
        ISLocalToGlobalMappingDestroy(&row_l2g);
        ISLocalToGlobalMappingDestroy(&col_l2g);
        ls->set_inequality(eq_data_->constraint_matrix,eq_data_->constraint_vec);

        constraint_assembly_ = new GenericAssembly< ConstraintAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    }
    ls->set_from_input( input_rec.val<Input::Record>("solver") );
    ls->set_solution(eq_fields_->output_field_ptr->vec().petsc_vec());
    bool use_feti = input_rec.val<bool>("feti");
    ls->use_feti(use_feti);
    eq_data_->ls = ls;

    stiffness_assembly_ = new GenericAssembly< StiffnessAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    rhs_assembly_ = new GenericAssembly< RhsAssemblyElasticity >(eq_fields_.get(), eq_data_.get());
    output_fields_assembly_ = new GenericAssembly< OutpuFieldsAssemblyElasticity >(eq_fields_.get(), eq_data_.get());

    eq_data_->dirichlet_by_eq = input_rec.val<bool>("dirichlet_by_eq");
    // Used for equality constraints and residual split over penalty Dirichlet rows.
    dirichlet_assembly_ = new GenericAssembly<DirichletAssemblyElasticity>(eq_fields_.get(), eq_data_.get());

    eq_data_->fix_nullspace = input_rec.val<bool>("fix_nullspace");
    view_matrices = input_rec.val<bool>("view_matrices");

    // initialization of balance object
//     balance_->allocate(eq_data_->dh_->distr()->lsize(),
//             max(feo->fe<1>()->n_dofs(), max(feo->fe<2>()->n_dofs(), feo->fe<3>()->n_dofs())));

}


Elasticity::~Elasticity()
{
//     delete time_;

    if (stiffness_assembly_!=nullptr) delete stiffness_assembly_;
    if (rhs_assembly_!=nullptr) delete rhs_assembly_;
    if (constraint_assembly_ != nullptr) delete constraint_assembly_;
    if (output_fields_assembly_!=nullptr) delete output_fields_assembly_;
    if (dirichlet_assembly_!=nullptr) delete dirichlet_assembly_;

    eq_data_.reset();
    eq_fields_.reset();
}



void Elasticity::update_output_fields()
{
    eq_fields_->set_time(time_->step(), LimitSide::right);
    
    // update ghost values of solution vector and prepare dependent fields
	eq_fields_->output_field_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_stress_ptr->vec().zero_entries();
    eq_fields_->output_von_mises_stress_ptr->vec().zero_entries();
    eq_fields_->output_mean_stress_ptr->vec().zero_entries();
	eq_fields_->output_cross_section_ptr->vec().zero_entries();
	eq_fields_->output_div_ptr->vec().zero_entries();
	eq_fields_->output_field_ptr->vec().local_to_ghost_end();

    // compute new output fields depending on solution (stress, divergence etc.)
    output_fields_assembly_->assemble(eq_data_->dh_);

    // finish assembly of vectors
    eq_fields_->output_stress_ptr->vec().assembly_begin();
    eq_fields_->output_stress_ptr->vec().assembly_end();
    eq_fields_->output_von_mises_stress_ptr->vec().assembly_begin();
    eq_fields_->output_von_mises_stress_ptr->vec().assembly_end();
    eq_fields_->output_mean_stress_ptr->vec().assembly_begin();
    eq_fields_->output_mean_stress_ptr->vec().assembly_end();
    eq_fields_->output_cross_section_ptr->vec().assembly_begin();
    eq_fields_->output_cross_section_ptr->vec().assembly_end();
    eq_fields_->output_div_ptr->vec().assembly_begin();
    eq_fields_->output_div_ptr->vec().assembly_end();

    // update ghost values of computed fields
    eq_fields_->output_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_von_mises_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_mean_stress_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_cross_section_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_div_ptr->vec().local_to_ghost_begin();
    eq_fields_->output_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_von_mises_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_mean_stress_ptr->vec().local_to_ghost_end();
    eq_fields_->output_cross_section_ptr->vec().local_to_ghost_end();
    eq_fields_->output_div_ptr->vec().local_to_ghost_end();
}




void Elasticity::zero_time_step()
{
	START_TIMER("Mechanics zero time step");
	eq_fields_->mark_input_times( *time_ );
	eq_fields_->set_time(time_->step(), LimitSide::right);
	std::stringstream ss; // print warning message with table of uninitialized fields
	if ( FieldCommon::print_message_table(ss, "mechanics") ) {
		WarningOut() << ss.str();
	}

    preallocate();
    

    // after preallocation we assemble the matrices and vectors required for balance of forces
//     for (auto subst_idx : eq_data_->balance_idx_)
//         balance_->calculate_instant(subst_idx, eq_data_->ls->get_solution());
    
//     update_solution();
    eq_data_->ls->start_add_assembly();
    MatSetOption(*eq_data_->ls->get_matrix(), MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    eq_data_->ls->mat_zero_entries();
    eq_data_->ls->rhs_zero_entries();
    stiffness_assembly_->assemble(eq_data_->dh_);
    rhs_assembly_->assemble(eq_data_->dh_);
	eq_data_->ls->finish_assembly();
	LinSys::SolveInfo si = eq_data_->ls->solve();
    ResidualSplit residual = compute_residual_split_by_dirichlet_dofs(eq_data_.get());
	MessageOut().fmt(
            "[mech solver] lin. it: {}, reason: {}, residual: {}, residual_free: {}, residual_support: {}, "
            "residual_free_interface: {}, residual_free_interior: {}, "
            "residual_filtered: {}, residual_filtered_free: {}, residual_filtered_support: {}, "
            "residual_filtered_free_interface: {}, residual_filtered_free_interior: {}, "
            "solution_interface_jump_norm: {}, solution_interface_jump_max: {}, solution_interface_jump_nnz: {}, "
            "lagrangian_residual: {}, lagrangian_residual_form: {}\n",
            si.n_iterations, si.converged_reason,
            residual.norm, residual.free_norm, residual.support_norm,
            residual.free_interface_norm, residual.free_interior_norm,
            residual.filtered_norm, residual.filtered_free_norm, residual.filtered_support_norm,
            residual.filtered_free_interface_norm, residual.filtered_free_interior_norm,
            residual.solution_interface_jump_norm, residual.solution_interface_jump_max,
            residual.solution_interface_jump_nnz,
            eq_data_->ls->get_lagrangian_residual_norm(),
            eq_data_->ls->get_lagrangian_residual_name());
    output_data();
    END_TIMER("Mechanics zero time step");
}



void Elasticity::preallocate()
{
    // preallocate system matrix
	eq_data_->ls->start_allocation();
	stiffness_assembly_->assemble(eq_data_->dh_);
	rhs_assembly_->assemble(eq_data_->dh_);

    eq_data_->dirichlet_dofs.clear();
    eq_data_->dirichlet_coefs.clear();
    eq_data_->dirichlet_values.clear();
    eq_data_->dirichlet_row_starts.clear();
    eq_data_->dirichlet_row_starts.push_back(0);
    dirichlet_assembly_->assemble(eq_data_->dh_);

    if (eq_data_->dirichlet_by_eq) {
        // preallocate dirichlet constraint matrix and vector
        unsigned int n_rows = eq_data_->dirichlet_values.size();
        PetscErrorCode ierr;
        ISLocalToGlobalMapping l2g_rows, l2g_cols;
        Distribution dir_row_ds(n_rows, PETSC_COMM_WORLD);
        vector<LongIdx> rows_l2g_indices;
        for (unsigned int i=dir_row_ds.begin(); i<dir_row_ds.end(); i++)
            rows_l2g_indices.push_back(i);
        ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, n_rows,
            rows_l2g_indices.data(), PETSC_USE_POINTER, &l2g_rows); CHKERRV(ierr);
        ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, eq_data_->dh_->get_local_to_global_map().size(),
            eq_data_->dh_->get_local_to_global_map().data(), PETSC_USE_POINTER, &l2g_cols); CHKERRV(ierr);
        ierr = MatCreateIS(PETSC_COMM_WORLD, 1, n_rows, eq_data_->dh_->lsize(), PETSC_DETERMINE, PETSC_DETERMINE,
                                  l2g_rows, l2g_cols, &eq_data_->dirichlet_matrix); CHKERRV( ierr );
        PetscInt *on_nz, *off_nz;
        on_nz  = new PetscInt[ n_rows ];
        off_nz = new PetscInt[ n_rows ];
        for (unsigned int i=0; i<n_rows; i++) {
            on_nz[i] = eq_data_->dirichlet_row_starts[i+1] - eq_data_->dirichlet_row_starts[i];
            off_nz[i] = 0;
        }
        ierr = MatISSetPreallocation(eq_data_->dirichlet_matrix, 0, on_nz, 0, off_nz);
        delete[] on_nz;
        delete[] off_nz;
        VecCreateMPI(PETSC_COMM_WORLD, n_rows, PETSC_DECIDE, &eq_data_->dirichlet_vec);

        // assemble dirichlet matrix and vector
        MatZeroEntries(eq_data_->dirichlet_matrix);
        VecZeroEntries(eq_data_->dirichlet_vec);
        for (unsigned int i=0; i<n_rows; i++) {
            for (unsigned int j=eq_data_->dirichlet_row_starts[i]; j<eq_data_->dirichlet_row_starts[i+1]; j++)
                MatSetValueLocal(eq_data_->dirichlet_matrix, i, eq_data_->dirichlet_dofs[j], eq_data_->dirichlet_coefs[j], INSERT_VALUES);
            VecSetValueLocal(eq_data_->dirichlet_vec, i, eq_data_->dirichlet_values[i], INSERT_VALUES);
        }
        MatAssemblyBegin(eq_data_->dirichlet_matrix, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(eq_data_->dirichlet_matrix, MAT_FINAL_ASSEMBLY);
        VecAssemblyBegin(eq_data_->dirichlet_vec);
        VecAssemblyEnd(eq_data_->dirichlet_vec);
        eq_data_->ls->set_equality(eq_data_->dirichlet_matrix, eq_data_->dirichlet_vec);
    }

    if (has_contact_)
        assemble_constraint_matrix();
}



void Elasticity::next_time()
{
    time_->next_time();
    time_->view("MECH");
    
}



void Elasticity::solve_linear_system()
{
    START_TIMER("Mechanics step");
    START_TIMER("data reinit");
    eq_fields_->set_time(time_->step(), LimitSide::right);
    END_TIMER("data reinit");
    
    // assemble stiffness matrix
    if (eq_data_->ls->get_matrix() == NULL
        || eq_fields_->subset(FieldFlag::in_main_matrix).changed())
    {
        DebugOut() << "Mechanics: Assembling matrix.\n";
        eq_data_->ls->start_add_assembly();
        eq_data_->ls->mat_zero_entries();
        stiffness_assembly_->assemble(eq_data_->dh_);
        eq_data_->ls->finish_assembly();
    }

    // assemble right hand side (due to sources and boundary conditions)
    if (eq_data_->ls->get_rhs() == NULL
        || eq_fields_->subset(FieldFlag::in_rhs).changed())
    {
        DebugOut() << "Mechanics: Assembling right hand side.\n";
        eq_data_->ls->start_add_assembly();
        eq_data_->ls->rhs_zero_entries();
        rhs_assembly_->assemble(eq_data_->dh_);
        eq_data_->ls->finish_assembly();
    }

    START_TIMER("solve");
    LinSys::SolveInfo si = eq_data_->ls->solve();
    ResidualSplit residual = compute_residual_split_by_dirichlet_dofs(eq_data_.get());
    MessageOut().fmt(
            "[mech solver] lin. it: {}, reason: {}, residual: {}, residual_free: {}, residual_support: {}, "
            "residual_free_interface: {}, residual_free_interior: {}, "
            "residual_filtered: {}, residual_filtered_free: {}, residual_filtered_support: {}, "
            "residual_filtered_free_interface: {}, residual_filtered_free_interior: {}, "
            "solution_interface_jump_norm: {}, solution_interface_jump_max: {}, solution_interface_jump_nnz: {}, "
            "lagrangian_residual: {}, lagrangian_residual_form: {}\n",
            si.n_iterations, si.converged_reason,
            residual.norm, residual.free_norm, residual.support_norm,
            residual.free_interface_norm, residual.free_interior_norm,
            residual.filtered_norm, residual.filtered_free_norm, residual.filtered_support_norm,
            residual.filtered_free_interface_norm, residual.filtered_free_interior_norm,
            residual.solution_interface_jump_norm, residual.solution_interface_jump_max,
            residual.solution_interface_jump_nnz,
            eq_data_->ls->get_lagrangian_residual_norm(),
            eq_data_->ls->get_lagrangian_residual_name());
    END_TIMER("solve");
    END_TIMER("Mechanics step");
}




void Elasticity::output_data()
{
    START_TIMER("Mechanics output");

    // gather the solution from all processors
    eq_fields_->output_fields.set_time( this->time().step(), LimitSide::left);
    //if (eq_fields_->output_fields.is_field_output_time(eq_fields_->output_field, this->time().step()) )
    update_output_fields();
    eq_fields_->output_fields.output(this->time().step());

    if (view_matrices) {
        eq_data_->dh_->view_dof_to_node_map("mechanics" + boost::lexical_cast<std::string>(this->time().step()));
        eq_data_->ls->view("permon");
    }

//     START_TIMER("Mechanics-balance");
//     balance_->calculate_instant(subst_idx, eq_data_->ls->get_solution());
//     balance_->output();
//     END_TIMER("Mechanics-balance");

    END_TIMER("Mechanics output");
}




void Elasticity::calculate_cumulative_balance()
{
//     if (balance_->cumulative())
//     {
//         balance_->calculate_cumulative(subst_idx, eq_data_->ls->get_solution());
//     }
}


void Elasticity::assemble_constraint_matrix()
{
    MatZeroEntries(eq_data_->constraint_matrix);
    VecZeroEntries(eq_data_->constraint_vec);
    constraint_assembly_->assemble(eq_data_->dh_);
    MatAssemblyBegin(eq_data_->constraint_matrix, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(eq_data_->constraint_matrix, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(eq_data_->constraint_vec);
    VecAssemblyEnd(eq_data_->constraint_vec);
}

const Vec &Elasticity::get_solution()
{
    return eq_data_->ls->get_solution();
}
