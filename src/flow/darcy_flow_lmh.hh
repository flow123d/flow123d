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
 * @file    darcy_flow_lmh.hh
 * @brief   Lumped mixed-hybrid model of linear Darcy flow, possibly unsteady.
 * @author  Jan Brezina
 *
 * Main object for mixed-hybrid discretization of the linear elliptic PDE (Laplace)
 * on  a multidimensional domain. Discretization of saturated Darcy flow using
 * RT0 approximation for the velocity
 */

/*
 * list of files dependent on this one:
 *
 * posprocess.cc
 * problem.cc
 * main.hh
 * transport.cc
 */


#ifndef DARCY_FLOW_LMH_HH
#define DARCY_FLOW_LMH_HH

#include <petscmat.h>                           // for Mat
#include <string.h>                             // for memcpy
#include <algorithm>                            // for swap
#include <boost/exception/info.hpp>             // for operator<<, error_inf...
#include <memory>                               // for shared_ptr, allocator...
#include <new>                                  // for operator new[]
#include <string>                               // for string, operator<<
#include <vector>                               // for vector, vector<>::con...
#include <armadillo>
#include "fields/bc_field.hh"                   // for BCField
#include "fields/field.hh"                      // for Field
#include "fields/field_fe.hh"                      // for FieldFE
#include "fields/field_set.hh"                  // for FieldSet
#include "fields/field_values.hh"               // for FieldValue<>::Scalar
#include "flow/darcy_flow_interface.hh"         // for DarcyFlowInterface
#include "input/input_exception.hh"             // for DECLARE_INPUT_EXCEPTION
#include "input/type_base.hh"                   // for Array
#include "input/type_generic.hh"                // for Instance
#include "mesh/mesh.h"                          // for Mesh
#include "petscvec.h"                           // for Vec, _p_Vec, VecScatter
#include "system/exceptions.hh"                 // for ExcStream, operator<<
#include "tools/time_governor.hh"               // for TimeGovernor
#include "la/vector_mpi.hh"                     // for VectorMPI

#include "flow/darcy_flow_mh.hh"                // for DarcyMH::EqData

class Balance;
class DarcyFlowMHOutput;
class Element;
class Intersection;
class LinSys;
// class LinSys_BDDC;
class LocalSystem;
namespace Input {
	class AbstractRecord;
	class Record;
	namespace Type {
		class Record;
		class Selection;
	}
}

template<int spacedim, class Value> class FieldAddPotential;
template<int spacedim, class Value> class FieldDivide;

/**
 * @brief Mixed-hybrid model of linear Darcy flow, possibly unsteady.
 *
 * Abstract class for various implementations of Darcy flow. In future there should be
 * one further level of abstraction for general time dependent problem.
 *
 * maybe TODO:
 * split compute_one_step to :
 * 1) prepare_next_timestep
 * 2) actualize_solution - this is for iterative nonlinear solvers
 *
 */


/**
 * @brief Mixed-hybrid of steady Darcy flow with sources and variable density.
 *
 * solve equations:
 * @f[
 *      q= -{\mathbf{K}} \nabla h -{\mathbf{K}} R \nabla z
 * @f]
 * @f[
 *      \mathrm{div} q = f
 * @f]
 *
 * where
 * - @f$ q @f$ is flux @f$[ms^{-1}]@f$ for 3d, @f$[m^2s^{-1}]@f$ for 2d and @f$[m^3s^{-1}]@f$ for 1d.
 * - @f$ \mathbf{K} @f$ is hydraulic tensor ( its orientation for 2d, 1d case is questionable )
 * - @f$ h = \frac{\pi}{\rho_0 g}+z @f$ is pressure head, @f$ \pi, \rho_0, g @f$ are the pressure, water density, and acceleration of gravity , respectively.
 *   Assumes gravity force acts counter to the direction of the @f$ z @f$ axis.
 * - @f$ R @f$ is destity or gravity variability coefficient. For density driven flow it should be
 * @f[
 *    R = \frac{\rho}{\rho_0} -1 = \rho_0^{-1}\sum_{i=1}^s c_i
 * @f]
 *   where @f$ c_i @f$ is concentration in @f$ kg m^{-3} @f$.
 *
 * The time key is optional, when not specified the equation is forced to steady regime. Using Steady TimeGovernor which have no dt constraints.
 *
 *
 * TODO:
 * Make solution regular field (need FeSeystem and parallel DofHandler for edge pressures), then remove get_solution_vector from
 * Equation interface.
 */
/**
 * Model for transition coefficients due to Martin, Jaffre, Roberts (see manual for full reference)
 *
 * TODO:
 * - how we can reuse field values computed during assembly
 *
 */

class DarcyLMH : public DarcyFlowInterface
{
public:
    TYPEDEF_ERR_INFO( EI_Reason, string);
    DECLARE_EXCEPTION(ExcSolverDiverge,
            << "Diverged nonlinear solver. Reason: " << EI_Reason::val
             );
    DECLARE_INPUT_EXCEPTION(ExcMissingTimeGovernor,
            << "Missing the key 'time', obligatory for the transient problems.");

    /** Class with all fields used in the equation DarcyFlow.
    * This is common to all implementations since this provides interface
    * to this equation for possible coupling.
    * 
    * This class is derived from DarcyMH::EqData especially due to the common output class DarcyFlowMHOutput.
    * This is the only dependence between DarcyMH and DarcyLMH classes.
    * It is also base class of RichardsLMH::EqData.
    * */
    class EqData : public DarcyMH::EqData {
    public:
        
        EqData();
        
        std::shared_ptr<SubDOFHandlerMultiDim> dh_p_;    ///< DOF handler represents DOFs of element pressure
        
        // Propagate test for the time term to the assembly.
        // This flag is necessary for switching BC to avoid setting zero neumann on the whole boundary in the steady case.
        bool use_steady_assembly_;
        
        // for time term assembly
        double time_step_;
        
        std::shared_ptr<LinSys> lin_sys_schur;  //< Linear system of the 2. Schur complement.
        VectorMPI p_edge_solution;               //< 2. Schur complement solution
        VectorMPI p_edge_solution_previous;      //< 2. Schur complement previous solution (iterative)
        VectorMPI p_edge_solution_previous_time; //< 2. Schur complement previous solution (time)

        std::map<LongIdx, LocalSystem> seepage_bc_systems;
    };

    /// Selection for enum MortarMethod.
    static const Input::Type::Selection & get_mh_mortar_selection();





    DarcyLMH(Mesh &mesh, const Input::Record in_rec, TimeGovernor *tm = nullptr);

    static const Input::Type::Record & type_field_descriptor();
    static const Input::Type::Record & get_input_type();

    double last_t() override {
        return time_->last_t();
    }

    std::shared_ptr< FieldFE<3, FieldValue<3>::VectorFixed> > get_velocity_field() override;

    void init_eq_data();
    void initialize() override;
    virtual void initialize_specific();
    void zero_time_step() override;
    void update_solution() override;

    /// Solve the problem without moving to next time and without output.
    void solve_time_step(bool output = true);
    
    /// postprocess velocity field (add sources)
    virtual void accept_time_step();
    virtual void postprocess();
    virtual void output_data() override;


    EqData &data() { return *data_; }

    /// Sets external storarivity field (coupling with other equation).
    void set_extra_storativity(const Field<3, FieldValue<3>::Scalar> &extra_stor)
    { data_->extra_storativity = extra_stor; }

    /// Sets external source field (coupling with other equation).
    void set_extra_source(const Field<3, FieldValue<3>::Scalar> &extra_src)
    { data_->extra_source = extra_src; }

    virtual ~DarcyLMH() override;


protected:
    /**
     * Returns true is the fields involved in the time term have values that makes the time term zero.
     * For time_global==true, it returns true if there are no field descriptors in the input list, so the
     * fields )of the time ter) have their default values for whole simulation.
     * If time_global==false (default), only the actual values are considered.
     */
    virtual bool zero_time_term(bool time_global=false);

    /// Solve method common to zero_time_step and update solution.
    void solve_nonlinear();

    /**
     * Create and preallocate MH linear system (including matrix, rhs and solution vectors)
     */
    void create_linear_system(Input::AbstractRecord rec);

    /**
     * Read initial condition into solution vector.
     * Must be called after create_linear_system.
     *
     * For the LMH scheme we have to be able to save edge pressures in order to
     * restart simulation or use results of one simulation as initial condition for other one.
     */
    void read_initial_condition();
    
    /**
     * In some circumstances, the intial condition must be processed.
     * It is called at the end of @p read_initial_condition().
     * This is used in Richards equation due the update of water content.
     */
    virtual void initial_condition_postprocess();
    
    /**
     * Allocates linear system matrix for MH.
     * TODO:
     * - use general preallocation methods in DofHandler
     */
    void allocate_mh_matrix();

    /**
     * Assembles linear system matrix for MH.
     * Element by element assembly is done using dim-template assembly class.
     * Assembles only steady part of the equation.
     * TODO:
     * - include time term - DONE
     * - add support for Robin type sources
     * - support for nonlinear solvers - assembly either residual vector, matrix, or both (using FADBAD++)
     */
    void assembly_mh_matrix(MultidimAssembly& assembler);
    
    void reconstruct_solution_from_schur(MultidimAssembly& assembler);

    /**
     * Assembly or update whole linear system.
     */
    virtual void assembly_linear_system();

//     void set_mesh_data_for_bddc(LinSys_BDDC * bddc_ls);
    /**
     * Return a norm of residual vector.
     * TODO: Introduce Equation::compute_residual() updating
     * residual field, standard part of EqData.
     */
    virtual double solution_precision() const;
    
    /// Print darcy flow matrix in matlab format into a file.
    void print_matlab_matrix(string matlab_file);

    /// Get vector of all DOF indices of given component (0..side, 1..element, 2..edge)
    std::vector<int> get_component_indices_vec(unsigned int component) const;

    /// Getter for the linear system of the 2. Schur complement.
    LinSys& lin_sys_schur()
    { return *(data_->lin_sys_schur); }

    std::shared_ptr<Balance> balance_;

    DarcyFlowMHOutput *output_object;

	int size;				    // global size of MH matrix

	bool data_changed_;

	// Setting of the nonlinear solver. TODO: Move to the solver class later on.
	double tolerance_;
	unsigned int min_n_it_;
	unsigned int max_n_it_;
	unsigned int nonlinear_iteration_; //< Actual number of completed nonlinear iterations, need to pass this information into assembly.

    // Temporary objects holding pointers to appropriate FieldFE
    // TODO remove after final fix of equations
    std::shared_ptr<FieldFE<3, FieldValue<3>::VectorFixed>> ele_flux_ptr;            ///< Field of flux in barycenter of every element.

	std::shared_ptr<EqData> data_;

    friend class DarcyFlowMHOutput;
    //friend class P0_CouplingAssembler;
    //friend class P1_CouplingAssembler;

private:
  /// Registrar of class to factory
  static const int registrar;
};

#endif  //DARCY_FLOW_LMH_HH
//-----------------------------------------------------------------------------
// vim: set cindent:

