/*
 * richards_lmh.hh
 *
 *  Created on: Sep 16, 2015
 *      Author: jb
 */

#ifndef SRC_FLOW_RICHARDS_LMH_HH_
#define SRC_FLOW_RICHARDS_LMH_HH_

#include "flow/darcy_flow_mh.hh"
#include "fields/vec_seq_double.hh"

/**
 * @brief Edge lumped mixed-hybrid solution of unsteady Darcy flow.
 *
 * The time term and sources are evenly distributed form an element to its edges.
 * This applies directly to the second Schur complement. After this system for pressure traces is solved we reconstruct pressures and side flows as follows:
 *
 * -# Element pressure is  average of edge pressure. This is in fact same as the MH for steady case so we let SchurComplement class do its job.
 *
 * -# We let SchurComplement to reconstruct fluxes and then account time term and sources which are evenly distributed from an element to its sides.
 *    It can be proved, that this keeps continuity of the fluxes over the edges.
 *
 * This lumping technique preserves discrete maximum principle for any time step provided one use acute mesh. But in practice even worse meshes are tractable.
 *
 * Ideas how to unify steady and unsteady flow:
 * zero_time_step:
 *
 * -# Set initial time.
 * -# Read initial condition. Reconstruct pressures.
 * -# Assembly system (possibly in matrix free way).
 * -# Reconstruct velocities (schur complement resolve).
 * -# Solve iteratively as regular time step if an input flag "steady_initial_time" is set.
 * -# (Detect that there is no time term. I such case use arbitrary long time step up to next change of data.
 *     Some kind of time step estimator would be nice.
 *
 * update solution:
 *
 * -# Move to the next time.
 * -# Update fields
 * -# Nonlinear solve.
 * -# In case of slow convergence, use shorter time-step, within estimated limits. Otherwise there is a different problem.
 */
class RichardsLMH : public DarcyMH
{
public:
    /// Class with all fields used in the equation DarcyFlow.
    /// This is common to all implementations since this provides interface
    /// to this equation for possible coupling.
    class EqData : public DarcyMH::EqData {
    public:
        EqData();
        // input fields
        Field<3, FieldValue<3>::Scalar > water_content_saturated;
        Field<3, FieldValue<3>::Scalar > water_content_residual;
        Field<3, FieldValue<3>::Scalar > genuchten_p_head_scale;
        Field<3, FieldValue<3>::Scalar > genuchten_n_exponent;

        //output fields
    };

    RichardsLMH(Mesh &mesh, const Input::Record in_rec);

    static const Input::Type::Record & get_input_type();
protected:
    /// Registrar of class to factory
    static const int registrar;

    void initialize_specific() override;
    //void local_assembly_specific(LocalAssemblyData &local_data) override;
    void assembly_source_term() override;

    void read_initial_condition() override;
    void assembly_linear_system() override;
    void setup_time_term();
    void postprocess() override;
private:

    std::shared_ptr<EqData> data_;
    /// PETSC scatter from the solution vector to the parallel edge vector with ghost values.
    VecScatter solution_2_edge_scatter_;

    /*
    Vec steady_diagonal;
    Vec steady_rhs;
    Vec new_diagonal;
    Vec previous_solution;
*/
    VectorMPI phead_edge_;
    //Vec time_term;
};




#endif /* SRC_FLOW_RICHARDS_LMH_HH_ */
