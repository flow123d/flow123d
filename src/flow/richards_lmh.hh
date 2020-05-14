/*
 * richards_lmh.hh
 *
 *  Created on: Sep 16, 2015
 *      Author: jb
 */

#ifndef SRC_FLOW_RICHARDS_LMH_HH_
#define SRC_FLOW_RICHARDS_LMH_HH_

#include <boost/exception/info.hpp>  // for operator<<, error_info::error_in...
#include <memory>                    // for shared_ptr
#include "fields/field.hh"           // for Field
#include "fields/field_fe.hh"           // for FieldFE
#include "fields/field_values.hh"    // for FieldValue<>::Scalar, FieldValue
#include "la/vector_mpi.hh"          // for VectorMPI
#include "flow/darcy_flow_lmh.hh"    // for DarcyLMH, DarcyLMH::EqData
#include "flow/darcy_flow_mh_output.hh" // for DarcyFlowMHOutput
#include "input/type_base.hh"        // for Array
#include "input/type_generic.hh"     // for Instance

class Mesh;
class SoilModelBase;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}

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
class RichardsLMH : public DarcyLMH
{
public:
    /// Class with all fields used in the equation DarcyFlow.
    /// This is common to all implementations since this provides interface
    /// to this equation for possible coupling.
    class EqData : public DarcyLMH::EqData {
    public:
        EqData();
        // input fields
        Field<3, FieldValue<3>::Scalar > water_content_saturated;   // corresponds to the porosity (theta_s = Vw/V = porosity)
        Field<3, FieldValue<3>::Scalar > water_content_residual;
        Field<3, FieldValue<3>::Scalar > genuchten_p_head_scale;
        Field<3, FieldValue<3>::Scalar > genuchten_n_exponent;

        //output fields
        Field<3, FieldValue<3>::Scalar > water_content;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> water_content_ptr;

        Field<3, FieldValue<3>::Scalar > conductivity_richards;
//         FieldFE<3, FieldValue<3>::Scalar > conductivity_richards;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar>> conductivity_ptr;

        // Auxiliary assembly fields.
        VectorMPI water_content_previous_time;
        VectorMPI capacity;

        // This is necessary in the assembly
        // TODO: store time information in the field set and in fields, is it ok also for more complex discretization methods?
        double time_step_;
        std::shared_ptr<SoilModelBase> soil_model_;
    };

    RichardsLMH(Mesh &mesh, const Input::Record in_rec, TimeGovernor *tm = nullptr);

    static const Input::Type::Record & get_input_type();
    
    void accept_time_step() override;
    
protected:
    /// Registrar of class to factory
    static const int registrar;

    bool zero_time_term(bool time_global=false) override;

    void initialize_specific() override;

    void initial_condition_postprocess() override;
    void assembly_linear_system() override;
private:

    std::shared_ptr<EqData> data_;
};




#endif /* SRC_FLOW_RICHARDS_LMH_HH_ */
