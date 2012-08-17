#ifndef TRANSPORT_OPERATOR_SPLITTING_HH_
#define TRANSPORT_OPERATOR_SPLITTING_HH_

#include "coupling/equation.hh"

#include <limits>
#include "io/output.h"


/// external types:
//class LinSys;
//struct Solver;
class Mesh;
//class SchurComplement;
//class Distribution;
//class SparseGraph;

class Reaction;
class Linear_reaction;
class Semchem_interface;
class ConvectionTransport;
class MaterialDatabase;


/**
 * @brief Specification of transport model interface.
 *
 * Here one has to specify methods for setting or getting data particular to
 * transport equations.
 */
class TransportBase : public EquationBase{
public:
    TransportBase(TimeMarks &marks, Mesh &mesh, MaterialDatabase &mat_base, const Input::Record in_rec)
    : EquationBase(marks, mesh, mat_base, in_rec )
    {}

    static Input::Type::AbstractRecord &get_input_type();
    static Input::Type::Record & get_input_type_output_record();

    /**
     * This method takes sequantial PETSc vector of side velocities and update
     * transport matrix. The ordering is same as ordering of sides in the mesh.
     *
     * TODO: We should pass whole velocity field object (description of base functions and dof numbering) and vector.
     */
    virtual void set_velocity_field(Vec &velocity_vector) =0;
    virtual void output_data() =0;
};



/**
 * @brief Empty transport class.
 */
class TransportNothing : public TransportBase {
public:
    TransportNothing(TimeMarks &marks, Mesh &mesh_in, MaterialDatabase &mat_base_in)
    : TransportBase(marks, mesh_in, mat_base_in, Input::Record() )
    {
        // make module solved for ever
        time_=new TimeGovernor();
    };

    virtual void get_solution_vector(double * &vector, unsigned int &size) {
        ASSERT( 0 , "Empty transport class do not provide solution!");
    }

    virtual void get_parallel_solution_vector(Vec &vector) {
        ASSERT( 0 , "Empty transport class do not provide solution!");
    };

    virtual void set_velocity_field(Vec &velocity_field) {};

    virtual void output_data() {};
};



/**
 * @brief Reaction transport implemented by operator splitting.
 */

class TransportOperatorSplitting : public TransportBase {
public:
	TransportOperatorSplitting(TimeMarks &marks,  Mesh &init_mesh, MaterialDatabase &material_database, const Input::Record &in_rec);
    virtual ~TransportOperatorSplitting();

    /**
     * @brief Declare input record type for the equation TransportOperatorSplittiong.
     *
     * TODO: The question is if this should be a general coupling class
     * (e.g. allow coupling TranportDG with reactions even if it is not good idea for numerical reasons.)
     * To make this a coupling class we should modify all main input files for transport problems.
     *
     */
    static Input::Type::Record &get_input_type();

    virtual void set_velocity_field(Vec &velocity_vector);
	virtual void update_solution();
	void read_simulation_step(double sim_step);
	//virtual void compute_one_step();
	//virtual void compute_until();
	void compute_internal_step();
	void output_data();
	 virtual void get_parallel_solution_vector(Vec &vc);
	 virtual void get_solution_vector(double* &vector, unsigned int &size);
	 void compute_until_save_time();
protected:

private:

    ConvectionTransport *convection;
    Reaction *decayRad;
    Semchem_interface *Semchem_reactions;
    //int steps;
    OutputTime *field_output;

    TimeMark::Type output_mark_type;
};

#endif // TRANSPORT_OPERATOR_SPLITTING_HH_
