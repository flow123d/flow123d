#ifndef TRANSPORT_OPERATOR_SPLITTING_HH_
#define TRANSPORT_OPERATOR_SPLITTING_HH_

#include "coupling/equation.hh"

#include <limits>

#include "io/output_time.hh"
//#include "flow/darcy_flow_mh.hh"
#include "flow/mh_dofhandler.hh"
#include "fields/field_algo_base.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "transport/substance.hh"
#include "transport/advection_process_base.hh"


/// external types:
class Mesh;
class ReactionTerm;
class ConvectionTransport;
class Semchem_interface;
class Balance;







/**
 * Abstract interface class for implementations of transport equation within TransportOperatorSplitting.
 */
class ConcentrationTransportBase : public EquationBase {
public:

    /**
     * Constructor.
     */
    ConcentrationTransportBase(Mesh &init_mesh, const Input::Record in_rec)
	: EquationBase(init_mesh, in_rec) {};


	/// Common specification of the input record for secondary equations.
    static Input::Type::AbstractRecord & get_input_type();


    virtual void initialize() = 0;

    /**
     * Set time interval which is considered as one time step by TransportOperatorSplitting.
     * In particular the velocity field dosn't change over this interval.
     *
     * Dependencies:
     *
     * velocity, porosity -> matrix, source_vector
     * matrix -> time_step
     *
     * data_read_times -> time_step (not necessary if we won't stick to jump times)
     * data -> source_vector
     * time_step -> scaling
     *
     *
     *
     */
    virtual void set_target_time(double target_time) = 0;

    /**
     * Use Balance object from upstream equation (e.g. in various couplings) instead of own instance.
     */
    virtual void set_balance_object(boost::shared_ptr<Balance> balance) = 0;

    virtual const vector<unsigned int> &get_subst_idx() = 0;

    /**
     * Calculate quantities necessary for cumulative balance (over time).
     * This method is called at each (sub)iteration of the time loop.
     */
    virtual void calculate_cumulative_balance() = 0;

    /**
     * Calculate instant quantities at output times.
     */
    virtual void calculate_instant_balance() = 0;


    virtual std::shared_ptr<OutputTime> &output_stream() = 0;

	virtual double **get_concentration_matrix() = 0;

	/// Return PETSc vector with solution for sbi-th substance.
	virtual const Vec &get_solution(unsigned int sbi) = 0;

	virtual void get_par_info(int * &el_4_loc, Distribution * &el_ds) = 0;

	virtual int *get_row_4_el() = 0;

    virtual void set_velocity_field(const MH_DofHandler &dh) = 0;

    /// Returns number of trnasported substances.
    virtual unsigned int n_substances() = 0;

    /// Returns reference to the vector of substnace names.
    virtual SubstanceList &substances() = 0;


protected:

    TimeGovernor &time() { return *time_; }

    friend class TransportOperatorSplitting;


};





/**
 * Class with fields that are common to all transport models.
 */
class TransportEqData : public FieldSet {
public:

	TransportEqData();
	inline virtual ~TransportEqData() {};

	/// Mobile porosity
	Field<3, FieldValue<3>::Scalar> porosity;

	/// Pointer to DarcyFlow field cross_section
	Field<3, FieldValue<3>::Scalar > cross_section;

	/// Concentration sources - density of substance source, only positive part is used.
	Field<3, FieldValue<3>::Vector> sources_density;
	/// Concentration sources - Robin type, in_flux = sources_sigma * (sources_conc - mobile_conc)
	Field<3, FieldValue<3>::Vector> sources_sigma;
	Field<3, FieldValue<3>::Vector> sources_conc;

};




/**
 * @brief Empty transport class.
 */
class TransportNothing : public AdvectionProcessBase {
public:
    inline TransportNothing(Mesh &mesh_in)
    : AdvectionProcessBase(mesh_in, Input::Record() )

    {
        // make module solved for ever
        time_=new TimeGovernor();
        time_->next_time();
    };

    inline virtual ~TransportNothing()
    {}

    inline void set_velocity_field(const MH_DofHandler &dh) override {};

    inline virtual void output_data() override {};


};



/**
 * @brief Coupling of a transport model with a reaction model by operator splitting.
 *
 * Outline:
 * Transport model is any descendant of TransportBase (even TransportOperatorSplitting itself). This
 * should perform the transport possibly with diffusion and usually without coupling between substances and phases.
 *
 * Reaction is any descendant of the ReactionBase class. This represents reactions in general way of any coupling that
 * happens between substances and phases on one element or more generally on one DoF.
 */

class TransportOperatorSplitting : public AdvectionProcessBase {
public:
    /**
     * @brief Declare input record type for the equation TransportOperatorSplittiong.
     *
     * TODO: The question is if this should be a general coupling class
     * (e.g. allow coupling TranportDG with reactions even if it is not good idea for numerical reasons.)
     * To make this a coupling class we should modify all main input files for transport problems.
     */
    static const Input::Type::Record & get_input_type();

    /// Constructor.
    TransportOperatorSplitting(Mesh &init_mesh, const Input::Record in_rec);
    /// Destructor.
    virtual ~TransportOperatorSplitting();

    virtual void set_velocity_field(const MH_DofHandler &dh) override;

    void zero_time_step() override;
    void update_solution() override;

    void compute_until_save_time();
    void compute_internal_step();
    void output_data() override;

   

private:
    /// Registrar of class to factory
    static const int registrar;

    std::shared_ptr<ConcentrationTransportBase> convection;
    std::shared_ptr<ReactionTerm> reaction;

    double *** semchem_conc_ptr;   //dumb 3-dim array (for phases, which are not supported any more) 
    Semchem_interface *Semchem_reactions;

};





#endif // TRANSPORT_OPERATOR_SPLITTING_HH_
