#ifndef TRANSPORT_OPERATOR_SPLITTING_HH_
#define TRANSPORT_OPERATOR_SPLITTING_HH_

#include "coupling/equation.hh"

#include <limits>
#include "io/output.h"
#include "flow/mh_dofhandler.hh"
#include "fields/field_base.hh"
#include "fields/field_values.hh"


/// external types:
//class LinSys;
//struct Solver;
class Mesh;
//class SchurComplement;
//class Distribution;
//class SparseGraph;

class Reaction;
class Linear_reaction;
//class Pade_approximant;
class Sorption;
class Semchem_interface;
class ConvectionTransport;


/**
 * @brief Specification of transport model interface.
 *
 * Here one has to specify methods for setting or getting data particular to
 * transport equations.
 */
class TransportBase : public EquationBase {
public:

	class TransportEqData : public EqDataBase {
	public:

		TransportEqData(const std::string& eq_name);
		virtual ~TransportEqData() {};

		Field<3, FieldValue<3>::Vector> init_conc; ///< Initial concentrations.
		Field<3, FieldValue<3>::Scalar> por_m;     ///< Mobile porosity

		/**
		 * Boundary conditions (Dirichlet) for concentrations.
		 * They are applied only in water inflow.
		 */
		BCField<3, FieldValue<3>::Vector> bc_conc;

		/// Pointer to DarcyFlow field cross_section
		Field<3, FieldValue<3>::Scalar > *cross_section;

		/// Concentration sources
		Field<3, FieldValue<3>::Vector> sources_density;
		Field<3, FieldValue<3>::Vector> sources_sigma;
		Field<3, FieldValue<3>::Vector> sources_conc;

	};

    TransportBase(Mesh &mesh, const Input::Record in_rec);
    virtual ~TransportBase();


    virtual TransportEqData *get_data() = 0;


    /**
     * This method takes sequential PETSc vector of side velocities and update
     * transport matrix. The ordering is same as ordering of sides in the mesh.
     *
     * TODO: We should pass whole velocity field object (description of base functions and dof numbering) and vector.
     */
    virtual void set_velocity_field(const MH_DofHandler &dh) {
        mh_dh=&dh;
    }
    virtual void output_data() =0;

    /**
     * Calculate mass balance: flux through boundary and volume sources
     */
    void mass_balance();

    virtual unsigned int n_substances() = 0;

    virtual vector<string> &substance_names() = 0;

    static Input::Type::AbstractRecord input_type;
    static Input::Type::Record input_type_output_record;

    const MH_DofHandler *mh_dh;

protected:

    virtual void calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance) = 0;
    virtual void calc_elem_sources(vector<vector<double> > &src_balance) = 0;

    FILE *balance_output_file;

};



/**
 * @brief Empty transport class.
 */
class TransportNothing : public TransportBase {
public:
   TransportNothing(Mesh &mesh_in)
    : TransportBase(mesh_in, Input::Record() )
    {
        // make module solved for ever
        time_=new TimeGovernor();
        time_->next_time();
    };

    virtual void get_solution_vector(double * &vector, unsigned int &size) {
        ASSERT( 0 , "Empty transport class do not provide solution!");
    }

    virtual void get_parallel_solution_vector(Vec &vector) {
        ASSERT( 0 , "Empty transport class do not provide solution!");
    };

    virtual void set_velocity_field(Vec &velocity_field) {};

    virtual void output_data() {};

    virtual TransportEqData *get_data() { return 0; };

    unsigned int n_substances() { return 0; };

    vector<string> &substance_names() {};


private:

    void calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance) {};
    void calc_elem_sources(vector<vector<double> > &src_balance) {};
};



/**
 * @brief Reaction transport implemented by operator splitting.
 */

class TransportOperatorSplitting : public TransportBase {
public:

	class EqData : public TransportBase::TransportEqData {
	public:

		EqData();
		virtual ~EqData() {};

        RegionSet read_boundary_list_item(Input::Record rec);

		Field<3, FieldValue<3>::Scalar> por_imm;   ///< Immobile porosity
		Field<3, FieldValue<3>::Vector> alpha;	   ///< Coefficients of non-equilibrium exchange
// TODO: sorp_type should be IntVector
		Field<3, FieldValue<3>::Vector> sorp_type;///< Type of sorption for each substance
		Field<3, FieldValue<3>::Vector> sorp_coef0;///< Coefficient of sorption for each substance
		Field<3, FieldValue<3>::Vector> sorp_coef1;///< Coefficient of sorption for each substance
		Field<3, FieldValue<3>::Scalar> phi;       ///< solid / solid mobile

	};

    TransportOperatorSplitting(Mesh &init_mesh, const Input::Record &in_rec);
    virtual ~TransportOperatorSplitting();

    /**
     * @brief Declare input record type for the equation TransportOperatorSplittiong.
     *
     * TODO: The question is if this should be a general coupling class
     * (e.g. allow coupling TranportDG with reactions even if it is not good idea for numerical reasons.)
     * To make this a coupling class we should modify all main input files for transport problems.
     *
     */
    static Input::Type::Record input_type;

    virtual void set_velocity_field(const MH_DofHandler &dh);
    virtual void update_solution();
    //virtual void compute_one_step();
    //virtual void compute_until();
    void compute_internal_step();
    void output_data();
    virtual void get_parallel_solution_vector(Vec &vc);
    virtual void get_solution_vector(double* &vector, unsigned int &size);
    void compute_until_save_time();

    unsigned int n_substances();
    vector<string> &substance_names();


   
    /**
     * @brief Sets pointer to data of other equations.
     * TODO: there should be also passed the sigma parameter between dimensions
     * @param cross_section is pointer to cross_section data of Darcy flow equation
     */
    void set_eq_data(Field<3, FieldValue<3>::Scalar > *cross_section);

    virtual EqData *get_data() { return &data; };



private:

    void calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance);
    void calc_elem_sources(vector<vector<double> > &src_balance);

    EqData data;

    ConvectionTransport *convection;
    Reaction *decayRad; //Linear_reaction *decayRad; //Reaction *decayRad;
    Sorption *sorptions;
    Semchem_interface *Semchem_reactions;
    //int steps;
    OutputTime *field_output;

    TimeMark::Type output_mark_type;

};



#endif // TRANSPORT_OPERATOR_SPLITTING_HH_
