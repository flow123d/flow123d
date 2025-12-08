#ifndef DG_MOCKUP_HH_
#define DG_MOCKUP_HH_


#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include <rev_num.h>

#include "fem/eval_points.hh"
#include "fem/integral_acc.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "fields/generic_field.hh"
#include "fields/multi_field.hh"
#include "fields/bc_multi_field.hh"
#include "fields/equation_output.hh"
#include "fields/field_model.hh"
#include "fields/field_constant.hh"
#include "coupling/equation.hh"
#include "tools/unit_si.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/fe_p.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "la/linsys.hh"
#include "la/linsys_PETSC.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "flow/assembly_lmh.hh"
#include "flow/assembly_models.hh"
#include "darcy_mockup_assembly.hh"
#include "balance_null.hh"

class GenericAssemblyBase;
template< template<IntDim...> class DimAssembly> class GenericAssembly;



namespace it = Input::Type;

namespace equation_data {

class EqFields : public FieldSet {
public:
    /**
     * For compatibility with old BCD file we have to assign integer codes starting from 1.
     */
    enum BC_Type {
        none=0,
        dirichlet=1,
        total_flux=4,
        seepage=5,
        river=6
    };

    /// Return a Selection corresponding to enum BC_Type.
    static const Input::Type::Selection & get_bc_type_selection() {
    	return it::Selection("Flow_Darcy_BC_Type")
            .add_value(none, "none",
                "Homogeneous Neumann boundary condition\n(zero normal flux over the boundary).")
            .add_value(dirichlet, "dirichlet",
                "Dirichlet boundary condition. "
                "Specify the pressure head through the ``bc_pressure`` field "
                "or the piezometric head through the ``bc_piezo_head`` field.")
            .add_value(total_flux, "total_flux", "Flux boundary condition (combines Neumann and Robin type). "
                "Water inflow equal to (($ \\delta_d(q_d^N + \\sigma_d (h_d^R - h_d) )$)). "
                "Specify the water inflow by the ``bc_flux`` field, the transition coefficient by ``bc_robin_sigma`` "
                "and the reference pressure head or piezometric head through ``bc_pressure`` or ``bc_piezo_head`` respectively.")
            .add_value(seepage, "seepage",
                "Seepage face boundary condition. Pressure and inflow bounded from above. Boundary with potential seepage flow "
                "is described by the pair of inequalities: "
                "(($h_d \\le h_d^D$)) and (($ -\\boldsymbol q_d\\cdot\\boldsymbol n \\le \\delta q_d^N$)), where the equality holds in at least one of them. "
                "Caution: setting (($q_d^N$)) strictly negative "
                "may lead to an ill posed problem since a positive outflow is enforced. "
                "Parameters (($h_d^D$)) and (($q_d^N$)) are given by the fields ``bc_switch_pressure`` (or ``bc_switch_piezo_head``) and ``bc_flux`` respectively."
                )
            .add_value(river, "river",
                "River boundary condition. For the water level above the bedrock, (($H_d > H_d^S$)), the Robin boundary condition is used with the inflow given by: "
                "(( $ \\delta_d(q_d^N + \\sigma_d(H_d^D - H_d) )$)). For the water level under the bedrock, constant infiltration is used: "
                "(( $ \\delta_d(q_d^N + \\sigma_d(H_d^D - H_d^S) )$)). Parameters: ``bc_pressure``, ``bc_switch_pressure``, "
                " ``bc_sigma``, ``bc_flux``."
                )
            .close();
    }

    /// Creation of all fields.
    EqFields() {
        *this += field_ele_pressure.name("pressure_p0")
                 .units(UnitSI().m())
                 .flags(FieldFlag::equation_result)
                 .input_default("1.0")
                 .description("Pressure solution - P0 interpolation.");

        *this += field_edge_pressure.name("pressure_edge")
                 .units(UnitSI().m())
                 .flags(FieldFlag::input_copy)
                 .input_default("1.0")
                 .description("Pressure solution - Crouzeix-Raviart interpolation.");

        *this += field_ele_piezo_head.name("piezo_head_p0")
    	         .units(UnitSI().m())
                 .flags(FieldFlag::equation_result)
                 .input_default("1.0")
                 .description("Piezo head solution - P0 interpolation.");

    	*this += field_ele_velocity.name("velocity_p0")
    	         .units(UnitSI().m().s(-1))
                 .flags(FieldFlag::equation_result)
                 .input_default("1.0")
                 .description("Velocity solution - P0 interpolation.");

        *this += flux.name("flux")
    	         .units(UnitSI().m().s(-1))
                 .flags(FieldFlag::equation_result)
                 .input_default("1.0")
                 .description("Darcy flow flux.");

        *this += anisotropy.name("anisotropy")
                .description("Anisotropy of the conductivity tensor.")
                .input_default("1.0")
                .units( UnitSI::dimensionless() );

        *this += cross_section.name("cross_section")
                .description("Complement dimension parameter (cross section for 1D, thickness for 2D).")
                .input_default("1.0")
                .units( UnitSI().m(3).md() );

        *this += conductivity.name("conductivity")
                .description("Isotropic conductivity scalar.")
                .input_default("1.0")
                .units( UnitSI().m().s(-1) )
                .set_limits(0.0);

        *this += sigma.name("sigma")
                .description("Transition coefficient between dimensions.")
                .input_default("1.0")
                .units( UnitSI::dimensionless() );

        *this += water_source_density.name("water_source_density")
                .description("Water source density.")
                .input_default("0.0")
                .units( UnitSI().s(-1) );

        *this += bc_type.name("bc_type")
                .description("Boundary condition type.")
                .input_selection( get_bc_type_selection() )
                .input_default("\"none\"")
                .units( UnitSI::dimensionless() );

        *this += bc_pressure
                .disable_where(bc_type, {none, seepage} )
                .name("bc_pressure")
                .description("Prescribed pressure value on the boundary. Used for all values of ``bc_type`` except ``none`` and ``seepage``. "
                    "See documentation of ``bc_type`` for exact meaning of ``bc_pressure`` in individual boundary condition types.")
                .input_default("0.0")
                .units( UnitSI().m() );

        *this += bc_flux
                .disable_where(bc_type, {none, dirichlet} )
                .name("bc_flux")
                .description("Incoming water boundary flux. Used for bc_types : ``total_flux``, ``seepage``, ``river``.")
                .input_default("0.0")
                .units( UnitSI().m().s(-1) );

        *this += bc_robin_sigma
                .disable_where(bc_type, {none, dirichlet, seepage} )
                .name("bc_robin_sigma")
                .description("Conductivity coefficient in the ``total_flux`` or the ``river`` boundary condition type.")
                .input_default("0.0")
                .units( UnitSI().s(-1) );

        *this += bc_switch_pressure
                .disable_where(bc_type, {none, dirichlet, total_flux} )
                .name("bc_switch_pressure")
                .description("Critical switch pressure for ``seepage`` and ``river`` boundary conditions.")
                .input_default("0.0")
                .units( UnitSI().m() );


        //these are for unsteady
        *this += init_pressure.name("init_pressure")
                .description("Initial condition for pressure in time dependent problems.")
                .input_default("0.0")
                .units( UnitSI().m() );

        *this += storativity.name("storativity")
                .description("Storativity (in time dependent problems).")
                .input_default("0.0")
                .units( UnitSI().m(-1) );

        *this += extra_storativity.name("extra_storativity")
                .description("Storativity added from upstream equation.")
                .units( UnitSI().m(-1) )
                .input_default("0.0")
                .flags( input_copy );

        *this += extra_source.name("extra_water_source_density")
                .description("Water source density added from upstream equation.")
                .input_default("0.0")
                .units( UnitSI().s(-1) )
                .flags( input_copy );

        *this += gravity_field.name("gravity")
                .description("Gravity vector.")
                .input_default("0.0")
                .units( UnitSI::dimensionless() );

        *this += bc_gravity.name("bc_gravity")
                .description("Boundary gravity vector.")
                .input_default("0.0")
                .units( UnitSI::dimensionless() );

        *this += init_piezo_head.name("init_piezo_head")
    	         .units(UnitSI().m())
                 .input_default("0.0")
                 .description("Init piezo head.");

        *this += bc_piezo_head.name("bc_piezo_head")
    	         .units(UnitSI().m())
                 .input_default("0.0")
                 .description("Boundary piezo head.");

        *this += bc_switch_piezo_head.name("bc_switch_piezo_head")
    	         .units(UnitSI().m())
                 .input_default("0.0")
                 .description("Boundary switch piezo head.");

        *this += ref_pressure.name("ref_pressure")
                 .units(UnitSI().m())
    			 .input_default("0.0")
                 .flags(FieldFlag::equation_result)
                 .description("Precomputed pressure of l2 difference output.");

    	*this += ref_velocity.name("ref_velocity")
    	         .units(UnitSI().m().s(-1))
    			 .input_default("0.0")
                 .flags(FieldFlag::equation_result)
                 .description("Precomputed velocity of l2 difference output.");

        *this += ref_divergence.name("ref_divergence")
                 .units(UnitSI().m())
    			 .input_default("0.0")
                 .flags(FieldFlag::equation_result)
                 .description("Precomputed divergence of l2 difference output.");

        this->set_default_fieldset();
        //time_term_fields = this->subset({"storativity"});
        //main_matrix_fields = this->subset({"anisotropy", "conductivity", "cross_section", "sigma", "bc_type", "bc_robin_sigma"});
        //rhs_fields = this->subset({"water_source_density", "bc_pressure", "bc_flux"});
    }

    /// Return coords field
    FieldCoords &X() {
        return this->X_;
    }


    /// Initialize selected fields as FieldModels
    void init_field_models()
    {
		field_ele_velocity.set(Model<3, FieldValue<3>::VectorFixed>::create(fn_mh_velocity(), flux, cross_section), 0.0);

        field_ele_piezo_head.set(
                Model<3, FieldValue<3>::Scalar>::create(fn_mh_piezohead(), gravity_field, this->X(), field_ele_pressure),
                0.0
        );
    }

    /// Initialize selected fields as FieldConstants
    void init_field_constants(double ele_piezo_head, arma::vec3 ele_velocity)
    {
        {
    	    FieldValue<3>::Scalar f_value(ele_piezo_head);
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
            field_algo->set_value(f_value);
            field_ele_piezo_head.set(field_algo, 0.0);
        }
        {
    	    FieldValue<3>::VectorFixed f_value(ele_velocity);
            auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::VectorFixed>>();
            field_algo->set_value(f_value);
            field_ele_velocity.set(field_algo, 0.0);
        }
    }


    Field<3, FieldValue<3>::TensorFixed > anisotropy;
    Field<3, FieldValue<3>::Scalar > conductivity;
    Field<3, FieldValue<3>::Scalar > cross_section;
    Field<3, FieldValue<3>::Scalar > water_source_density;
    Field<3, FieldValue<3>::Scalar > sigma;

    BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
    BCField<3, FieldValue<3>::Scalar > bc_pressure;
    BCField<3, FieldValue<3>::Scalar > bc_flux;
    BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;
    BCField<3, FieldValue<3>::Scalar > bc_switch_pressure;

    Field<3, FieldValue<3>::Scalar > init_pressure;
    Field<3, FieldValue<3>::Scalar > storativity;
    Field<3, FieldValue<3>::Scalar > extra_storativity; /// Externally added storativity.
    Field<3, FieldValue<3>::Scalar > extra_source; /// Externally added water source.

    Field<3, FieldValue<3>::Scalar> field_ele_pressure;
    Field<3, FieldValue<3>::Scalar> field_ele_piezo_head;
    Field<3, FieldValue<3>::VectorFixed > field_ele_velocity;
    Field<3, FieldValue<3>::VectorFixed > flux;
    Field<3, FieldValue<3>::Scalar> field_edge_pressure;

    Field<3, FieldValue<3>::VectorFixed > gravity_field; /// Holds gravity vector acceptable in FieldModel
    BCField<3, FieldValue<3>::VectorFixed > bc_gravity; /// Same as previous but used in boundary fields
    Field<3, FieldValue<3>::Scalar> init_piezo_head;
    BCField<3, FieldValue<3>::Scalar> bc_piezo_head;
    BCField<3, FieldValue<3>::Scalar> bc_switch_piezo_head;

    Field<3, FieldValue<3>::Scalar> ref_pressure; /// Precompute l2 difference outputs
    Field<3, FieldValue<3>::VectorFixed> ref_velocity;
    Field<3, FieldValue<3>::Scalar> ref_divergence;
};

class EqData {
public:
    EqData() {}

    void init() {
        auto size = dh_p_->get_local_to_global_map().size();
        save_local_system_.resize(size);
        bc_fluxes_reconstruted.resize(size);
        loc_system_.resize(size);
        postprocess_solution_.resize(size);
    }

    void reset() {
        std::fill(save_local_system_.begin(), save_local_system_.end(), false);
        std::fill(bc_fluxes_reconstruted.begin(), bc_fluxes_reconstruted.end(), false);
    }

    /**
     * Gravity vector and constant shift of pressure potential. Used to convert piezometric head
     * to pressure head and vice versa.
     */
    arma::vec4 gravity_;
    arma::vec3 gravity_vec_;

    // Mirroring the following members of DarcyLMH:
    Mesh *mesh;
    std::shared_ptr<DOFHandlerMultiDim> dh_;         ///< full DOF handler represents DOFs of sides, elements and edges
    std::shared_ptr<SubDOFHandlerMultiDim> dh_cr_;   ///< DOF handler represents DOFs of edges
    std::shared_ptr<DOFHandlerMultiDim> dh_cr_disc_; ///< DOF handler represents DOFs of sides
    std::shared_ptr<SubDOFHandlerMultiDim> dh_p_;    ///< DOF handler represents DOFs of element pressure


    uint water_balance_idx;

    int is_linear;              ///< Hack fo BDDC solver.
    bool force_no_neumann_bc;       ///< auxiliary flag for switchting Dirichlet like BC

    /// Idicator of dirichlet or neumann type of switch boundary conditions.
    std::vector<char> bc_switch_dirichlet;

	VectorMPI full_solution;     //< full solution [vel,press,lambda] from 2. Schur complement

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

    /// Shared Balance object
	std::shared_ptr<BalanceNull> balance_;

	unsigned int nonlinear_iteration_; //< Actual number of completed nonlinear iterations, need to pass this information into assembly.

    std::vector<LocalSystem> loc_system_;
    std::vector<LocalConstraint> loc_constraint_;
    std::vector<arma::vec> postprocess_solution_;
    std::array<std::vector<unsigned int>, 3> loc_side_dofs;
    std::array<std::vector<unsigned int>, 3> loc_edge_dofs;
    std::array<unsigned int, 3> loc_ele_dof;

    std::vector<bool> save_local_system_;       ///< Flag for saving the local system. Currently used only in case of seepage BC.
    std::vector<bool> bc_fluxes_reconstruted;   ///< Flag indicating whether the fluxes for seepage BC has been reconstructed already.
    std::array<unsigned int, 3> schur_offset_;  ///< Index offset in the local system for the Schur complement (of dim = 1,2,3).
};

} // end of namespace equation_data


/// Test class
class DarcyMockupTest : public testing::Test {
public:
    template<unsigned int dim> using MHMatrixAssemblyLMHDim = MHMatrixAssemblyLMH<dim, equation_data::EqFields, equation_data::EqData>;
    template<unsigned int dim> using MHMatrixEvalFieldsDim = MHMatrixEvalFields<dim, equation_data::EqFields, equation_data::EqData>;

    DarcyMockupTest()
    {
		string root_dir=string(UNIT_TESTS_BIN_DIR) + "/coupling";
		string build = string(__DATE__) + ", " + string(__TIME__)
	            + " flags: (unknown compiler flags)";

        FilePath::set_io_dirs(".",root_dir,"",".");
        Profiler::instance();
        Profiler::instance()->set_program_info("Flow123d",
                string(FLOW123D_VERSION_NAME_), string(FLOW123D_GIT_BRANCH_), string(FLOW123D_GIT_REVISION_), build);
        Profiler::set_memory_monitoring(false);
    }

    ~DarcyMockupTest() {}

    /// Run assembly algorithms with different type of assembly and type of field
    void run_fullassembly_const(const string &eq_data_input, const std::string &mesh_file);
    void run_fullassembly_model(const string &eq_data_input, const std::string &mesh_file);
    void run_computelocal_const(const string &eq_data_input, const std::string &mesh_file);
    void run_computelocal_model(const string &eq_data_input, const std::string &mesh_file);
    void run_evalfields_const(const string &eq_data_input, const std::string &mesh_file);
    void run_evalfields_model(const string &eq_data_input, const std::string &mesh_file);

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}
};


/*******************************************************************************
 * Equivalent to TransportDG class
 */
template<template<IntDim...> class MhMatrix>
class DarcyMockup : public EquationBase {
public:

    const it::Record & type_field_descriptor() {
        std::string equation_name = "TestEquation";
        const it::Record &field_descriptor =
        it::Record(equation_name + "_Data",FieldCommon::field_descriptor_record_description(equation_name + "_Data") )
            .copy_keys( equation_data::EqFields().make_field_descriptor_type(equation_name + "_Data_aux") )
            .declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Boundary piezometric head for BC types: dirichlet, robin, and river." )
            .declare_key("bc_switch_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Boundary switch piezometric head for BC types: seepage, river." )
            .declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(),
                    "Initial condition for the pressure given as the piezometric head." )
            .close();
        return field_descriptor;
    }

    it::Record & get_input_type() {
        std::string equation_name = "TestEquation";
        it::Record ns_rec = Input::Type::Record("NonlinearSolver", "Non-linear solver settings.")
            .declare_key("linear_solver", LinSys::get_input_type(), it::Default("{}"),
                "Linear solver for MH problem.")
            .declare_key("tolerance", it::Double(0.0), it::Default("1E-6"),
                "Residual tolerance.")
            .declare_key("min_it", it::Integer(0), it::Default("1"),
                "Minimum number of iterations (linear solutions) to use.\nThis is usefull if the convergence criteria "
                "does not characterize your goal well enough so it converges prematurely, possibly even without a single linear solution."
                "If greater then 'max_it' the value is set to 'max_it'.")
            .declare_key("max_it", it::Integer(0), it::Default("100"),
                "Maximum number of iterations (linear solutions) of the non-linear solver.")
            .declare_key("converge_on_stagnation", it::Bool(), it::Default("false"),
                "If a stagnation of the nonlinear solver is detected the solver stops. "
                "A divergence is reported by default, forcing the end of the simulation. By setting this flag to 'true', the solver "
                "ends with convergence success on stagnation, but it reports warning about it.")
            .close();

        equation_data::EqFields eq_fields;

        return it::Record(equation_name, "Lumped Mixed-Hybrid solver for saturated Darcy flow.")
            .derive_from(DarcyFlowInterface::get_input_type())
            .copy_keys(EquationBase::record_template())
    		.copy_keys(EquationBase::user_fields_template(equation_name))
            .declare_key("gravity", it::Array(it::Double(), 3,3), it::Default("[ 0, 0, -1]"),
                    "Vector of the gravity force. Dimensionless.")
    		.declare_key("input_fields", it::Array( type_field_descriptor() ), it::Default::obligatory(),
                    "Input data for Darcy flow model.")
            .declare_key("nonlinear_solver", ns_rec, it::Default("{}"),
                    "Non-linear solver for MH problem.")
//            .declare_key("output_stream", OutputTime::get_input_type(), it::Default("{}"),
//                    "Output stream settings.\n Specify file format, precision etc.")

//            .declare_key("output", DarcyFlowMHOutput::get_input_type(eq_fields, equation_name()),
//                    IT::Default("{ \"fields\": [ \"pressure_p0\", \"velocity_p0\" ] }"),
//                    "Specification of output fields and output times.")
//            .declare_key("output_specific", DarcyFlowMHOutput::get_input_type_specific(), it::Default::optional(),
//                    "Output settings specific to Darcy flow model.\n"
//                    "Includes raw output and some experimental functionality.")
//            .declare_key("balance", Balance::get_input_type(), it::Default("{}"),
//                    "Settings for computing mass balance.")
//    	      .declare_key("mortar_method", get_mh_mortar_selection(), it::Default("\"None\""),
//        	    	"Method for coupling Darcy flow between dimensions on incompatible meshes. [Experimental]" )
    		.close();
    }

    DarcyMockup(bool use_linsys)
    : use_linsys_(use_linsys)
    {
        eq_data_ = make_shared<equation_data::EqData>();
        eq_fields_ = make_shared<equation_data::EqFields>();
        this->eq_fieldset_ = eq_fields_;
    }

    ~DarcyMockup() {
//        if (eq_data_->dif_coef.size() > 0) {
//            // initialize called
//
//            for (unsigned int i=0; i<eq_data_->n_substances(); i++)
//            {
//                if (eq_data_->ls != nullptr) {
//                    delete eq_data_->ls[i];
//                    delete eq_data_->ls_dt[i];
//                }
//
//                if (stiffness_matrix.size() > 0) {
//                    if (stiffness_matrix[i])
//                        chkerr(MatDestroy(&stiffness_matrix[i]));
//                    if (mass_matrix[i])
//                        chkerr(MatDestroy(&mass_matrix[i]));
//                    if (rhs[i])
//                        chkerr(VecDestroy(&rhs[i]));
//                    if (mass_vec[i])
//                        chkerr(VecDestroy(&mass_vec[i]));
//                    if (eq_data_->ret_vec[i])
//                        chkerr(VecDestroy(&eq_data_->ret_vec[i]));
//                }
//            }
//            if (eq_data_->ls != nullptr) {
//                delete[] eq_data_->ls;
//                delete[] eq_data_->ls_dt;
//                eq_data_->ls = nullptr;
//            }
//
//            if (mass_assembly_ != nullptr) {
//                delete mass_assembly_;
//                delete stiffness_assembly_;
//                delete sources_assembly_;
//            }
//        }
//
    }

    void create_and_set_mesh(const std::string &mesh_file) {
        std::string input_str = "{ mesh_file=\"" + mesh_file + "\" }";
        this->mesh_ = mesh_full_constructor(input_str);
        START_TIMER("n_mesh_elements");
        uint n_elements = this->mesh_->n_elements();
        // for(i =0; i<n_elements; i++) i+1;
        ADD_CALLS(n_elements);
        END_TIMER("n_mesh_elements");

        // Set up physical parameters.
        eq_fields_->set_mesh(*mesh_);

        { // init DOF handler for pressure fields
    // 		std::shared_ptr< FiniteElement<0> > fe0_rt = std::make_shared<FE_RT0_disc<0>>();
    		std::shared_ptr< FiniteElement<1> > fe1_rt = std::make_shared<FE_RT0_disc<1>>();
    		std::shared_ptr< FiniteElement<2> > fe2_rt = std::make_shared<FE_RT0_disc<2>>();
    		std::shared_ptr< FiniteElement<3> > fe3_rt = std::make_shared<FE_RT0_disc<3>>();
    		std::shared_ptr< FiniteElement<0> > fe0_disc = std::make_shared<FE_P_disc<0>>(0);
    		std::shared_ptr< FiniteElement<1> > fe1_disc = std::make_shared<FE_P_disc<1>>(0);
    		std::shared_ptr< FiniteElement<2> > fe2_disc = std::make_shared<FE_P_disc<2>>(0);
    		std::shared_ptr< FiniteElement<3> > fe3_disc = std::make_shared<FE_P_disc<3>>(0);
    		std::shared_ptr< FiniteElement<0> > fe0_cr = std::make_shared<FE_CR<0>>();
    		std::shared_ptr< FiniteElement<1> > fe1_cr = std::make_shared<FE_CR<1>>();
    		std::shared_ptr< FiniteElement<2> > fe2_cr = std::make_shared<FE_CR<2>>();
    		std::shared_ptr< FiniteElement<3> > fe3_cr = std::make_shared<FE_CR<3>>();
    // 	    static FiniteElement<0> fe0_sys = FE_P_disc<0>(0); //TODO fix and use solution with FESystem<0>( {fe0_rt, fe0_disc, fe0_cr} )
    		FESystem<0> fe0_sys( {fe0_disc, fe0_disc, fe0_cr} );
    		FESystem<1> fe1_sys( {fe1_rt, fe1_disc, fe1_cr} );
    		FESystem<2> fe2_sys( {fe2_rt, fe2_disc, fe2_cr} );
    		FESystem<3> fe3_sys( {fe3_rt, fe3_disc, fe3_cr} );
    	    MixedPtr<FESystem> fe_sys( std::make_shared<FESystem<0>>(fe0_sys), std::make_shared<FESystem<1>>(fe1_sys),
    	                                    std::make_shared<FESystem<2>>(fe2_sys), std::make_shared<FESystem<3>>(fe3_sys) );
    		std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_sys);
    		eq_data_->dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    		eq_data_->dh_->distribute_dofs(ds);
        }
    }

    bool zero_time_term(bool time_global=false) {
        if (time_global) {
            return (eq_fields_->storativity.input_list_size() == 0);
        } else {
            return eq_fields_->storativity.field_result(mesh_->region_db().get_region_set("BULK")) == result_zeros;
        }
    }

    void accept_time_step() {
    	eq_data_->p_edge_solution_previous_time.copy_from(eq_data_->p_edge_solution);
    	eq_data_->p_edge_solution_previous_time.local_to_ghost_begin();
    	eq_data_->p_edge_solution_previous_time.local_to_ghost_end();
    }

     /// Initialize equation
    void initialize(const string &input);

    /// Execute zero time step. Do not call method directly, use run_simulation
    void zero_time_step();

    /// Update equation solution. Do not call method directly, use run_simulation
    void update_solution();

    /// Solve method common to zero_time_step and update solution.
    void solve_nonlinear();

    /// Assembly or update whole linear system.
    void assembly_linear_system();

    /// Call zero_time_step and update_solution for 4 time steps.
    void run_simulation() {
        this->zero_time_step();
        for (uint i=0; i<4; ++i) {
            this->update_solution();
        }
    }

    void allocate_mh_matrix()
    {
        // to make space for second schur complement, max. 10 neighbour edges of one el.
        double zeros[100000];
        for(int i=0; i<100000; i++) zeros[i] = 0.0;

        std::vector<LongIdx> tmp_rows;
        tmp_rows.reserve(200);

        std::vector<LongIdx> dofs, dofs_ngh;
        dofs.reserve(eq_data_->dh_cr_->max_elem_dofs());
        dofs_ngh.reserve(eq_data_->dh_cr_->max_elem_dofs());

        // DebugOut() << "Allocate new schur\n";
        for ( DHCellAccessor dh_cell : eq_data_->dh_cr_->own_range() ) {
            ElementAccessor<3> ele = dh_cell.elm();

            const uint ndofs = dh_cell.n_dofs();
            dofs.resize(dh_cell.n_dofs());
            dh_cell.get_dof_indices(dofs);

            int* dofs_ptr = dofs.data();
            eq_data_->lin_sys_schur->mat_set_values(ndofs, dofs_ptr, ndofs, dofs_ptr, zeros);

            tmp_rows.clear();

            // compatible neighborings rows
            unsigned int n_neighs = ele->n_neighs_vb();
            for ( DHCellSide neighb_side : dh_cell.neighb_sides() ) {
                // every compatible connection adds a 2x2 matrix involving
                // current element pressure  and a connected edge pressure

                // read neighbor dofs (dh_cr dofhandler)
                // neighbor cell owning neighb_side
                DHCellAccessor dh_neighb_cell = neighb_side.cell();

                const uint ndofs_ngh = dh_neighb_cell.n_dofs();
                dofs_ngh.resize(ndofs_ngh);
                dh_neighb_cell.get_dof_indices(dofs_ngh);

                // local index of pedge dof on neighboring cell
                tmp_rows.push_back(dofs_ngh[neighb_side.side().side_idx()]);
            }

            eq_data_->lin_sys_schur->mat_set_values(ndofs, dofs_ptr, n_neighs, tmp_rows.data(), zeros); // (edges)  x (neigh edges)
            eq_data_->lin_sys_schur->mat_set_values(n_neighs, tmp_rows.data(), ndofs, dofs_ptr, zeros); // (neigh edges) x (edges)
            eq_data_->lin_sys_schur->mat_set_values(n_neighs, tmp_rows.data(), n_neighs, tmp_rows.data(), zeros);  // (neigh edges) x (neigh edges)

            tmp_rows.clear();

            eq_data_->lin_sys_schur->mat_set_values(ndofs, dofs_ptr, tmp_rows.size(), tmp_rows.data(), zeros);   // master edges x slave edges
            eq_data_->lin_sys_schur->mat_set_values(tmp_rows.size(), tmp_rows.data(), ndofs, dofs_ptr, zeros);   // slave edges  x master edges
            eq_data_->lin_sys_schur->mat_set_values(tmp_rows.size(), tmp_rows.data(), tmp_rows.size(), tmp_rows.data(), zeros);  // slave edges  x slave edges
        }
    }


	int size;				    // global size of MH matrix

	bool data_changed_;

	// Setting of the nonlinear solver. TODO: Move to the solver class later on.
	double tolerance_;
	unsigned int min_n_it_;
	unsigned int max_n_it_;

	std::shared_ptr<equation_data::EqFields> eq_fields_;
	std::shared_ptr<equation_data::EqData> eq_data_;
    Input::Record in_rec_;

    /// general assembly objects, hold assembly objects of appropriate dimension
    GenericAssemblyBase * mh_matrix_assembly_;
//    GenericAssemblyBase * reconstruct_schur_assembly_;
    bool use_linsys_;
};


#endif /* DG_MOCKUP_HH_ */
