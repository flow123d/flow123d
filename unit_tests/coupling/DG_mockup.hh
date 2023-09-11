#ifndef DG_MOCKUP_HH_
#define DG_MOCKUP_HH_


#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
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

template<unsigned int dim> class MassAssembly;
template<unsigned int dim> class StiffnessAssembly;
template<unsigned int dim> class SourcesAssembly;
template< template<IntDim...> class DimAssembly> class GenericAssembly;


/*******************************************************************************
 * Functors of FieldModels
 */
using Sclr = double;
using Vect = arma::vec3;
using Tens = arma::mat33;

// Functor computing velocity norm
struct fn_conc_v_norm {
    inline Sclr operator() (Vect vel) {
        return arma::norm(vel, 2);
    }
};

// Functor computing mass matrix coefficients (cross_section * water_content)
struct fn_conc_mass_matrix {
	inline Sclr operator() (Sclr csec, Sclr wcont) {
        return csec * wcont;
    }
};

// Functor computing retardation coefficients:
// (1-porosity) * rock_density * sorption_coefficient * rock_density
struct fn_conc_retardation {
    inline Sclr operator() (Sclr csec, Sclr por_m, Sclr rho_s, Sclr sorp_mult) {
        return (1.-por_m)*rho_s*sorp_mult*csec;
    }
};

// Functor computing sources density output (cross_section * sources_density)
struct fn_conc_sources_dens {
    inline Sclr operator() (Sclr csec, Sclr sdens) {
        return csec * sdens;
    }
};

// Functor computing sources sigma output (cross_section * sources_sigma)
struct fn_conc_sources_sigma {
    inline Sclr operator() (Sclr csec, Sclr ssigma) {
        return csec * ssigma;
    }
};

// Functor computing sources concentration output (sources_conc)
struct fn_conc_sources_conc {
    inline Sclr operator() (Sclr sconc) {
        return sconc;
    }
};

// Functor computing advection coefficient (velocity)
struct fn_conc_ad_coef {
    inline Vect operator() (Vect velocity) {
        return velocity;
    }
};

// Functor computing diffusion coefficient (see notes in function)
struct fn_conc_diff_coef {
    inline Tens operator() (Tens diff_m, Vect velocity, Sclr v_norm, Sclr alphaL, Sclr alphaT, Sclr water_content, Sclr porosity, Sclr c_sec) {

        // used tortuosity model dues to Millington and Quirk(1961) (should it be with power 10/3 ?)
        // for an overview of other models see: Chou, Wu, Zeng, Chang (2011)
        double tortuosity = pow(water_content, 7.0 / 3.0)/ (porosity * porosity);

        // result
        Tens K;

        // Note that the velocity vector is in fact the Darcian flux,
        // so we need not to multiply vnorm by water_content and cross_section.
	    //K = ((alphaL-alphaT) / vnorm) * K + (alphaT*vnorm + Dm*tortuosity*cross_cut*water_content) * arma::eye(3,3);

        if (fabs(v_norm) > 0) {
            /*
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    K(i,j) = (velocity[i]/vnorm)*(velocity[j]);
            */
            K = ((alphaL - alphaT) / v_norm) * arma::kron(velocity.t(), velocity);

            //arma::mat33 abs_diff_mat = arma::abs(K -  kk);
            //double diff = arma::min( arma::min(abs_diff_mat) );
            //ASSERT_PERMANENT(  diff < 1e-12 )(diff)(K)(kk);
        } else
            K.zeros();

        // Note that the velocity vector is in fact the Darcian flux,
        // so to obtain |v| we have to divide vnorm by porosity and cross_section.
        K += alphaT*v_norm*arma::eye(3,3) + diff_m*(tortuosity*c_sec*water_content);

        return K;
    }
};


class DGMockupTest : public testing::Test {
public:
	DGMockupTest()
    {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        Profiler::set_memory_monitoring(false, false);
    }

    ~DGMockupTest() {}

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		static ofstream os( FilePath("benchmark_" + file_name + "_test.log", FilePath::output_file) );
		Profiler::instance()->output(MPI_COMM_WORLD, os);
	}
};

/*******************************************************************************
 * Test class
 */
class DGMockup : public EquationBase {
public:
    typedef std::vector<std::shared_ptr<FieldFE< 3, FieldValue<3>::Scalar>>> FieldFEScalarVec;

    enum Abstract_bc_types {
//        abc_none,
        abc_inflow,
        abc_dirichlet,
        abc_total_flux,
        abc_diffusive_flux
    };

    static Input::Type::Record & get_input_type() {
        std::string equation_name = "TestEquation";
        return Input::Type::Record(equation_name, "Benchmark test equation record.")
            .declare_key("solver", LinSys_PETSC::get_input_type(), Input::Type::Default("{}"),
                    "Solver for the linear system.")
            .declare_key("input_fields", Input::Type::Array(
                    DGMockup::EqFields()
                        .make_field_descriptor_type(equation_name)),
                    IT::Default::obligatory(),
                    "Input fields of the equation.")
            //.declare_key("output",
            //        EqFields().output_fields.make_output_type(equation_name, ""),
            //        Input::Type::Default("{ \"fields\": [ " + Model::ModelEqData::default_output_field() + "] }"),
            //        "Specification of output fields and output times.")
            .close();
    }

    class EqFields : public FieldSet {
    public:
        enum Concentration_bc_types {
            bc_inflow,
            bc_dirichlet,
            bc_total_flux,
            bc_diffusive_flux
        };

        static const Input::Type::Selection & get_bc_type_selection() {
        	return Input::Type::Selection("Solute_AdvectionDiffusion_BC_Type", "Types of boundary conditions for advection-diffusion solute transport model.")
                      .add_value(bc_inflow, "inflow",
                              "Default transport boundary condition.")
                      .add_value(bc_dirichlet, "dirichlet",
                              "Dirichlet boundary condition (($ c = c_D $)).")
                      .add_value(bc_total_flux, "total_flux",
                              "Total mass flux boundary condition.")
                      .add_value(bc_diffusive_flux, "diffusive_flux",
                              "Diffusive flux boundary condition.")
        			  .close();
        }

        /// Constructor
        EqFields()
        {
            *this += porosity.name("porosity")
                    .description("Porosity of the mobile phase.")
                    .input_default("1.0")
                    .units( UnitSI::dimensionless() )
                    .flags_add(in_main_matrix & in_rhs);
            *this += water_content.name("water_content")
                    .description("INTERNAL. Water content passed from unsaturated Darcy flow model.")
					.input_default("0.0")
                    .units( UnitSI::dimensionless() )
                    .flags_add(input_copy & in_time_term & in_main_matrix & in_rhs);
            *this += cross_section.name("cross_section")
                    .description("Complement dimension parameter (cross section for 1D, thickness for 2D).")
                    .input_default("1.0")
                    .units( UnitSI().m(3).md() )
                    .flags_add(in_time_term & in_main_matrix & in_rhs);
            *this += flow_flux.name("flow_flux")
	    	         .units(UnitSI().m().s(-1))
                     .input_default("1.0")
	                 .flags(FieldFlag::equation_result)
                     .flags_add(in_time_term & in_main_matrix & in_rhs);
            *this += sources_density.name("sources_density")
                    .description("Density of concentration sources.")
                    .input_default("0.0")
                    .units( UnitSI().kg().m(-3).s(-1) )
                    .flags_add(in_rhs);
            *this += sources_sigma.name("sources_sigma")
                    .description("Concentration flux.")
                    .input_default("0.0")
                    .units( UnitSI().s(-1) )
                    .flags_add(in_main_matrix & in_rhs);
            *this += sources_conc.name("sources_conc")
                    .description("Concentration sources threshold.")
                    .input_default("0.0")
                    .units( UnitSI().kg().m(-3) )
                    .flags_add(in_rhs);

            *this+=bc_type
                    .name("bc_type")
                    .description(
                    "Type of boundary condition.")
                    .units( UnitSI::dimensionless() )
                    .input_default("\"inflow\"")
                    .input_selection( get_bc_type_selection() )
                    .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
            *this+=bc_dirichlet_value
                    .name("bc_conc")
                    .units( UnitSI().kg().m(-3) )
                    .description("Dirichlet boundary condition (for each substance).")
                    .input_default("0.0")
                    .flags_add( in_rhs );
        	*this+=bc_flux
        	        .disable_where(bc_type, { bc_inflow, bc_dirichlet })
        			.name("bc_flux")
        			.description("Flux in Neumann boundary condition.")
        			.units( UnitSI().kg().m().s(-1).md() )
        			.input_default("0.0")
        			.flags_add(FieldFlag::in_rhs);
        	*this+=bc_robin_sigma
        	        .disable_where(bc_type, { bc_inflow, bc_dirichlet })
        			.name("bc_robin_sigma")
        			.description("Conductivity coefficient in Robin boundary condition.")
        			.units( UnitSI().m(4).s(-1).md() )
        			.input_default("0.0")
        			.flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
            *this+=init_condition
                    .name("init_conc")
                    .units( UnitSI().kg().m(-3) )
                    .description("Initial values for concentration of substances.")
                    .input_default("0.0");
            *this+=disp_l
                    .name("disp_l")
                    .description("Longitudinal dispersivity in the liquid (for each substance).")
                    .units( UnitSI().m() )
                    .input_default("0.0")
                    .flags_add( in_main_matrix & in_rhs );
            *this+=disp_t
                    .name("disp_t")
                    .description("Transverse dispersivity in the liquid (for each substance).")
                    .units( UnitSI().m() )
                    .input_default("0.0")
                    .flags_add( in_main_matrix & in_rhs );
            *this+=diff_m
                    .name("diff_m")
                    .description("Molecular diffusivity in the liquid (for each substance).")
                    .units( UnitSI().m(2).s(-1) )
                    .input_default("0.0")
                    .flags_add( in_main_matrix & in_rhs );
            *this+=rock_density
            		.name("rock_density")
        			.description("Rock matrix density.")
        			.units(UnitSI().kg().m(-3))
        			.input_default("0.0")
        			.flags_add( in_time_term );
            *this+=sorption_coefficient
            		.name("sorption_coefficient")
        			.description("Coefficient of linear sorption.")
        			.units(UnitSI().m(3).kg(-1))
        			.input_default("0.0")
        			.flags_add( in_time_term );
        	*this+=output_field
        	        .name("conc")
                    .description("Concentration solution.")
        	        .units( UnitSI().kg().m(-3) )
        	        .flags( equation_result );
        	// initiaization of FieldModels
            *this += v_norm.name("v_norm")
                    .description("Velocity norm field.")
                    .input_default("0.0")
                    .units( UnitSI().m().s(-1) );
            *this += mass_matrix_coef.name("mass_matrix_coef")
                    .description("Matrix coefficients computed by model in mass assemblation.")
                    .input_default("0.0")
                    .units( UnitSI().m(3).md() );
            *this += retardation_coef.name("retardation_coef")
                    .description("Retardation coefficients computed by model in mass assemblation.")
                    .input_default("0.0")
                    .units( UnitSI().m(3).md() );
            *this += sources_density_out.name("sources_density_out")
                    .description("Concentration sources output - density of substance source, only positive part is used..")
                    .input_default("0.0")
                    .units( UnitSI().kg().s(-1).md() );
            *this += sources_sigma_out.name("sources_sigma_out")
                    .description("Concentration sources - Robin type, in_flux = sources_sigma * (sources_conc - mobile_conc).")
                    .input_default("0.0")
                    .units( UnitSI().s(-1).m(3).md() );
            *this += sources_conc_out.name("sources_conc_out")
                    .description("Concentration sources output.")
                    .input_default("0.0")
                    .units( UnitSI().kg().m(-3) );
            *this += advection_coef.name("advection_coef")
                    .description("Advection coefficients model.")
                    .input_default("0.0")
                    .units( UnitSI().m().s(-1) );
            *this += diffusion_coef.name("diffusion_coef")
                    .description("Diffusion coefficients model.")
                    .input_default("0.0")
                    .units( UnitSI().m(2).s(-1) );

            *this+=fracture_sigma
                    .name("fracture_sigma")
                    .description(
                    "Coefficient of diffusive transfer through fractures (for each substance).")
                    .units( UnitSI::dimensionless() )
                    .input_default("1.0")
                    .flags_add(FieldFlag::in_main_matrix);
            *this+=dg_penalty
                    .name("dg_penalty")
                    .description(
                    "Penalty parameter influencing the discontinuity of the solution (for each substance). "
                    "Its default value 1 is sufficient in most cases. Higher value diminishes the inter-element jumps.")
                    .units( UnitSI::dimensionless() )
                    .input_default("1.0")
                    .flags_add(FieldFlag::in_rhs & FieldFlag::in_main_matrix);
            *this += region_id.name("region_id")
                        .units( UnitSI::dimensionless())
                        .flags(FieldFlag::equation_external_output)
                        .description("Region ids.");
            *this += subdomain.name("subdomain")
              .units( UnitSI::dimensionless() )
              .flags(FieldFlag::equation_external_output)
              .description("Subdomain ids of the domain decomposition.");

            // add all input fields to the output list
            output_fields += *this;

            this->add_coords_field();
            this->set_default_fieldset();
        }

        /// Setup components of MultiFields defined as FieldModels, internal method
        void setup_mf_components()
        {
            // initialize multifield components
        	sorption_coefficient.setup_components();
            sources_conc.setup_components();
            sources_density.setup_components();
            sources_sigma.setup_components();
            diff_m.setup_components();
            disp_l.setup_components();
            disp_t.setup_components();

        }

        /// Initialize selected fields as FieldModels
        void init_field_models()
        {
        	setup_mf_components();

            v_norm.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_v_norm(), flow_flux), 0.0);
            mass_matrix_coef.set(Model<3, FieldValue<3>::Scalar>::create(fn_conc_mass_matrix(), cross_section, water_content), 0.0);
            retardation_coef.set(
                Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_retardation(), cross_section, porosity, rock_density, sorption_coefficient),
        		0.0
            );
            sources_density_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_dens(), cross_section, sources_density), 0.0);
            sources_sigma_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_sigma(), cross_section, sources_sigma), 0.0);
            sources_conc_out.set(Model<3, FieldValue<3>::Scalar>::create_multi(fn_conc_sources_conc(), sources_conc), 0.0);
            std::vector<typename Field<3, FieldValue<3>::VectorFixed>::FieldBasePtr> ad_coef_ptr_vec;
            for (unsigned int sbi=0; sbi<sorption_coefficient.size(); sbi++)
                ad_coef_ptr_vec.push_back( Model<3, FieldValue<3>::VectorFixed>::create(fn_conc_ad_coef(), flow_flux) );
            advection_coef.set(ad_coef_ptr_vec, 0.0);
            diffusion_coef.set(
                Model<3, FieldValue<3>::TensorFixed>::create_multi(
                    fn_conc_diff_coef(), diff_m, flow_flux, v_norm, disp_l, disp_t, water_content, porosity, cross_section
                ),
                0.0
            );
        }

        /// Initialize selected fields as FieldConstants
        void init_field_constants(double vnorm, double mass_matrix, double ret_coef, double sources_dens, double sourc_sigma,
                double sourc_conc, arma::vec3 advec_coef, arma::mat33 diff_coef)
        {
        	setup_mf_components();

            {
        	    FieldValue<3>::Scalar f_value(vnorm);
                auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
                field_algo->set_value(f_value);
                v_norm.set(field_algo, 0.0);
            }
            {
        	    FieldValue<3>::Scalar f_value(mass_matrix);
                auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
                field_algo->set_value(f_value);
                mass_matrix_coef.set(field_algo, 0.0);
            }
            {
                FieldValue<3>::Scalar f_value(ret_coef);
                std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
                for (unsigned int sbi=0; sbi<sorption_coefficient.size(); sbi++) {
                    auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
                    field_algo->set_value(f_value);
                    field_vec.push_back(field_algo);
                }
                retardation_coef.set(field_vec, 0.0);
            }
            {
                FieldValue<3>::Scalar f_value(sources_dens);
                std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
                for (unsigned int sbi=0; sbi<sources_density.size(); sbi++) {
                    auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
                    field_algo->set_value(f_value);
                    field_vec.push_back(field_algo);
                }
                sources_density_out.set(field_vec, 0.0);
            }
            {
                FieldValue<3>::Scalar f_value(sourc_sigma);
                std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
                for (unsigned int sbi=0; sbi<sources_sigma.size(); sbi++) {
                    auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
                    field_algo->set_value(f_value);
                    field_vec.push_back(field_algo);
                }
                sources_sigma_out.set(field_vec, 0.0);
            }
            {
                FieldValue<3>::Scalar f_value(sourc_conc);
                std::vector<typename Field<3, FieldValue<3>::Scalar>::FieldBasePtr> field_vec;
                for (unsigned int sbi=0; sbi<sources_conc.size(); sbi++) {
                    auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::Scalar>>();
                    field_algo->set_value(f_value);
                    field_vec.push_back(field_algo);
                }
                sources_conc_out.set(field_vec, 0.0);
            }
            {
                FieldValue<3>::TensorFixed f_value(diff_coef);
                std::vector<typename Field<3, FieldValue<3>::TensorFixed>::FieldBasePtr> field_vec;
                for (unsigned int sbi=0; sbi<diff_m.size(); sbi++) {
                    auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::TensorFixed>>();
                    field_algo->set_value(f_value);
                    field_vec.push_back(field_algo);
                }
                diffusion_coef.set(field_vec, 0.0);
            }
        }

        // from TransportEqFields
        Field<3, FieldValue<3>::Scalar> porosity;             ///< Mobile porosity - usually saturated water content in the case of unsaturated flow model
        Field<3, FieldValue<3>::Scalar> water_content;        ///v Water content - result of unsaturated water flow model or porosity
        Field<3, FieldValue<3>::Scalar > cross_section;       ///< Pointer to DarcyFlow field cross_section
        Field<3, FieldValue<3>::VectorFixed > flow_flux;      ///< Flow flux, can be result of water flow model.
        MultiField<3, FieldValue<3>::Scalar> sources_density; ///< Concentration sources - density of substance source, only positive part is used.
        MultiField<3, FieldValue<3>::Scalar> sources_sigma;   ///< Concentration sources - Robin type, in_flux = sources_sigma * (sources_conc - mobile_conc)
        MultiField<3, FieldValue<3>::Scalar> sources_conc;

        // from ConcentrationModel
        BCMultiField<3, FieldValue<3>::Enum > bc_type;               ///< Type of boundary condition (see also BC_Type)
	    BCMultiField<3, FieldValue<3>::Scalar> bc_dirichlet_value;   ///< Prescribed concentration for Dirichlet/reference concentration for flux b.c.
	    BCMultiField<3, FieldValue<3>::Scalar > bc_flux;             ///< Flux value in total/diffusive flux b.c.
	    BCMultiField<3, FieldValue<3>::Scalar > bc_robin_sigma;      ///< Transition coefficient in total/diffusive flux b.c.
	    MultiField<3, FieldValue<3>::Scalar> init_condition;         ///< Initial concentrations.
	    MultiField<3, FieldValue<3>::Scalar> disp_l;                 ///< Longitudal dispersivity (for each substance).
	    MultiField<3, FieldValue<3>::Scalar> disp_t;                 ///< Transversal dispersivity (for each substance).
	    MultiField<3, FieldValue<3>::TensorFixed> diff_m;            ///< Molecular diffusivity (for each substance).
	    Field<3, FieldValue<3>::Scalar > rock_density;               ///< Rock matrix density.
	    MultiField<3, FieldValue<3>::Scalar > sorption_coefficient;  ///< Coefficient of linear sorption.
	    MultiField<3, FieldValue<3>::Scalar> output_field;
	    // from ConcentrationModel - instances of FieldModel
	    Field<3, FieldValue<3>::Scalar > mass_matrix_coef;           ///< Field represents coefficients of mass matrix.
        MultiField<3, FieldValue<3>::Scalar> retardation_coef;       ///< Field represents retardation coefficients due to sorption.
        MultiField<3, FieldValue<3>::Scalar> sources_density_out;    ///< Concentration sources - density output
        MultiField<3, FieldValue<3>::Scalar> sources_sigma_out;      ///< Concentration sources - sigma output
        MultiField<3, FieldValue<3>::Scalar> sources_conc_out;       ///< Concentration sources - concentration output
	    MultiField<3, FieldValue<3>::VectorFixed> advection_coef;    ///< Advection coefficients.
	    MultiField<3, FieldValue<3>::TensorFixed> diffusion_coef;    ///< Diffusion coefficients.
        Field<3, FieldValue<3>::Scalar > v_norm;                     ///< Velocity norm field.

        // from TransportDG
        MultiField<3, FieldValue<3>::Scalar> fracture_sigma;    ///< Transition parameter for diffusive transfer on fractures (for each substance).
	    MultiField<3, FieldValue<3>::Scalar> dg_penalty;        ///< Penalty enforcing inter-element continuity of solution (for each substance).
        Field<3, FieldValue<3>::Scalar> region_id;
        Field<3, FieldValue<3>::Scalar> subdomain;

        EquationOutput output_fields;
    };

    class EqData {
	public:
        EqData() {}


        int dg_variant;                           ///< DG variant ((non-)symmetric/incomplete
        unsigned int dg_order;                    ///< Polynomial order of finite elements.
        std::shared_ptr<DOFHandlerMultiDim> dh_;  ///< Object for distribution of dofs.
        std::vector<std::string> substances_;     ///< Names of substances

        vector<vector<arma::mat33> > dif_coef;    ///< Diffusion coefficients.
        unsigned int max_edg_sides;               ///< Maximal number of edge sides (evaluate from dim 1,2,3)

        std::vector<Vec> ret_vec;  ///< Auxiliary vectors for calculation of sources in balance due to retardation (e.g. sorption).
        LinSys **ls;               ///< Linear algebra system for the transport equation.
        LinSys **ls_dt;            ///< Linear algebra system for the time derivative (actually it is used only for handling the matrix structures).

        std::vector<VectorMPI> output_vec;        ///< Vector of solution data.

        FieldFEScalarVec conc_fe;
        std::shared_ptr<DOFHandlerMultiDim> dh_p0;

        inline unsigned int n_substances() const {
            return substances_.size();
        }

        inline const std::vector<std::string> &subst_names() const {
            return substances_;
        }

        double elem_anisotropy(ElementAccessor<3> e) const
        {
            double h_max = 0, h_min = numeric_limits<double>::infinity();
            for (unsigned int i=0; i<e->n_nodes(); i++)
                for (unsigned int j=i+1; j<e->n_nodes(); j++)
                {
                    double dist = arma::norm(*e.node(i) - *e.node(j));
                    h_max = max(h_max, dist);
                    h_min = min(h_min, dist);
                }
            return h_max/h_min;
        }

        void set_DG_parameters_boundary(Side side,
                    const int K_size,
                    const vector<arma::mat33> &K,
                    const double flux,
                    const arma::vec3 &normal_vector,
                    const double alpha,
                    double &gamma)
        {
            double delta = 0, h = 0;

            // calculate the side diameter
            if (side.dim() == 0)
            {
                h = 1;
            }
            else
            {
                for (unsigned int i=0; i<side.n_nodes(); i++)
                    for (unsigned int j=i+1; j<side.n_nodes(); j++) {
                        double dist = arma::norm(*side.node(i) - *side.node(j));
                        h = max(h, dist);
                    }

            }

            // delta is set to the average value of Kn.n on the side
            for (int k=0; k<K_size; k++)
                delta += dot(K[k]*normal_vector,normal_vector);
            delta /= K_size;

            gamma = 0.5*fabs(flux) + alpha/h*delta*elem_anisotropy(side.element());
        }
    };

	enum DGVariant {
		non_symmetric = -1,  // Non-symmetric weighted interior penalty DG
		incomplete = 0,      // Incomplete weighted interior penalty DG
		symmetric = 1        // Symmetric weighted interior penalty DG
	};

	DGMockup()
    {
        eq_data_ = make_shared<EqData>();
        eq_fields_ = make_shared<EqFields>();
        this->eq_fieldset_ = eq_fields_;

        // DG data parameters
        eq_data_->dg_variant = DGVariant::non_symmetric;
        eq_data_->dg_order = 1;
    }

    ~DGMockup() {
        if (eq_data_->dif_coef.size() > 0) {
            // initialize called

            for (unsigned int i=0; i<eq_data_->n_substances(); i++)
            {
                if (eq_data_->ls != nullptr) {
                    delete eq_data_->ls[i];
                    delete eq_data_->ls_dt[i];
                }

                if (stiffness_matrix.size() > 0) {
                    if (stiffness_matrix[i])
                        chkerr(MatDestroy(&stiffness_matrix[i]));
                    if (mass_matrix[i])
                        chkerr(MatDestroy(&mass_matrix[i]));
                    if (rhs[i])
                        chkerr(VecDestroy(&rhs[i]));
                    if (mass_vec[i])
                        chkerr(VecDestroy(&mass_vec[i]));
                    if (eq_data_->ret_vec[i])
                        chkerr(VecDestroy(&eq_data_->ret_vec[i]));
                }
            }
            if (eq_data_->ls != nullptr) {
                delete[] eq_data_->ls;
                delete[] eq_data_->ls_dt;
                eq_data_->ls = nullptr;
            }

            if (mass_assembly_ != nullptr) {
                delete mass_assembly_;
                delete stiffness_assembly_;
                delete sources_assembly_;
            }
        }

    }

    void create_and_set_mesh(const std::string &mesh_file) {
        std::string input_str = "{ mesh_file=\"" + mesh_file + "\" }";
        this->mesh_ = mesh_full_constructor(input_str);

        // Set up physical parameters.
        eq_fields_->set_mesh(*mesh_);
        eq_fields_->region_id = GenericField<3>::region_id(*this->mesh_);
        eq_fields_->subdomain = GenericField<3>::subdomain(*this->mesh_);

        MixedPtr<FE_P_disc> fe(eq_data_->dg_order);
        shared_ptr<DiscreteSpace> ds = make_shared<EqualOrderDiscreteSpace>(this->mesh_, fe);
        eq_data_->dh_ = make_shared<DOFHandlerMultiDim>(*this->mesh_);
        eq_data_->dh_->distribute_dofs(ds);
    }

    /// Initialize equation
    void initialize(const string &input, std::vector<std::string> substances);

    /// Execute zero time step. Do not call method directly, use run_simulation
    void zero_time_step();

    /// Update equation solution. Do not call method directly, use run_simulation
    void update_solution();

    /// Call zero_time_step and update_solution for 4 time steps.
    void run_simulation() {
        this->zero_time_step();
        for (uint i=0; i<4; ++i) {
            this->update_solution();
        }
    }


    std::shared_ptr<EqFields> eq_fields_;  ///< Fields for model parameters.
    std::shared_ptr<EqData> eq_data_;      ///< Data for model parameters.
    Input::Record in_rec_;

	std::vector<Vec> rhs;               ///< Vector of right hand side.
	std::vector<Mat> stiffness_matrix;  ///< The stiffness matrix.
	std::vector<Mat> mass_matrix;       ///< The mass matrix.
	std::vector<Vec> mass_vec;          ///< Mass from previous time instant (necessary when coefficients of mass matrix change in time).

	GenericAssembly< MassAssembly > * mass_assembly_;
    GenericAssembly< StiffnessAssembly > * stiffness_assembly_;
    GenericAssembly< SourcesAssembly > * sources_assembly_;
};


#endif /* DG_MOCKUP_HH_ */
