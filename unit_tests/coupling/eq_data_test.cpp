/*
 * eq_data_test.cpp
 *
 *  Created on: Dec 18, 2012
 *      Author: jb
 *
 *
 *  This unit test is meant as proof of concept for EqDataBase + Fields.
 *  We demonstrate and test:
 *  - inheritance of fields from some base class that do not perform any action with the data
 *  - declaration of own fields and their registration in EqData constructor
 *  - reading fields from the input
 *  - accessing field data
 */

/*
 * TODO:
 * - Meli bychom mit vsechny rovnice pojmenovane? Jake schema mit pro pojmenovani recordu EqData, ma mit kazda instance konkretni rovnice vlastni
 *   Record (nesmysl)? Tedy v nazvu recordu pouzit nazev tridy
 */

#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include <vector>

#include "system/sys_profiler.hh"

#include "input/input_type.hh"
#include "input/type_output.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"

#include "fields/field_set.hh"
#include "fields/field_add_potential.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
#include "fields/multi_field.hh"
#include "fields/field_constant.hh"
#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"
#include "coupling/equation.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"

#include "mesh/mesh.h"
#include "mesh/bc_mesh.hh"
#include "io/msh_gmshreader.h"
#include "mesh/region.hh"
#include <armadillo>

using namespace std;
namespace IT=Input::Type;





class SomeEquationBase : public EquationBase {
protected:
    class EqData : public FieldSet {
    public:
        enum BC_type {
            none=0,
            dirichlet=1,
            neumann=2,
            robin=3,
            total_flux=4
        };

        static const IT::Selection & get_bc_type_selection();

        EqData()
        {

            *this += anisotropy.name("anisotropy")
                    .description("Anisotropy of the conductivity tensor.")
                    .input_default("1.0")
                    .units( UnitSI::dimensionless() );
            
            *this += bc_type.name("bc_type")
                    .description("Boundary condition type.")
                    .input_selection( get_bc_type_selection() )
                    .input_default("\"none\"")
                    .units( UnitSI::dimensionless() );

            *this += bc_pressure
                    .disable_where(bc_type, {none, neumann} )
                    .name("bc_pressure")
                    .description("Prescribed pressure value on the boundary. Used for all values of ``bc_type`` except ``none`` and ``seepage``. "
                        "See documentation of ``bc_type`` for exact meaning of ``bc_pressure`` in individual boundary condition types.")
                    .input_default("0.0")
                    .units( UnitSI().m() );

            *this += bc_flux
                    .disable_where(bc_type, {none, dirichlet, robin} )
                    .name("bc_flux")
                    .description("Incoming water boundary flux. Used for bc_types : ``total_flux``, ``seepage``, ``river``.")
                    .input_default("0.0")
                    .units( UnitSI().m().s(-1) );

            *this += bc_robin_sigma
                    .disable_where(bc_type, {none, dirichlet, neumann} )
                    .name("bc_robin_sigma")
                    .description("Conductivity coefficient in the ``total_flux`` or the ``river`` boundary condition type.")
                    .input_default("0.0")
                    .units( UnitSI().s(-1) );

            *this += bc_conc.name("bc_conc")
                    .description("Boundary condition for concentration of substances.")
                    .input_default("0.0")
                    .units( UnitSI().kg().m(-3) );
        }

        /// Return coords field
        FieldCoords &X() {
            return this->X_;
        }

        Field<3, FieldValue<3>::TensorFixed > anisotropy;
        BCField<3, FieldValue<3>::Enum > bc_type; // Discrete need Selection for initialization
        BCField<3, FieldValue<3>::Scalar > bc_pressure;
        BCField<3, FieldValue<3>::Scalar > bc_flux;
        BCField<3, FieldValue<3>::Scalar > bc_robin_sigma;
        BCField<3, FieldValue<3>::VectorFixed > bc_conc;

    };

    void output_data() override {}
};

const IT::Selection & SomeEquationBase::EqData::get_bc_type_selection() {
	return IT::Selection("EqData_bc_Type")
             .add_value(none, "none")
             .add_value(dirichlet, "dirichlet")
             .add_value(neumann, "neumann")
             .add_value(robin, "robin")
             .add_value(total_flux, "total_flux")
			 .close();
}



class SomeEquation : public testing::Test, SomeEquationBase {
public:

    class EqData : public SomeEquationBase::EqData, public ElementCacheMap {
    public:

        EqData() : SomeEquationBase::EqData() {
            *this += init_pressure.name("init_pressure").description("Initial condition as pressure.").input_default("0.0");
            *this += init_conc.name("init_conc").description("Initial condition for the concentration.").input_default("0.0");
            *this += bulk_set_field.name("bulk_set_field").input_default("0.0");
            *this += conc_mobile.name("conc_mobile").input_default("0.0");
            *this += bc_gravity.name("bc_gravity").description("Boundary gravity vector.").input_default("0.0");
            *this += bc_piezo_head.name("bc_piezo_head").description("Boundary piezo head.").input_default("0.0");

            init_pressure.units( UnitSI::dimensionless() );
            init_conc.units( UnitSI::dimensionless() );
            bulk_set_field.units( UnitSI::dimensionless() );
            conc_mobile.units( UnitSI::dimensionless() );
            bc_gravity.units( UnitSI::dimensionless() );
            bc_piezo_head.units(UnitSI().m());

            // Asumme following types:
            eval_points_ = std::make_shared<EvalPoints>();
            Quadrature *q_bulk_1d = new QGauss(1, 0);
            Quadrature *q_bulk_2d = new QGauss(2, 0);
            Quadrature *q_bulk_3d = new QGauss(3, 0);
            Quadrature *q_bdr = new QGauss(2, 0);
            mass_eval[0] = eval_points_->add_bulk<1>(*q_bulk_1d );
            mass_eval[1] = eval_points_->add_bulk<2>(*q_bulk_2d );
            mass_eval[2] = eval_points_->add_bulk<3>(*q_bulk_3d );
            bdr_eval = eval_points_->add_boundary<3>(*q_bdr );
            this->init(eval_points_);
        }

        void register_eval_points(bool bdr=false) {
            unsigned int reg_idx = computed_dh_cell_.elm().region_idx().idx();
            unsigned int dim = computed_dh_cell_.dim();
            for (auto p : mass_eval[dim-1]->points(this->position_in_cache(computed_dh_cell_.elm_idx()), this) ) {
                this->eval_point_data_.emplace_back(reg_idx, computed_dh_cell_.elm_idx(), p.eval_point_idx(), computed_dh_cell_.local_idx());
            }

            if (bdr)
                for (DHCellSide cell_side : computed_dh_cell_.side_range()) {
                    if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                        unsigned int bdr_reg = cell_side.cond().element_accessor().region_idx().idx(); // region of boundary element
                        //std::cout << "Bulk elm: " << computed_dh_cell_.elm_idx() << ", boundary element: " << cell_side.cond().bc_ele_idx() << ", region: " << bdr_reg << std::endl;
                        for (auto p : bdr_eval->points(cell_side, this) ) {
                            // point on side of bulk element
                            this->eval_point_data_.emplace_back(reg_idx, cell_side.elem_idx(), p.eval_point_idx(), cell_side.cell().local_idx());
                            // point on boundary element
                            BulkPoint p_bdr = p.point_bdr(cell_side.cond().element_accessor()); // equivalent point on boundary element
                            this->eval_point_data_.emplace_back(bdr_reg, cell_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
                        }
                    }
                }
            this->eval_point_data_.make_permanent();
        }

        void update_cache(bool bdr=false) {
        	this->start_elements_update();
            this->register_eval_points(bdr);
            this->create_patch();
            this->cache_update(*this);
            this->finish_elements_update();
        }


        Field<3, FieldValue<3>::Scalar > init_pressure;
        Field<3, FieldValue<3>::VectorFixed > init_conc;
        Field<3, FieldValue<3>::Scalar > bulk_set_field;
        BCField<3, FieldValue<3>::VectorFixed > bc_gravity;      // Holds gravity vector
        BCField<3, FieldValue<3>::Scalar> bc_piezo_head;
        MultiField<3, FieldValue<3>::Scalar > conc_mobile;
        std::shared_ptr<EvalPoints> eval_points_;
        std::array<std::shared_ptr<BulkIntegral>, 3> mass_eval;  // dim 1,2,3
        std::shared_ptr<BoundaryIntegral> bdr_eval;              // only dim 2
        DHCellAccessor computed_dh_cell_;
    };

protected:
    static Input::Type::Record & get_input_type();
    static MultiField<3, FieldValue<3>::Scalar> empty_mf;
    std::shared_ptr<EqData> data_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    std::vector<string> component_names;

    SomeEquation() {
        Profiler::instance();
        //data.gravity_=arma::vec4("3.0 2.0 1.0 -5.0");
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        data_->add_coords_field();
        data_->set_default_fieldset();
        mesh = mesh_full_constructor("{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh);
        component_names = { "comp_0", "comp_1", "comp_2" };
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        TimeGovernor tg(0.0, 1.0);

        data_->set_components(component_names);        // set number of substances posibly read from elsewhere
        //data_->bc_conc.set_components(component_names);

        /* Regions in the test mesh:
         * $PhysicalNames
            6
            1       37      "1D diagonal"
            2       38      "2D XY diagonal"
            2       101     ".top side"
            2       102     ".bottom side"
            3       39      "3D back"
            3       40      "3D front"
            $EndPhysicalNames
         */

        DebugOut() << "init\n";

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data_->set_mesh(*mesh);

        arma::vec3 gravity_vec("3 2 1");
        FieldValue<3>::VectorFixed gvalue(gravity_vec);
        auto field_algo=std::make_shared<FieldConstant<3, FieldValue<3>::VectorFixed>>();
        field_algo->set_value(gvalue);
        data_->bc_gravity.set(field_algo, 0.0);
        data_->bc_pressure.add_factory(
                std::make_shared<AddPotentialFactory<3, FieldValue<3>::Scalar> >
                (data_->bc_gravity, data_->X(), data_->bc_piezo_head) );

        data_->set_input_list( inputs[input_last], tg );
        data_->set_time(tg.step(), LimitSide::right);
        data_->cache_reallocate( *(data_.get()), *(data_.get()) );
    }

    ~SomeEquation() {
        delete mesh;
    };

    void set_dh_cell(unsigned int elm_idx, unsigned int reg_id, bool updt_bdr=false) {
        ElementAccessor<3> elm = mesh->element_accessor(elm_idx);
        EXPECT_EQ(reg_id, elm.region().id()); // check region id
        data_->computed_dh_cell_ = this->dh_->cell_accessor_from_element(elm.idx());
        data_->update_cache(updt_bdr);
    }


    Mesh * mesh;
};



MultiField<3, FieldValue<3>::Scalar> SomeEquation::empty_mf = MultiField<3, FieldValue<3>::Scalar>();

IT::Record & SomeEquation::get_input_type() {
	return IT::Record("SomeEquation","")
	        .declare_key("data", IT::Array(
	        		IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
	                .copy_keys( SomeEquation::EqData().make_field_descriptor_type("SomeEquation") )
	                .declare_key("bc_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
	                .declare_key("init_piezo_head", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
					.close()
	                ), IT::Default::obligatory(), ""  )
			.close();
}



TEST_F(SomeEquation, values) {
    // Test input for 'values' test
    string eq_data_input = R"YAML(
    data:
      - rid: 37
        init_pressure: !FieldConstant
          value: 1.1
        init_conc: [ 1, 2, 3 ]
        # MultiField
        conc_mobile: !FieldConstant 
          value: [1, 2, 3]
      - region: ["2D XY diagonal", "3D back"]
        init_pressure: 2.2
        conc_mobile: [1, 2, 3]
      - region: "BULK"
        bulk_set_field: 5.7
        conc_mobile: [1, 2, 3]

      # boundary          
      - rid: 101
        bc_type: !FieldConstant
          value: "dirichlet"
        bc_pressure: !FieldConstant
          value: 1.23
        bc_conc: !FieldFormula
          value: "[X[0], 10+X[0], 20+X[0]]"
        conc_mobile: !FieldConstant 
          value: [5, 6, 7]
      - rid: 102
        bc_type: "dirichlet"
        bc_piezo_head: 1.23
        bc_conc: !FieldFormula
          value: "[X[0], 10+X[0], 20+X[0]]"
        conc_mobile: [5, 6, 7]
    )YAML";

    read_input(eq_data_input);

    DebugOut().fmt("elements size: {} {}\n", mesh->n_elements(), mesh->bc_mesh()->n_elements());

    {
        // region 37 "1D diagonal"
        this->set_dh_cell(0, 37);

        auto p = *( data_->mass_eval[0]->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()).begin() );

        EXPECT_DOUBLE_EQ(1.1, data_->init_pressure(p) );            // scalar
        for (unsigned int i=0; i<data_->conc_mobile.size(); ++i) {  // multifield
            auto conc_mobile_val = data_->conc_mobile[i](p);
            EXPECT_DOUBLE_EQ( 1.0 + i, conc_mobile_val );
        }
        auto anisotropy_value = data_->anisotropy(p);               // tensor
        arma::mat33 expected_aniso = {1.0, 0, 0, 0, 1, 0, 0, 0, 1};
        EXPECT_ARMA_EQ(expected_aniso, anisotropy_value);
        auto init_conc_value = data_->init_conc(p);                 // vector init_conc
        arma::vec3 expected_conc = {1.0, 2.0, 3.0};
        EXPECT_ARMA_EQ(expected_conc, init_conc_value);
        EXPECT_EQ( 5.7, data_->bulk_set_field(p) );                 // bulk_set_filed - test setting on region set, test setting field without default value
    }

    {
        // region 38 "2D XY diagonal"
        this->set_dh_cell(1, 38);

        auto p = *( data_->mass_eval[0]->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()).begin() );

        EXPECT_DOUBLE_EQ(2.2, data_->init_pressure(p) );            // scalar
        for (unsigned int i=0; i<data_->conc_mobile.size(); ++i) {  // multifield
            auto conc_mobile_val = data_->conc_mobile[i](p);
            EXPECT_DOUBLE_EQ( 1.0 + i, conc_mobile_val );
        }
        EXPECT_EQ( 5.7, data_->bulk_set_field(p) );                 // bulk_set_filed
    }

    {
        // region 39 "3D back"
        this->set_dh_cell(3, 39);

        auto p = *( data_->mass_eval[0]->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()).begin() );

        EXPECT_DOUBLE_EQ(2.2, data_->init_pressure(p) );
        EXPECT_EQ( 5.7, data_->bulk_set_field(p) );
    }

    {
    	// bulk element of boundary element in region 101 ".top side"
    	this->set_dh_cell(8, 40, true);

        bool found_bdr = false;
        for (DHCellSide cell_side : data_->computed_dh_cell_.side_range()) {
            if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                if (cell_side.cond().bc_ele_idx() == 1) {
                    found_bdr = true;
                    EXPECT_EQ(101, cell_side.cond().element_accessor().region().id());

                    auto p_side = *( data_->bdr_eval->points(cell_side, data_.get()).begin() );
                    auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
                    EXPECT_EQ( EqData::dirichlet, data_->bc_type(p_bdr) );
                    EXPECT_DOUBLE_EQ(1.23, data_->bc_pressure(p_bdr) );
                }
            }
        }
        EXPECT_TRUE( found_bdr );
    }

    {
        // bulk element of boundary element in region 102 ".bottom side"
        this->set_dh_cell(6, 40, true);

        bool found_bdr = false;
        for (DHCellSide cell_side : data_->computed_dh_cell_.side_range()) {
            if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                if (cell_side.cond().bc_ele_idx() == 3) {
                    found_bdr = true;
                    EXPECT_EQ(102, cell_side.cond().element_accessor().region().id());

                    auto p_side = *( data_->bdr_eval->points(cell_side, data_.get()).begin() );
                    auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
                    arma::vec3 elm_center = cell_side.cond().element_accessor().centre();
                    EXPECT_EQ( EqData::dirichlet, data_->bc_type(p_bdr) );
                    EXPECT_DOUBLE_EQ(1.23 + (3*elm_center(0) + 2*elm_center(1) + elm_center(2)), data_->bc_pressure(p_bdr) );    // piezo_head

                    arma::vec bc_value = data_->bc_conc(p_bdr);
                    EXPECT_DOUBLE_EQ(elm_center(0), bc_value(0) );
                    EXPECT_DOUBLE_EQ(elm_center(0)+10, bc_value(1) );
                    EXPECT_DOUBLE_EQ(elm_center(0)+20, bc_value(2) );
                }
            }
        }
        EXPECT_TRUE( found_bdr );
    }

    std::stringstream ss;
    if ( FieldCommon::print_message_table(ss, "test") ) {
        WarningOut() << ss.str();
    }
}



TEST_F(SomeEquation, wrong_time_order) {
    // Test input for wrong_time_order
    string eq_data = R"YAML(
    data:
      - region: "BULK"
        time: 1.0
        bulk_set_field: 0.0
      - region: "BULK"
        time: 0.0
        bulk_set_field: 1.0
    )YAML";

    EXPECT_THROW({read_input(eq_data);}, FieldCommon::ExcNonascendingTime);
}
