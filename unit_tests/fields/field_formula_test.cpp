/*
 * field_formula_test.cpp
 *
 *  Created on: Jan 8, 2013
 *      Author: jb
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_formula.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"



class FieldEvalFormulaTest : public testing::Test {

public:
    class EqData : public FieldSet, public ElementCacheMap {
    public:
        EqData() {
            *this += vector_field
                        .name("vector_field")
                        .description("Velocity vector.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().m().s(-1) );
            *this += density_unit_conversion
                        .name("density_unit_conversion")
                        .description("Density with unit conversion.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg().m(-3) );
            *this += scalar_field
                        .name("scalar_field")
                        .description("Pressure head")
                        .units( UnitSI().m() );
            *this += scalar_z
                        .name("scalar_z")
                        .description("Pressure head")
                        .input_default("0.0")
                        .units( UnitSI().m() );
            *this += scalar_with_depth
                        .name("scalar_with_depth")
                        .description("FieldFormula with depth")
                        .input_default("0.0")
                        .units( UnitSI().m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += const_scalar
                        .name("const_scalar")
                        .input_default("0.0")
                        .description("")
                        .units( UnitSI::dimensionless() );
            *this += integer_scalar
                        .name("integer_scalar")
                        .input_default("0")
                        .description("")
                        .units( UnitSI::dimensionless() );

            // Asumme following types:
            eval_points_ = std::make_shared<EvalPoints>();
            Quadrature *q_bulk = new QGauss(3, 2);
            Quadrature *q_side = new QGauss(2, 2);
            mass_eval = eval_points_->add_bulk<3>(*q_bulk );
            side_eval = eval_points_->add_edge<3>(*q_side );
            // ngh_side_eval = ...
            this->init(eval_points_);

            this->add_coords_field();
            this->set_default_fieldset();
        }

        void register_eval_points() {
            unsigned int reg_idx = computed_dh_cell_.elm().region_idx().idx();
            for (auto p : mass_eval->points(this->position_in_cache(computed_dh_cell_.elm_idx()), this) ) {
                this->eval_point_data_.emplace_back(reg_idx, computed_dh_cell_.elm_idx(), p.eval_point_idx(), computed_dh_cell_.local_idx());
            }

            for (DHCellSide cell_side : computed_dh_cell_.side_range()) {
            	for( DHCellSide edge_side : cell_side.edge_sides() ) {
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : side_eval->points(edge_side, this) ) {
                        this->eval_point_data_.emplace_back(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    }
                }
            }
            this->eval_point_data_.make_permanent();
        }

        void update_cache() {
            this->register_eval_points();
            this->create_patch();
            this->cache_update(*this);
            this->finish_elements_update();
        }

        void reallocate_cache() {
            this->cache_reallocate(*this, *this);
        }


        // fields
        Field<3, FieldValue<3>::Scalar > scalar_field;                 ///< Coordinate 'x' or 'y', reference value of other fields
        Field<3, FieldValue<3>::Scalar > scalar_z;                     ///< Coordinate 'z', reference of depth computing
        Field<3, FieldValue<3>::Scalar > scalar_with_depth;            ///< Tests depth 'd' value
        Field<3, FieldValue<3>::VectorFixed > vector_field;            ///< Tests formula vector
        Field<3, FieldValue<3>::VectorFixed > density_unit_conversion; ///< Tests unit conversion
        Field<3, FieldValue<3>::TensorFixed > tensor_field;            ///< Tests formula tensor
        Field<3, FieldValue<3>::Scalar > const_scalar;                 ///< Tests field dependency
        Field<3, FieldValue<0>::Integer > integer_scalar;
        std::shared_ptr<EvalPoints> eval_points_;
        std::shared_ptr<BulkIntegral> mass_eval;
        std::shared_ptr<EdgeIntegral> side_eval;
        //std::shared_ptr<CouplingIntegral> ngh_side_eval;
        DHCellAccessor computed_dh_cell_;
    };

    FieldEvalFormulaTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        data_ = std::make_shared<EqData>();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/cube_2x1.msh\", optimize_mesh=false }");
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~FieldEvalFormulaTest() {}

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalFormulaTest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("density_unit_conversion", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("scalar_z", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("scalar_with_depth", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("tensor_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
                        .declare_key("const_scalar", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("integer_scalar", FieldAlgorithmBase< 3, FieldValue<0>::Integer >::get_input_type_instance(), "" )
                        .close()
                        ), IT::Default::obligatory(), ""  )
                .close();
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        TimeGovernor tg(0.0, 1.0);

        //data.set_components(component_names);        // set number of substances posibly read from elsewhere

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        data_->set_mesh(*mesh_);
        data_->set_input_list( inputs[input_last], tg );
        data_->set_time(tg.step(), LimitSide::right);
    }


    std::shared_ptr<EqData> data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};


TEST_F(FieldEvalFormulaTest, evaluate) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldFormula
          value: x
        scalar_z: !FieldFormula
          value: z
        scalar_with_depth: !FieldFormula
          value: d
          surface_region: ".2D top"
        vector_field: !FieldFormula
          value: "[x, 2*x, 0.5]"
        density_unit_conversion: !FieldFormula
          value: "[x, x**2, 2*x+t]"
          unit: g*cm^-3
        tensor_field: !FieldFormula
          value: "[ [x, 0.2, 0.3], [0.2, 0.4, 0.5], [0.3, 0.5, 0.6] ]"
        integer_scalar: 1
      - region: 3D right
        time: 0.0
        scalar_field: !FieldFormula
          value: y
        scalar_z: !FieldFormula
          value: z
        scalar_with_depth: !FieldFormula
          value: d
          surface_region: ".2D top"
        vector_field:  !FieldFormula
          value: "[y, 2*y, 0.5]"
        density_unit_conversion: !FieldFormula
          value: "[y, y**2, 2*y+t]"
          unit: g*cm^-3
        tensor_field: !FieldFormula
          value: "[ [y, 2.2, 2.3], [2.2, 2.4, 2.5], [2.3, 2.5, 2.6] ]"
        integer_scalar: 1
    )YAML";
	this->read_input(eq_data_input);
    data_->reallocate_cache();

    std::vector<unsigned int> cell_idx = {3, 4, 5, 9};
    std::vector<arma::mat33>  expected_tensor = {{0, 0, 0, 0, 0, 0, 0, 0, 0},
                                                 {0.1, 0.2, 0.3, 0.2, 0.4, 0.5, 0.3, 0.5, 0.6},  //region 1
                                                 {0, 0, 0, 0, 0, 0, 0, 0, 0},
                                                 {2.1, 2.2, 2.3, 2.2, 2.4, 2.5, 2.3, 2.5, 2.6}};  //region 3
    for (uint i=0; i<cell_idx.size(); ++i) {
        DebugOut() << "TEST CELL: i=" << i;
        uint test_point = 0; // index to expected vals
    	data_->start_elements_update();
    	data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), cell_idx[i]);  // element ids stored to cache: (3 -> 2,3,4), (4 -> 3,4,5,10), (5 -> 0,4,5,11), (10 -> 8,9,10)
        data_->update_cache();

        uint r_idx = data_->computed_dh_cell_.elm().region().idx(); // element regions: {1,1,1,3} for elements {3,4,5,9}

        // Bulk integral, no sides.
        for( BulkPoint q_point: data_->mass_eval->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()) ) {
            double coord = data_->scalar_field(q_point); // X coord on reg 1,  Y coord on reg 3

            double depth = data_->scalar_with_depth(q_point);
            EXPECT_DOUBLE_EQ(1-data_->scalar_z(q_point), depth);

            // Evaluation of the scalar field tested elseewhere. we can
            // Can activqte again only after we support taking values from FieldCommon.

//            auto coords = data_->X_(q_point);
//            double expect_scalar;
//            if (r_idx == 1) expect_scalar = coords[0];
//            else expect_scalar = coords[1];
//            EXPECT_DOUBLE_EQ( expect_scalar, coord);

            arma::vec3 expected_vector;
            expected_vector(0) = coord;
            expected_vector(1) = 2*coord;
            expected_vector(2) = 0.5;
            EXPECT_ARMA_EQ(expected_vector, data_->vector_field(q_point));

            auto density = data_->density_unit_conversion(q_point);
            expected_vector(0) = 1000*coord;
            expected_vector(1) = expected_vector(0)*coord;
            expected_vector(2) = 2*expected_vector(0);
            EXPECT_ARMA_EQ(expected_vector, data_->density_unit_conversion(q_point));

            auto exp_tensor = expected_tensor[r_idx];
            exp_tensor(0,0) = coord;
            EXPECT_ARMA_EQ(exp_tensor, data_->tensor_field(q_point));
            test_point++;
            /* // Extracting the cached values.
            double cs = cross_section(q_point);

            // Following would be nice to have. Not clear how to
            // deal with more then single element as fe_values have its own cache that has to be updated.
            auto base_fn_grad = presssure_field_fe.base_value(q_point);
            loc_matrix += outer_product((cs * base_fn_grad),  base_fn_grad) */
        }
//        depth /= 4;
//        std::cout << "Element: " << cell_idx[i] << ", average value: " << depth << std::endl;

        // Side integrals.
        // FieldFE<..> conc;
        for (DHCellSide side : data_->computed_dh_cell_.side_range()) {
        	for(DHCellSide edg_side : side.edge_sides()) {
           	    // vector of local side quadrature points
        	    Range<EdgePoint> side_points = data_->side_eval->points(side, data_.get());
        	    for (EdgePoint side_p : side_points) {

        	        //uint r_idx = edg_side.element().region().idx();
        	        //DebugOut() << "ele region: " << r_idx;

        	        double coord = data_->scalar_field(side_p);
                    //EXPECT_DOUBLE_EQ( dbl_scalar[ expected_scalar[i][test_point] ], data_->scalar_field(side_p));
                    arma::vec3 expected_vector;
                    expected_vector(0) = coord;
                    expected_vector(1) = 2*coord;
                    expected_vector(2) = 0.5;
                    EXPECT_ARMA_EQ(expected_vector, data_->vector_field(side_p));

                    auto exp_tensor = expected_tensor[r_idx];
                    exp_tensor(0,0) = coord;
                    EXPECT_ARMA_EQ(exp_tensor, data_->tensor_field(side_p));
                    test_point++;
                    EdgePoint ngh_p = side_p.point_on(edg_side);
                    //EXPECT_DOUBLE_EQ( coord, data_->scalar_field(ngh_p));
                    test_point++;
        	        //loc_mat += cross_section(side_p) * sigma(side_p) *
        		    //    (conc.base_value(side_p) * velocity(side_p)
        		    //    + conc.base_value(ngh_p) * velocity(ngh_p)) * side_p.normal() / 2;
                }
            }
        }
    }

}

TEST_F(FieldEvalFormulaTest, field_dependency) {
    string eq_data_input = R"YAML(
    data:
      - region: 3D left
        time: 0.0
        scalar_field: !FieldFormula
          value: const_scalar * const_scalar
        vector_field: !FieldFormula
          value: "[scalar_field, 2*scalar_field, 0.5]"
        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        const_scalar: 0.5
        integer_scalar: 1
      - region: 3D right
        time: 0.0
        scalar_field: !FieldFormula
          value: 0.5 * const_scalar
        vector_field:  !FieldFormula
          value: "[scalar_field, 2*scalar_field, 0.5]"
        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        const_scalar: 1.5
        integer_scalar: 1
    )YAML";
    this->read_input(eq_data_input);
    data_->reallocate_cache();

    std::vector<unsigned int> cell_idx = {3, 4, 5, 9};
    // scalar_value 0.25 on region 1, scalar_value 0.75 on region 3
    std::vector<double> region_value = { 0, 0.25, 0, 0.75 };

    for (uint i=0; i<cell_idx.size(); ++i) {
        data_->start_elements_update();
        data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), cell_idx[i]);  // element ids stored to cache: (3 -> 2,3,4), (4 -> 3,4,5,10), (5 -> 0,4,5,11), (10 -> 8,9,10)
        data_->update_cache();

        double expected_val = region_value[data_->computed_dh_cell_.elm().region().idx()];
        arma::vec3 expected_vector;
        expected_vector(0) = expected_val;
        expected_vector(1) = 2*expected_val;
        expected_vector(2) = 0.5;

        // Bulk integral, no sides.
        for( BulkPoint q_point: data_->mass_eval->points(data_->position_in_cache(data_->computed_dh_cell_.elm_idx()), data_.get()) ) {
            EXPECT_DOUBLE_EQ( expected_val, data_->scalar_field(q_point));
            EXPECT_ARMA_EQ(expected_vector, data_->vector_field(q_point));
        }

        // Side integrals.
        // FieldFE<..> conc;
        for (DHCellSide side : data_->computed_dh_cell_.side_range()) {
            for(DHCellSide edg_side : side.edge_sides()) {
                //DebugOut() << "ele region: " << edg_side.element().region().idx();
                // vector of local side quadrature points
           	    Range<EdgePoint> side_points = data_->side_eval->points(side, data_.get());
           	    double expected_edg_side_val = region_value[edg_side.element().region().idx()];
                for (EdgePoint side_p : side_points) {

                    EXPECT_DOUBLE_EQ( expected_val, data_->scalar_field(side_p));
                    EXPECT_ARMA_EQ(expected_vector, data_->vector_field(side_p));
                    EdgePoint ngh_p = side_p.point_on(edg_side);
                    EXPECT_DOUBLE_EQ( expected_edg_side_val, data_->scalar_field(ngh_p));

                    //DebugOut() << "el side val: " << data_->scalar_field(side_p) << ", edge side val: " << data_->scalar_field(ngh_p) << "";
                }
            }
        }
    }

}


TEST_F(FieldEvalFormulaTest, dependency_unknown_field_exc) {
    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldFormula
          value: 2 * unknown_scalar
        vector_field: [0.1, 0.2, 0.3]
        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        const_scalar: 0.5
        integer_scalar: 1
    )YAML";
    this->read_input(eq_data_input);
    EXPECT_THROW_WHAT( { data_->reallocate_cache(); }, FieldSet::ExcUnknownField,
        "Unknown field 'unknown_scalar' in the formula");
}


TEST_F(FieldEvalFormulaTest, dependency_notdouble_field_exc) {
	typedef FieldFormula<3, FieldValue<3>::Scalar > ScalarFormula;

    string eq_data_input = R"YAML(
    data:
      - region: BULK
        time: 0.0
        scalar_field: !FieldFormula
          value: 2 * integer_scalar
        vector_field: [0.1, 0.2, 0.3]
        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        const_scalar: 0.5
        integer_scalar: 1
    )YAML";
    this->read_input(eq_data_input);
    EXPECT_THROW_WHAT( { data_->reallocate_cache(); }, ScalarFormula::ExcNotDoubleField,
        "Can not use integer valued field 'integer_scalar' in the formula");
}


// TODO Prepared test of BParser::Exception
//  - uncomment after removing of FParser
//  - use EXPECT_THROW_WHAT instead of try-catch
//TEST_F(FieldEvalFormulaTest, formula_invalid_expr) {
//    string eq_data_input = R"YAML(
//    data:
//      - region: BULK
//        time: 0.0
//        scalar_field: !FieldFormula
//          value: 2 * max(3, 5
//        vector_field: [0.1, 0.2, 0.3]
//        tensor_field: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
//        const_scalar: 0.5
//        integer_scalar: 1
//    )YAML";
//    this->read_input(eq_data_input);
//    try {
//        data_->reallocate_cache();
//    } catch (Input::Exception &e) {
//        std::string s( e.what() );
//        int pos = (int)s.find("Error: Expected \")\" at");
//        EXPECT_GE(pos, 0); // substring "Error: Expected ..." found
//    }
//}



string set_time_input = R"INPUT(
[ 
      { TYPE="FieldFormula",  value="[x, x*y, y+t]" },
      { TYPE="FieldFormula",  value="[x, x*y, y]" },
      { TYPE="FieldFormula",  value="[x+t, x*y+t, y+t]" },
      { TYPE="FieldFormula",  value="[x, x*y, y]" }
]

)INPUT";

TEST(FieldFormula, set_time) {
    typedef FieldAlgorithmBase<3, FieldValue<3>::VectorFixed > VectorField;

    Profiler::instance();

    // setup FilePath directories
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

    Input::Type::Array  input_type(VectorField::get_input_type_instance());
    input_type.finish();

    // read input string
    Input::ReaderToStorage reader( set_time_input, input_type, Input::FileFormat::format_JSON );
    Input::Array in_array=reader.get_root_interface<Input::Array>();

    auto it = in_array.begin<Input::AbstractRecord>();
    FieldAlgoBaseInitData init_data("formula", 3, UnitSI::dimensionless());

    {
        auto field=VectorField::function_factory(*it, init_data);
        EXPECT_TRUE( field->set_time(1.0) );
        EXPECT_TRUE( field->set_time(2.0) );
    }
    ++it;

    {
        auto field=VectorField::function_factory(*it, init_data);
        TimeGovernor tg(3.0, 1.0);
        auto step0 = tg.step();
        EXPECT_TRUE( field->set_time(step0) );
        tg.next_time();
        auto step1 = tg.step();
        EXPECT_TRUE( field->set_time(step1) ); // should be false, needs revision of FieldFormula::set_time return value
    }
    ++it;

    {
        auto field=VectorField::function_factory(*it, init_data);
        EXPECT_TRUE( field->set_time(1.5) );
        EXPECT_TRUE( field->set_time(2.5) );
    }
    ++it;

    {
        auto field=VectorField::function_factory(*it, init_data);
        TimeGovernor tg(0.0, 2.0);
        auto step0 = tg.step();
        EXPECT_TRUE( field->set_time(step0) );
        tg.next_time();
        auto step1 = tg.step();
        EXPECT_TRUE( field->set_time(step1) ); // should be false, needs revision of FieldFormula::set_time return value
    }

}


TEST(SurfaceDepth, base_test) {
    Profiler::instance();
	FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");

	std::string mesh_in_string = "{mesh_file=\"fields/surface_reg.msh\"}";
	Mesh * mesh = mesh_constructor(mesh_in_string);
    auto reader = reader_constructor(mesh_in_string);
	reader->read_physical_names(mesh);
	reader->read_raw_mesh(mesh);

	SurfaceDepth sd(mesh, ".top side", "0 0 1");
	EXPECT_DOUBLE_EQ( sd.compute_distance( arma::vec3("1 0.5 -0.9") ), 1.9 );
	EXPECT_DOUBLE_EQ( sd.compute_distance( arma::vec3("-1 0.5 0.9") ), 0.1 );

	delete mesh;
}
