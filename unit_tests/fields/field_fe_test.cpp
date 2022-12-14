/*
 * field_fe_test.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: jb
 */



#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include <limits>
#include "arma_expect.hh"


#include "fields/eval_points.hh"
#include "fields/eval_subset.hh"
#include "fields/field_value_cache.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
#include "fields/field_fe.hh"
#include "tools/unit_si.hh"
#include "fields/bc_field.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "input/reader_to_storage.hh"
#include "system/sys_profiler.hh"
#include "la/vector_mpi.hh"
#include "tools/mixed.hh"




//class FieldFETest : public testing::Test {
//public:
//    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
//    typedef FieldFE<3, FieldValue<3>::VectorFixed > VecField;
//
//    virtual void SetUp() {
//    	this->mesh = nullptr;
//        // setup FilePath directories
//        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
//
//        Profiler::instance();
//        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
//    }
//
//    virtual void TearDown() {
//    	dh.reset();
//    	if (mesh != nullptr) delete mesh;
//    }
//
//    void create_mesh(std::string mesh_file_str) {
//        mesh = mesh_full_constructor("{ mesh_file=\"" + mesh_file_str + "\", optimize_mesh=false }");
//    }
//
//    void create_dof_handler(double val1, double val2, double val3) {
//        dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
//        v.resize(3);
//        v.set(0, val1);
//        v.set(1, val2);
//        v.set(2, val3);
//        dof_values[0] = val1;
//        dof_values[1] = val2;
//        dof_values[2] = val3;
//    }
//
//    const FieldAlgoBaseInitData& init_data(std::string field_name) {
//    	static const FieldAlgoBaseInitData init_data(field_name, 0, UnitSI::dimensionless());
//    	return init_data;
//    }
//
//    static Input::Type::Record &get_input_type() {
//        return Input::Type::Record("Test","")
//            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("native_data", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .close();
//    }
//
//    Mesh *mesh;
//    std::shared_ptr<DOFHandlerMultiDim> dh;
//    double dof_values[3];
//    VectorMPI v;
//
//};
//
//
//TEST_F(FieldFETest, scalar) {  // moved to FieldEvalFETest, set_fe_data_scalar
//    create_mesh("fields/one_element_2d.msh");
//    create_dof_handler(1, 2, 3);
//
//	MixedPtr<FE_P_disc> fe(1);
//    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
//    ScalarField field;
//
//    dh->distribute_dofs(ds);
//    field.set_fe_data(dh, v);
//    field.set_time(0.0);
//
//    Armor::array pts(3, 1);
//    pts.reinit(3);
//    pts.append(Armor::vec<3>({ 1, 1, 5 }));
//    pts.append(Armor::vec<3>({ 4, 0, 5 }));
//    pts.append(Armor::vec<3>({ 2, 3, 5 }));
//    vector<double> values(3);
//
//    // test values at vertices of the triangle
//    field.value_list( pts, mesh->element_accessor(0), values );
//    EXPECT_DOUBLE_EQ( dof_values[0], values[0] );
//    EXPECT_DOUBLE_EQ( dof_values[1], values[1] );
//    EXPECT_DOUBLE_EQ( dof_values[2], values[2] );
//
//    // test value at barycenter
//    EXPECT_DOUBLE_EQ( (dof_values[0]+dof_values[1]+dof_values[2])/3, field.value({ 7./3, 4./3, 5 }, mesh->element_accessor(0)) );
//}
//
//
//TEST_F(FieldFETest, vector) {  // moved to FieldEvalFETest, set_fe_data_vector
//    create_mesh("fields/one_element_2d.msh");
//    create_dof_handler(0, 0, 1);
//
//	MixedPtr<FE_RT0> fe;
//    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh, fe);
//    VecField field;
//
//    dh->distribute_dofs(ds);
//    field.set_fe_data(dh, v);
//    field.set_time(0.0);
//
//    // The Raviart-Thomas function given by the following dofs
//    // is 3/7*(x-7/3, y-4/3, 0).
//
//    arma::vec3 result = { 2./7, 1./14, 0 };
//
//    EXPECT_NEAR( 0, arma::norm(result - field.value({ 3, 1.5, 5 }, mesh->element_accessor(0)), 2), 1e-15 );
//}
//
//
//string input = R"INPUT(
//{
//   scalar={
//       TYPE="FieldFE",
//       mesh_data_file="fields/simplest_cube_data.msh",
//       field_name="scalar"
//   }
//   native_data={
//       TYPE="FieldFE",
//       mesh_data_file="output/test_output_vtk_ascii_ref.vtu",
//       field_name="flow_data"
//   }
//}
//)INPUT";
//
//
//
//TEST_F(FieldFETest, scalar_from_input) { // equivalent with TEST(FieldFENewTest, scalar) - moved to FieldEvalFETest, input_msh
//    create_mesh("fields/simplest_cube_data.msh");
//
//    Input::ReaderToStorage reader( input, FieldFETest::get_input_type(), Input::FileFormat::format_JSON );
//    Input::Record rec=reader.get_root_interface<Input::Record>();
//
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
//    field.set_mesh(mesh,false);
//    field.set_time(0.0);
//
//    Space<3>::Point point;
//    for(unsigned int i=0; i < mesh->n_elements(); i++) {
//        EXPECT_DOUBLE_EQ( (i+1)*0.1 , field.value(point, mesh->element_accessor(i)) );
//    }
//}
//
//
//TEST_F(FieldFETest, native_data) { // moved to FieldEvalFETest, native_data
//    create_mesh("fields/simplest_cube_3d.msh");
//
//    Input::ReaderToStorage reader( input, FieldFETest::get_input_type(), Input::FileFormat::format_JSON );
//    Input::Record rec=reader.get_root_interface<Input::Record>();
//
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("native_data"), init_data("native_data"));
//    field.set_mesh(mesh,false);
//    field.set_time(0.0);
//
//    Space<3>::Point point;
//    for(unsigned int i=0; i < mesh->n_elements(); i++) {
//        EXPECT_DOUBLE_EQ( i*0.2 , field.value(point, mesh->element_accessor(i)) );
//    }
//}


/**
 * Main class of unit test.
 *
 * Simulates simplified behavior of assembly algorithm. Class contains instance of inner class EqData
 * (descendant of FieldSet) with instances of different fields (different shapes and boundary or bulk).
 * These fields are not obligatory and only some instances can be used in test cases.
 *
 * Methods used in test cases:
 *  1. create_mesh - creates mesh and appropriate DOF handler
 *  2. read_input - reads user input in YAML format and set it to EqData instance
 *     (alternatively can by used set_dof_values method in tests of FieldFE::set_fe_data method)
 *  3. EqData::reallocate_cache - set time and reallocate cache, method can be called repeatedly
 *     for different times
 *  4. eval_bulk_field - evaluates values on all elements of bulk mesh and compare them with reference
 *     field or value
 *  5. eval_boundary_field - evaluates value on one specified element of boundary mesh and compare them
 *     with reference field or value. Needs set idx of boundary element on boundary mesh and idx of
 *     appropriate bulk element on bulk mesh
 */
class FieldEvalFETest : public testing::Test {
public:
    typedef Field<3, FieldValue<3>::Scalar > ScalarField;
    typedef Field<3, FieldValue<3>::Enum > EnumField;
    typedef Field<3, FieldValue<3>::VectorFixed > VectorField;
    typedef Field<3, FieldValue<3>::TensorFixed > TensorField;
    typedef BCField<3, FieldValue<3>::Scalar > BcScalarField;
    typedef BCField<3, FieldValue<3>::Enum > BcEnumField;
    typedef BCField<3, FieldValue<3>::VectorFixed > BcVectorField;
    typedef BCField<3, FieldValue<3>::TensorFixed > BcTensorField;

    class EqData : public FieldSet, public ElementCacheMap {
    public:
        enum enum_type {
            none=0,
            dirichlet=1,
            neumann=2,
            robin=3,
            total_flux=4
        };

        static const IT::Selection & get_enum_selection() {
        	return IT::Selection("EqData_enum_Type")
                     .add_value(none, "none")
                     .add_value(dirichlet, "dirichlet")
                     .add_value(neumann, "neumann")
                     .add_value(robin, "robin")
                     .add_value(total_flux, "total_flux")
        			 .close();
        }

        /**
         * Descendant of FieldSet. Contains different fields for evaluation and referenced.
         */
        EqData() : tg_(0.0, 1.0) {
            *this += scalar_field
                        .name("scalar_field")
                        .description("Scalar field.")
						.input_default("0.0")
                        .units( UnitSI().m() );
            *this += vector_field
                        .name("vector_field")
                        .description("Vector field.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += tensor_field
                        .name("tensor_field")
                        .description("Tensor field.")
						.input_default("0.0")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += enum_field
                        .name("enum_field")
                        .description("Enum field.")
						.input_selection( get_enum_selection() )
						.input_default("\"none\"")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += bc_scalar_field
                        .name("bc_scalar_field")
                        .description("Boundary scalar field.")
						.input_default("0.0")
                        .units( UnitSI().m() );
            *this += bc_vector_field
                        .name("bc_vector_field")
                        .description("Boundary vector field.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += bc_tensor_field
                        .name("bc_tensor_field")
                        .description("Boundary tensor field.")
						.input_default("0.0")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += bc_enum_field
                        .name("bc_enum_field")
                        .description("Boundary enum field.")
						.input_selection( get_enum_selection() )
						.input_default("\"none\"")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += scalar_ref
                        .name("scalar_ref")
                        .description("Reference scalar field.")
						.input_default("0.0")
                        .units( UnitSI().m() );
            *this += vector_ref
                        .name("vector_ref")
                        .description("Reference vector field.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += tensor_ref
                        .name("tensor_ref")
                        .description("Reference tensor field.")
						.input_default("0.0")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);
            *this += bc_scalar_ref
                        .name("bc_scalar_ref")
                        .description("Reference boundary scalar field.")
						.input_default("0.0")
                        .units( UnitSI().m() );
            *this += bc_vector_ref
                        .name("bc_vector_ref")
                        .description("Reference boundary vector field.")
                        .input_default("0.0")
                        .flags_add(in_main_matrix)
                        .units( UnitSI().kg(3).m() );
            *this += bc_tensor_ref
                        .name("bc_tensor_ref")
                        .description("Reference boundary tensor field.")
						.input_default("0.0")
                        .units( UnitSI::dimensionless() )
                        .flags_add(in_main_matrix);

            // Asumme following types:
            eval_points_ = std::make_shared<EvalPoints>();
            Quadrature *q_bulk1 = new QGauss(1, 0);
            Quadrature *q_bulk2 = new QGauss(2, 0);
            Quadrature *q_bulk3 = new QGauss(3, 0);
            Quadrature *q_bdr1 = new QGauss(0, 0);
            Quadrature *q_bdr2 = new QGauss(1, 0);
            Quadrature *q_bdr3 = new QGauss(2, 0);
            mass_integral[0] = eval_points_->add_bulk<1>(*q_bulk1 );
            mass_integral[1] = eval_points_->add_bulk<2>(*q_bulk2 );
            mass_integral[2] = eval_points_->add_bulk<3>(*q_bulk3 );
            bdr_integral[0] = eval_points_->add_boundary<1>(*q_bdr1 );
            bdr_integral[1] = eval_points_->add_boundary<2>(*q_bdr2 );
            bdr_integral[2] = eval_points_->add_boundary<3>(*q_bdr3 );
            this->init(eval_points_);

            this->add_coords_field();
            this->set_default_fieldset();
        }

        void register_eval_points(bool bdr=false) {
            unsigned int reg_idx = computed_dh_cell_.elm().region_idx().idx();
            for (auto p : mass_integral[computed_dh_cell_.dim()-1]->points(this->position_in_cache(computed_dh_cell_.elm_idx()), this) ) {
                this->eval_point_data_.emplace_back(reg_idx, computed_dh_cell_.elm_idx(), p.eval_point_idx(), computed_dh_cell_.local_idx());
            }

            if (bdr)
                for (DHCellSide cell_side : computed_dh_cell_.side_range()) {
                    if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                        unsigned int bdr_reg = cell_side.cond().element_accessor().region_idx().idx(); // region of boundary element
                        //DebugOut() << "Bulk elm: " << computed_dh_cell_.elm_idx() << ", boundary element: " << cell_side.cond().bc_ele_idx() << ", region: " << bdr_reg << std::endl;
                        for (auto p : bdr_integral[computed_dh_cell_.dim()-1]->points(cell_side, this) ) {
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
            this->clear_element_eval_points_map();
            this->start_elements_update();
            this->register_eval_points(bdr);
            this->create_patch();
            this->cache_update(*this);
            this->finish_elements_update();
        }

        void reallocate_cache() {
            this->set_time(tg_.step(), LimitSide::right);
            this->cache_reallocate(*this, *this);
        }


        // fields
        ScalarField scalar_field;
        VectorField vector_field;
        TensorField tensor_field;
        EnumField enum_field;
        BcScalarField bc_scalar_field;
        BcVectorField bc_vector_field;
        BcTensorField bc_tensor_field;
        BcEnumField bc_enum_field;
        ScalarField scalar_ref;
        VectorField vector_ref;
        TensorField tensor_ref;
        BcScalarField bc_scalar_ref;
        BcVectorField bc_vector_ref;
        BcTensorField bc_tensor_ref;

        std::shared_ptr<EvalPoints> eval_points_;
        std::array< std::shared_ptr<BulkIntegral>, 3> mass_integral;
        std::array< std::shared_ptr<BoundaryIntegral>, 3> bdr_integral;
        DHCellAccessor computed_dh_cell_;
        TimeGovernor tg_;
    };

    FieldEvalFETest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        eq_data_ = std::make_shared<EqData>();
    }

    ~FieldEvalFETest() {}

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalFETest::EqData().make_field_descriptor_type("SomeEquation") )
                        .declare_key("scalar_field", FieldAlgorithmBase< 3, FieldValue<3>::Scalar >::get_input_type_instance(), "" )
                        .declare_key("vector_field", FieldAlgorithmBase< 3, FieldValue<3>::VectorFixed >::get_input_type_instance(), "" )
                        .declare_key("tensor_field", FieldAlgorithmBase< 3, FieldValue<3>::TensorFixed >::get_input_type_instance(), "" )
                        .close()
                        ), IT::Default::obligatory(), ""  )
                .close();
    }

    void create_mesh(const std::string &mesh_file) {
        std::string input_str = "{ mesh_file=\"" + mesh_file + "\", optimize_mesh=false }";
        mesh_ = mesh_full_constructor(input_str);
        dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    /// Specialization of set_fe_data unit tests
    void set_dof_values(std::vector<double> vals) {
        v.resize(vals.size());
        dof_values.resize(vals.size());
        for (unsigned int i=0; i<vals.size(); ++i) {
            v.set(i, vals[i]);
            dof_values[i] = vals[i];
        }

        eq_data_->set_mesh(*mesh_);
    }

    void read_input(const string &input) {
        // read input string
        Input::ReaderToStorage reader( input, get_input_type(), Input::FileFormat::format_YAML );
        Input::Record in_rec=reader.get_root_interface<Input::Record>();

        //data.set_components(component_names);        // set number of substances posibly read from elsewhere

        static std::vector<Input::Array> inputs;
        unsigned int input_last = inputs.size(); // position of new item
        inputs.push_back( in_rec.val<Input::Array>("data") );

        eq_data_->set_mesh(*mesh_);
        eq_data_->set_input_list( inputs[input_last], eq_data_->tg_ );
    }


    /// Compares field value with ref_val. Internal method, called only from check_point_value method.
    template<class Val>
    bool compare_vals(const Val &field_val, const Val &ref_val);


    /// Internal method of eval_bulk_field and eval_boundary_field
    template<class EvalField, class RefVal>
    bool check_point_value(EvalField &field, BulkPoint &point, const RefVal &ref_val)
    {
	    try {
            auto field_val = field( point );
            return this->compare_vals(field_val, ref_val);
        } catch (ExceptionBase &e) {
            std::cout << e.what() << std::endl;
            return false;
        }
        return true;
    }


    /// Evaluates and compare bulk field, uses FieldFormula as reference
    template<class FieldType>
    bool eval_bulk_field(FieldType &eval_field, FieldType &ref_field)
    {
        bool is_passed = true;
	    for(unsigned int i=0; i < mesh_->n_elements(); i++) {
            eq_data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), i);
            uint dim = eq_data_->computed_dh_cell_.dim()-1;
            eq_data_->update_cache();

            auto p = *( eq_data_->mass_integral[dim]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );
            is_passed = is_passed & check_point_value(eval_field, p, ref_field(p) );
        }
	    return is_passed;
    }


    /// Evaluates and compare bulk field, uses vector of value as reference
    template<class FieldType, class RefVal>
    bool eval_bulk_field(FieldType &eval_field, std::vector<RefVal> ref_val)
    {
        static_assert( std::is_same<typename FieldType::ValueType::return_type, RefVal>::value, "Wrong type of RefVal.");
        ASSERT_EQ( mesh_->n_elements(), ref_val.size() );

        bool is_passed = true;
	    for(unsigned int i=0; i < mesh_->n_elements(); i++) {
            eq_data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), i);
            uint dim = eq_data_->computed_dh_cell_.dim()-1;
            eq_data_->update_cache();

            auto p = *( eq_data_->mass_integral[dim]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );
            is_passed = is_passed & check_point_value(eval_field, p, ref_val[i] );
        }
	    return is_passed;
    }


    /// Evaluates and compare bulk field, uses single value as reference on all evaluated elements
    template<class FieldType, class RefVal>
    bool eval_bulk_field(FieldType &eval_field, RefVal ref_val)
    {
    	std::vector<RefVal> ref_val_vec(mesh_->n_elements(), ref_val);
	    return eval_bulk_field(eval_field, ref_val_vec);
    }


    /// Evaluates and compare boundary field, uses FieldFormula as reference
    template<class FieldType>
    bool eval_boundary_field(FieldType &eval_field, FieldType &ref_field, unsigned int i_bulk_elem, unsigned int i_bdr_elem)
    {
        eq_data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), i_bulk_elem);
        uint dim = eq_data_->computed_dh_cell_.dim()-1;
        eq_data_->update_cache(true);

        for (DHCellSide cell_side : eq_data_->computed_dh_cell_.side_range()) {
            if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                if (cell_side.cond().bc_ele_idx() == i_bdr_elem) {
                    auto p_side = *( eq_data_->bdr_integral[dim]->points(cell_side, eq_data_.get()).begin() );
                    auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
                    //DebugOut() << "Input_id: " << cell_side.cond().element_accessor().input_id() << ", value: " << eval_field(p_bdr) << std::endl;
                    // Boundary found - field is evaluated and checked
                    return check_point_value(eval_field, p_bdr, ref_field(p_bdr));
                }
            }
        }
        // Boundary not found - test doesn'ù pass
        std::cout << "Wrong indices pair of bulk and boundary element!" << std::endl;
        std::cout << "Boundary element of idx: " << i_bdr_elem << " is not relevant to bulk element of idx: " << i_bulk_elem << std::endl;
	    return false;
    }


    /// Evaluates and compare boundary field, uses value as reference
    template<class FieldType, class RefVal>
    bool eval_boundary_field(FieldType &eval_field, RefVal ref_val, unsigned int i_bulk_elem, unsigned int i_bdr_elem)
    {
        static_assert( std::is_same<typename FieldType::ValueType::return_type, RefVal>::value, "Wrong type of RefVal.");

        eq_data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), i_bulk_elem);
        uint dim = eq_data_->computed_dh_cell_.dim()-1;
        eq_data_->update_cache(true);

        for (DHCellSide cell_side : eq_data_->computed_dh_cell_.side_range()) {
            if ( (cell_side.side().edge().n_sides() == 1) && (cell_side.side().is_boundary()) ) {
                if (cell_side.cond().bc_ele_idx() == i_bdr_elem) {
                    auto p_side = *( eq_data_->bdr_integral[dim]->points(cell_side, eq_data_.get()).begin() );
                    auto p_bdr = p_side.point_bdr( cell_side.cond().element_accessor() );
                    // Boundary found - field is evaluated and checked
                    return check_point_value(eval_field, p_bdr, ref_val);
                }
            }
        }
        // Boundary not found - test doesn'ù pass
        std::cout << "Wrong indices pair of bulk and boundary element!" << std::endl;
        std::cout << "Boundary element of idx: " << i_bdr_elem << " is not relevant to bulk element of idx: " << i_bulk_elem << std::endl;
	    return false;
    }


    std::shared_ptr<EqData> eq_data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    std::vector<double> dof_values;           ///< used in test set_fe_data
    VectorMPI v;                              ///< used in test set_fe_data
};


/// General function compares value of Vector or Tensor fields with ref_val
template<class Val>
bool FieldEvalFETest::compare_vals(const Val &field_val, const Val &ref_val)
{
    try {
        auto ref_shape = mat_shape(ref_val);
        auto shape = mat_shape(field_val);
        if (ref_shape[0] != shape[0] || ref_shape[1] != shape[1]) {
            if (ref_shape[0] != shape[0]) std::cout << "Different number of rows of field value and ref value!" << std::endl;
            if (ref_shape[1] != shape[1]) std::cout << "Different number of cols of field value and ref value!" << std::endl;
            return false;
        }
        double magnitude = std::max( arma::norm(ref_val, 1), arma::norm(field_val, 1) );
        // abs criterium
        if (magnitude < 8*std::numeric_limits<double>::epsilon()) return true;

        double error = arma::norm(ref_val - field_val, 1)/magnitude;
        // rel criterium
        if (error > 8*std::numeric_limits<double>::epsilon()) {
            unsigned int w = 11* field_val.n_cols;
            std::cout << std::setw(w) << "Expected" << std::setw(w) << "Result"
                      << " rel. error: " << error << std::endl;
            for(unsigned int i_row = 0; i_row < field_val.n_rows; i_row++) {
                std::cout << std::setw(11) << ref_val.row(i_row)
                          << std::setw(11) << field_val.row(i_row)
                          << std::setw(20) << ref_val.row(i_row) - field_val.row(i_row) << std::endl;
            }
            return false;
        }
    } catch (ExceptionBase &e) {
        std::cout << e.what() << std::endl;
        return false;
    }
    return true;
}


/// Template specialization for Enum fields.
template<>
bool FieldEvalFETest::compare_vals(const unsigned int &field_val, const unsigned int &ref_val)
{
    try {
        if ( ref_val != field_val ) {
            std::cout << "\nEvaluated val: " << field_val << std::endl;
            std::cout << "Expected val:  " << ref_val << std::endl;
            return false;
        }
    } catch (ExceptionBase &e) {
        std::cout << e.what() << std::endl;
        return false;
    }
    return true;
}


/// Template specialization for Scalar fields.
template<>
bool FieldEvalFETest::compare_vals(const double &field_val, const double &ref_val)
{
    try {
        if ( std::abs(ref_val - field_val) > 8*std::numeric_limits<double>::epsilon() ) {
            std::cout << "\nEvaluated val: " << field_val << std::endl;
            std::cout << "Expected val:  " << ref_val << std::endl;
            return false;
        }
    } catch (ExceptionBase &e) {
        std::cout << e.what() << std::endl;
        return false;
    }
    return true;
}





TEST_F(FieldEvalFETest, input_msh) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          default_value: 0.0
        vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
          default_value: 0.0
        tensor_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: tensor_fixed
          default_value: 0.0
        enum_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: enum
          default_value: 0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          default_value: 0.0
        bc_vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: vector_fixed
          default_value: 0.0
        bc_tensor_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: tensor_fixed
          default_value: 0.0
        bc_enum_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: enum
          default_value: 0
        scalar_ref: !FieldFormula
          value: x+2*y+t
        vector_ref: !FieldFormula
          value: [x+2*y, y+2*z+0.5*t, z+2*x+t]
        tensor_ref: !FieldFormula
          value: [2*x+y, 2*y+z+0.5*t, 2*z+x+t]
        bc_scalar_ref: !FieldFormula
          value: x+y+z+t
        bc_vector_ref: !FieldFormula
          value: [3*x, 3*y+0.5*t, 3*z+t]
        bc_tensor_ref: !FieldFormula
          value: [x+y+z, 2*(x+y+z)+0.5*t, 3*(x+y+z)+t]
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
    	eq_data_->reallocate_cache();

        // BULK fields
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, eq_data_->scalar_ref) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, eq_data_->vector_ref) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, eq_data_->tensor_ref) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->enum_field, j) );

        // BOUNDARY fields
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, eq_data_->bc_scalar_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, eq_data_->bc_vector_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_tensor_field, eq_data_->bc_tensor_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_enum_field, j+1, 3, 0) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, input_vtk) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/vtk_ascii_data.vtu
          field_name: scalar_field
        vector_field: !FieldFE
          mesh_data_file: fields/vtk_ascii_data.vtu
          field_name: vector_field
        tensor_field: !FieldFE
          mesh_data_file: fields/vtk_ascii_data.vtu
          field_name: tensor_field
        scalar_ref: !FieldFormula
          value: x+2*y
        vector_ref: !FieldFormula
          value: [x+2*y, y+2*z, z+2*x]
        tensor_ref: !FieldFormula
          value: [2*x+y, 2*y+z, 2*z+x]
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    eq_data_->reallocate_cache();
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, eq_data_->scalar_ref) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, eq_data_->vector_ref) );
    EXPECT_TRUE( eval_bulk_field(eq_data_->tensor_field, eq_data_->tensor_ref) );
}


TEST_F(FieldEvalFETest, time_shift) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          default_value: 0.0
          read_time_shift: 1.0
        scalar_ref: !FieldFormula
          value: x+2*y+(t+1)
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, eq_data_->scalar_ref) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, default_values) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        vector_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: unset_bulk
          default_value: 0.1
    )YAML";

    string expected_vector = "0.1 0.1 0.1";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, arma::vec3(expected_vector)) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, unit_conversion) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          unit: "const; const=0.1*m"
          default_value: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/simplest_cube_data.msh
          field_name: scalar
          unit: "const; const=0.1*m"
          default_value: 0.0
        scalar_ref: !FieldFormula
          value: 0.1*(x+2*y+t)
        bc_scalar_ref: !FieldFormula
          value: 0.1*(x+y+z+t)
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();

        // BULK field
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, eq_data_->scalar_ref) );

        // BOUNDARY field
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, eq_data_->bc_scalar_ref, 3, 0) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, identic_mesh) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: scalar
          default_value: 0.0
          interpolation: identic_mesh
        vector_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: vector_fixed
          default_value: 0.0
          interpolation: identic_mesh
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: scalar
          default_value: 0.0
          interpolation: identic_mesh
        bc_vector_field: !FieldFE
          mesh_data_file: fields/identic_mesh_data.msh
          field_name: vector_fixed
          default_value: 0.0
          interpolation: identic_mesh
        scalar_ref: !FieldFormula
          value: x+2*y+t
        vector_ref: !FieldFormula
          value: [x+2*y, y+2*z+0.5*t, z+2*x+t]
        bc_scalar_ref: !FieldFormula
          value: x+y+z+t
        bc_vector_ref: !FieldFormula
          value: [3*x, 3*y+0.5*t, 3*z+t]
    )YAML";

    this->create_mesh("mesh/simplest_cube.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();

        // BULK field
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, eq_data_->scalar_ref) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, eq_data_->vector_ref) );

        // BOUNDARY field
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, eq_data_->bc_scalar_ref, 3, 0) );
        EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, eq_data_->bc_vector_ref, 3, 0) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, interpolation_gauss) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: scalar
          default_value: 0.0
          #interpolation: P0_gauss
        vector_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: vector_fixed
          default_value: 0.0
          #interpolation: P0_gauss
    )YAML";

    std::vector< std::vector<double> >     expected_scalars = { {0.25, 0.15, 0.25, 0.35}, {0.75, 0.65, 0.75, 0.85} };
    std::vector< std::vector<arma::vec3> > expected_vectors = { {"2.5 3.5 4.5", "1.5 2.5 3.5", "2.5 3.5 4.5", "3.5 4.5 5.5"},
                                                                {"5.5 6.5 7.5", "4.5 5.5 6.5", "5.5 6.5 7.5", "6.5 7.5 8.5"} };

    this->create_mesh("fields/interpolation_rect_small.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {  // time loop
        eq_data_->reallocate_cache();
        EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, expected_scalars[j]) );
        EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, expected_vectors[j]) );
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, interpolation_1d_2d) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        bc_scalar_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: scalar
          default_value: 0.0
          #interpolation: P0_intersection
        bc_vector_field: !FieldFE
          mesh_data_file: fields/interpolation_rectangle.msh
          field_name: vector_fixed
          default_value: 0.0
          #interpolation: P0_intersection
    )YAML";

    std::vector< std::vector<double> >     expected_scalars = { {0.25, 0.15, 0.25, 0.35}, {0.75, 0.65, 0.75, 0.85} };
    std::vector< std::vector<arma::vec3> > expected_vectors = {{"2.5 3.5 4.5", "1.5 2.5 3.5", "2.5 3.5 4.5", "3.5 4.5 5.5"},
                                                               {"5.5 6.5 7.5", "4.5 5.5 6.5", "5.5 6.5 7.5", "6.5 7.5 8.5"} };

    this->create_mesh("fields/interpolation_rect_small.msh");
    this->read_input(eq_data_input);

    for (unsigned int j=0; j<2; j++) {
        eq_data_->reallocate_cache();

        for(unsigned int i=0; i < mesh_->n_elements(); i++) {  // time loop
            EXPECT_TRUE( eval_boundary_field(eq_data_->bc_scalar_field, expected_scalars[j][i], i, i) );
            EXPECT_TRUE( eval_boundary_field(eq_data_->bc_vector_field, expected_vectors[j][i], i, i) );
        }
        eq_data_->tg_.next_time();
    }
}


TEST_F(FieldEvalFETest, native) {
    string eq_data_input = R"YAML(
    data:
      - region: ALL
        time: 0.0
        scalar_field: !FieldFE
          mesh_data_file: output/test_output_vtk_ascii_ref.vtu
          field_name: flow_data
    )YAML";

    std::vector<double> expected_scalars = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };

    this->create_mesh("fields/simplest_cube_3d.msh");
    this->read_input(eq_data_input);

    eq_data_->reallocate_cache();
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, expected_scalars) );
}


TEST_F(FieldEvalFETest, set_fe_data_scalar) {
    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarFieldFE;
    this->create_mesh("fields/one_element_2d.msh");
    this->set_dof_values( {0.5} );

    MixedPtr<FE_P_disc> fe(0);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe);
    dh_->distribute_dofs(ds);

    std::shared_ptr<ScalarFieldFE> fe_field = std::make_shared<ScalarFieldFE>();
    fe_field->set_fe_data(dh_, v);
    eq_data_->scalar_field.set(fe_field, 0.0);

    eq_data_->reallocate_cache();
    EXPECT_TRUE( eval_bulk_field(eq_data_->scalar_field, 0.5) );
}


TEST_F(FieldEvalFETest, set_fe_data_vector) {
    typedef FieldFE<3, FieldValue<3>::VectorFixed > VectorFieldFE;
    this->create_mesh("fields/one_element_2d.msh");
    this->set_dof_values( {0.5, 1.5, 2.5} );

    MixedPtr<FE_RT0> fe;
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>(mesh_, fe);
    dh_->distribute_dofs(ds);

    std::shared_ptr<VectorFieldFE> fe_field = std::make_shared<VectorFieldFE>();
    fe_field->set_fe_data(dh_, v);
    eq_data_->vector_field.set(fe_field, 0.0);

    eq_data_->reallocate_cache();
    arma::vec3 expected = { 1./7, 2./7, 0.0 };
    EXPECT_TRUE( eval_bulk_field(eq_data_->vector_field, expected) );
}


/*******************************************************************************
 *                                                                             *
 *     New tests of Elementwise and Interpolation P0 replaced with FieldFE     *
 *     TODO: Remove - replaced by FieldEvalFETest                              *
 *                                                                             *
 *******************************************************************************/
//
//string elem_input = R"YAML(
//##### tests of elementwise
//scalar: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: scalar
//  default_value: 0.0
//scalar_unit_conversion: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: scalar
//  unit: "const; const=100*m^0"
//  default_value: 0.0
//scalar_time_shift: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: scalar
//  default_value: 0.0
//  read_time_shift: 1.0
//enum: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: enum
//  default_value: 0
//vector_fixed: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: vector_fixed
//  default_value: 0.0
//tensor_fixed: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: tensor_fixed
//  default_value: 0.0
//vtk_scalar: !FieldFE
//  mesh_data_file: fields/vtk_ascii_data.vtu
//  field_name: scalar_field
//vtk_vector: !FieldFE
//  mesh_data_file: fields/vtk_binary_data.vtu
//  field_name: vector_field
//vtk_tensor: !FieldFE
//  mesh_data_file: fields/vtk_compressed_data.vtu
//  field_name: tensor_field
//default_values: !FieldFE
//  mesh_data_file: fields/simplest_cube_data.msh
//  field_name: porosity
//  default_value: 0.1
//scalar_identic_mesh: !FieldFE
//  mesh_data_file: fields/identic_mesh_data.msh
//  field_name: scalar
//  default_value: 0.0
//  interpolation: identic_mesh
//vector_identic_mesh: !FieldFE
//  mesh_data_file: fields/identic_mesh_data.msh
//  field_name: vector_fixed
//  default_value: 0.0
//  interpolation: identic_mesh
//##### tests of intersection interpolation P0
//interp_scalar_intersect: !FieldFE
//  mesh_data_file: fields/interpolate_boundary_data.msh
//  field_name: scalar
//  default_value: 0.0
//  interpolation: P0_intersection
//interp_vector_fixed_intersect: !FieldFE
//  mesh_data_file: fields/interpolate_boundary_data.msh
//  field_name: vector_fixed
//  default_value: 0.0
//  interpolation: P0_intersection
//##### tests of gauss interpolation P0
//interp_scalar_gauss: !FieldFE
//  mesh_data_file: fields/interpolate_boundary_data.msh
//  field_name: scalar
//  default_value: 0.0
//interp_vector_fixed_gauss: !FieldFE
//  mesh_data_file: fields/interpolate_boundary_data.msh
//  field_name: vector_fixed
//  default_value: 0.0
//#interp_scalar_gauss_large: !FieldFE
//#  mesh_data_file: fields/bigger_3d_cube_0.5.msh
//#  field_name: scalar
//#interp_tensor_gauss_fixed: !FieldFE
//#  mesh_data_file: fields/interpolate_boundary_data.msh
//#  field_name: tensor_fixed
//)YAML";
//
//
//class FieldFENewTest : public testing::Test {
//public:
//    typedef FieldFE<3, FieldValue<3>::Scalar > ScalarField;
//    typedef FieldFE<3, FieldValue<3>::Enum > EnumField;
//    typedef FieldFE<3, FieldValue<3>::VectorFixed > VecFixField;
//    typedef FieldFE<3, FieldValue<3>::TensorFixed > TensorField;
//
//    virtual void SetUp() {
//        // setup FilePath directories
//        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
//
//        Profiler::instance();
//        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
//
//        mesh = mesh_full_constructor("{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }");
//
//        Input::Type::Record rec_type = Input::Type::Record("Test","")
//            .declare_key("scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("scalar_unit_conversion", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("scalar_time_shift", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("enum", EnumField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("vector_fixed", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("tensor_fixed", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("vtk_scalar", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("vtk_vector", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("vtk_tensor", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("default_values", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("scalar_identic_mesh", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("vector_identic_mesh", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("interp_scalar_intersect", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("interp_vector_fixed_intersect", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("interp_scalar_gauss", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .declare_key("interp_vector_fixed_gauss", VecFixField::get_input_type(), Input::Type::Default::obligatory(),"" )
////            .declare_key("interp_scalar_gauss_large", ScalarField::get_input_type(), Input::Type::Default::obligatory(),"" )
////            .declare_key("interp_tensor_gauss_fixed", TensorField::get_input_type(), Input::Type::Default::obligatory(),"" )
//            .close();
//
//        Input::ReaderToStorage reader( elem_input, rec_type, Input::FileFormat::format_YAML );
//        rec=reader.get_root_interface<Input::Record>();
//
//        test_time[0] = 0.0;
//        test_time[1] = 1.0;
//        test_time[2] = 2.0;
//
//    }
//    virtual void TearDown() {
//    	delete mesh;
//    }
//
//    const FieldAlgoBaseInitData& init_data(std::string field_name) {
//    	static const FieldAlgoBaseInitData init_data(field_name, 0, UnitSI::dimensionless());
//    	return init_data;
//    }
//
//    Mesh * mesh;
//    Input::Record rec;
//    Space<3>::Point point;
//    double test_time[3];
//
//};
//
//
//
//TEST_F(FieldFENewTest, scalar) { // moved to FieldEvalFETest, input_msh
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
//    field.set_mesh(mesh,false);
//
//    for (unsigned int j=0; j<2; j++) {
//        field.set_time(test_time[j]);
//        for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_DOUBLE_EQ( j*0.1+(i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, bc_scalar) { // moved to FieldEvalFETest, input_msh
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar"), init_data("scalar"));
//    BCMesh *bc_mesh = mesh->bc_mesh();
//    field.set_mesh(mesh,true);
//
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//
//        // only 4 BC elements are compatible with the comp mesh
//        for(unsigned int i=0; i < 4; i++) {
//            auto ele = bc_mesh->element_accessor(i);
//            EXPECT_DOUBLE_EQ( j*0.1 + (ele.input_id()+1)*0.1 , field.value(point,ele) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, scalar_unit_conv) { // moved to FieldEvalFETest, unit_conversion
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar_unit_conversion"), init_data("scalar_unit_conversion"));
//    field.set_mesh(mesh,false);
//
//    for (unsigned int j=0; j<2; j++) {
//        field.set_time(test_time[j]);
//        for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_DOUBLE_EQ( j*10.0+(i+1)*10.0 , field.value(point,mesh->element_accessor(i)) );
//        }
//    }
//    field.set_time(test_time[2]); //temporary solution, remove this line after fix 'bc_scalar_unit_conv' test
//}
//
//
//TEST_F(FieldFENewTest, bc_scalar_unit_conv) { // moved to FieldEvalFETest, unit_conversion
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar_unit_conversion"), init_data("scalar_unit_conversion"));
//    field.set_mesh(mesh,true);
//    for (unsigned int j=0; j<3; j++) {
//    	field.set_time(test_time[j]);
//
//        BCMesh *bc_mesh = mesh->bc_mesh();
//        // only 4 BC elements are compatible with the comp mesh
//        for(unsigned int i=0; i < 4; i++) {
//            auto ele = bc_mesh->element_accessor(i);
//            EXPECT_DOUBLE_EQ(j*10.0+(ele.input_id()+1)*10.0 , field.value(point,ele) );
//        }
//    }
//
//}
//
//
//TEST_F(FieldFENewTest, scalar_time_shift) { // moved to FieldEvalFETest, time_shift
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar_time_shift"), init_data("scalar_time_shift"));
//    field.set_mesh(mesh,false);
//
//    for (unsigned int j=0; j<2; j++) {
//        field.set_time(test_time[j]);
//        for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_DOUBLE_EQ( j*0.1+(i+2)*0.1 , field.value(point,mesh->element_accessor(i)) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, vector_fixed) { // moved to FieldEvalFETest, input_msh
//	string expected_vals[2] = {"1 2 3", "2 3 4"};
//    VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data("vector_fixed"));
//    field.set_mesh(mesh,false);
//     for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//         for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, bc_vector_fixed) { // moved to FieldEvalFETest, input_msh
//	string expected_vals[2] = {"4 5 6", "5 6 7"};
//    VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("vector_fixed"), init_data("vector_fixed"));
//    field.set_mesh(mesh,true);
//    BCMesh *bc_mesh = mesh->bc_mesh();
//     for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//        // only 6 BC elements are compatible with the comp mesh
//        for(unsigned int i=0; i < 4; i++) {
//            auto ele = bc_mesh->element_accessor(i);
//            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,ele)) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, tensor_fixed) { // moved to FieldEvalFETest, input_msh
//	string expected_vals[2] = {"1 2 3; 4 5 6; 7 8 9", "2 3 4; 5 6 7; 8 9 10"};
//    TensorField field;
//    field.init_from_input(rec.val<Input::Record>("tensor_fixed"), init_data("tensor_fixed"));
//    field.set_mesh(mesh,false);
//     for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//     	for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,mesh->element_accessor(i)) );
//            EXPECT_TRUE( match.min() );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, bc_tensor_fixed) { // moved to FieldEvalFETest, input_msh
//	string expected_vals[2] = {"4 5 6; 7 8 9; 10 11 12", "5 6 7; 8 9 10; 11 12 13"};
//    TensorField field;
//    field.init_from_input(rec.val<Input::Record>("tensor_fixed"), init_data("tensor_fixed"));
//    field.set_mesh(mesh, true);
//    BCMesh *bc_mesh = mesh->bc_mesh();
//     for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//        // only 4 BC elements are compatible with the comp mesh
//        for(unsigned int i=0; i < 4; i++) {
//            auto ele = bc_mesh->element_accessor(i);
//            arma::umat match = ( arma::mat33(expected_vals[j]) == field.value(point,ele) );
//            EXPECT_TRUE( match.min() );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, vtk_scalar) { // moved to FieldEvalFETest, input_vtk
//	ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("vtk_scalar"), init_data("vtk_scalar"));
//    field.set_mesh(mesh, false);
//	field.set_time(0.0);
//
//	for(unsigned int i=0; i<mesh->n_elements(); i++) {
//		EXPECT_DOUBLE_EQ( (i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
//	}
//}
//
//
//
//TEST_F(FieldFENewTest, vtk_vector) { // moved to FieldEvalFETest, input_vtk
//	string expected_vals = "0.5 1 1.5";
//    VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("vtk_vector"), init_data("vtk_vector"));
//    field.set_mesh(mesh, false);
//   	field.set_time(0.0);
//    for(unsigned int i=0; i < mesh->n_elements(); i++) {
//    	EXPECT_TRUE( arma::min(arma::vec3(expected_vals) == field.value(point,mesh->element_accessor(i))) );
//    }
//}
//
//
//TEST_F(FieldFENewTest, vtk_tensor) { // moved to FieldEvalFETest, input_vtk
//	string expected_vals = "1 2 3; 4 5 6; 7 8 9";
//    TensorField field;
//    field.init_from_input(rec.val<Input::Record>("vtk_tensor"), init_data("vtk_tensor"));
//    field.set_mesh(mesh,false);
//   	field.set_time(0.0);
//    for(unsigned int i=0; i < mesh->n_elements(); i++) {
//    	arma::umat match = ( arma::mat33(expected_vals) == field.value(point,mesh->element_accessor(i)) );
//        EXPECT_TRUE( match.min() );
//    }
//}
//
//
//TEST_F(FieldFENewTest, scalar_enum) { // moved to FieldEvalFETest, input_msh
//    EnumField field;
//    field.init_from_input(rec.val<Input::Record>("enum"), init_data("enum"));
//    field.set_mesh(mesh,false);
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//     	for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_EQ( j, field.value(point,mesh->element_accessor(i)) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, bc_scalar_enum) { // moved to FieldEvalFETest, input_msh
//    EnumField field;
//    field.init_from_input(rec.val<Input::Record>("enum"), init_data("enum"));
//    field.set_mesh(mesh, true);
//    for (unsigned int j=0; j<2; j++) {
//		field.set_time(test_time[j]);
// 		for(unsigned int i=0; i < 4; i++) {
//			EXPECT_EQ( j+1, field.value(point,mesh->bc_mesh()->element_accessor(i)) );
//		}
//    }
//}
//
//
//TEST_F(FieldFENewTest, default_values) { // moved to FieldEvalFETest, default_values
//	string expected_vals = "0.1 0.1 0.1";
//    VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("default_values"), init_data("default_values"));
//    field.set_mesh(mesh,true);
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//     	for(unsigned int i=0; i < 4; i++) {
//            EXPECT_TRUE( arma::min(arma::vec3(expected_vals) == field.value(point,mesh->bc_mesh()->element_accessor(i))) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, scalar_identic_mesh) { // moved to FieldEvalFETest, identic_mesh
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar_identic_mesh"), init_data("scalar_identic_mesh"));
//    field.set_mesh(mesh,false);
//
//    for (unsigned int j=0; j<2; j++) {
//        field.set_time(test_time[j]);
//        for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_DOUBLE_EQ( 1.0+j*0.1+(i+1)*0.1 , field.value(point,mesh->element_accessor(i)) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, bc_scalar_identic_mesh) { // moved to FieldEvalFETest, identic_mesh
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("scalar_identic_mesh"), init_data("scalar_identic_mesh"));
//    field.set_mesh(mesh,true);
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//
//        for(unsigned int i=9; i < 13; i++) {
//            EXPECT_DOUBLE_EQ( 2.0+j*0.1+(i-8)*0.1 , field.value(point,mesh->bc_mesh()->element_accessor(i-9)) );
//        }
//    }
//
//}
//
//
//TEST_F(FieldFENewTest, vector_fixed_identic_mesh) { // moved to FieldEvalFETest, identic_mesh
//	string expected_vals[2] = {"3 4 5", "6 7 8"};
//    VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("vector_identic_mesh"), init_data("vector_identic_mesh"));
//    field.set_mesh(mesh,false);
//     for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//         for(unsigned int i=0; i < mesh->n_elements(); i++) {
//            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->element_accessor(i))) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, bc_vector_fixed_identic_mesh) { // moved to FieldEvalFETest, identic_mesh
//	string expected_vals[2] = {"1 2 3", "4 5 6"};
//    VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("vector_identic_mesh"), init_data("vector_identic_mesh"));
//    field.set_mesh(mesh,true);
//     for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//     	for(unsigned int i=0; i < 4; i++) {
//            EXPECT_TRUE( arma::min(arma::vec3(expected_vals[j]) == field.value(point,mesh->bc_mesh()->element_accessor(i))) );
//        }
//    }
//}
//
//
//TEST_F(FieldFENewTest, intersection_1d_2d_elements_small_scalar) {
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("interp_scalar_intersect"), init_data("interp_scalar_intersect"));
//    field.set_mesh(mesh, true);
//    //std::vector<unsigned int> expected_vals = {4,9,6,7};
//
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//    	std::cout << "Time: " << test_time[j] << std::endl;
//
//    	for (unsigned int i=0; i<4; ++i) {
//    		ElementAccessor<3> elm = mesh->bc_mesh()->element_accessor(i);
//    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
//    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
//    	}
//    }
//
//}
//
//
//TEST_F(FieldFENewTest, intersection_1d_2d_elements_small_vector) {
//	VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("interp_vector_fixed_intersect"), init_data("interp_vector_fixed_intersect"));
//    field.set_mesh(mesh, true);
//    //std::vector<unsigned int> expected_vals = {4,9,6,7};
//
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//    	std::cout << "Time: " << test_time[j] << std::endl;
//
//    	for (unsigned int i=0; i<4; ++i) {
//    		ElementAccessor<3> elm = mesh->bc_mesh()->element_accessor(i);
//    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
//    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
//    	}
//    }
//
//}
//
//
//TEST_F(FieldFENewTest, gauss_1d_2d_elements_small_scalar) { // moved to FieldEvalFETest, interpolation_gauss
//    ScalarField field;
//    field.init_from_input(rec.val<Input::Record>("interp_scalar_gauss"), init_data("interp_scalar_gauss"));
//    field.set_mesh(mesh, true);
//    std::vector<unsigned int> expected_vals = {4,9,6,7};
//
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//    	std::cout << "Time: " << test_time[j] << std::endl;
//
//    	for (unsigned int i=0; i<4; ++i) {
//    		ElementAccessor<3> elm = mesh->bc_mesh()->element_accessor(i);
//    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
//    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
//    	}
//    }
//
//}
//
//
//TEST_F(FieldFENewTest, gauss_1d_2d_elements_small_vector) { // moved to FieldEvalFETest, interpolation_gauss
//	VecFixField field;
//    field.init_from_input(rec.val<Input::Record>("interp_vector_fixed_gauss"), init_data("interp_vector_fixed_gauss"));
//    field.set_mesh(mesh, true);
//    std::vector<unsigned int> expected_vals = {4,9,6,7};
//
//    for (unsigned int j=0; j<2; j++) {
//    	field.set_time(test_time[j]);
//    	std::cout << "Time: " << test_time[j] << std::endl;
//
//    	for (unsigned int i=0; i<4; ++i) {
//    		ElementAccessor<3> elm = mesh->bc_mesh()->element_accessor(i);
//    		std::cout << " - " << field.value(elm.centre(), elm) << std::endl;
//    		//EXPECT_DOUBLE_EQ( 0.1*(j+expected_vals[i]), field.value(point, mesh->element_accessor(i+9)) );
//    	}
//    }
//
//}
