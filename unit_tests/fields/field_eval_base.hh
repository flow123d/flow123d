/*!
 *
ï»¿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 *
 * @file    field_eval_base.hh
 * @brief
 */

#ifndef FIELD_EVAL_BASE_HH_
#define FIELD_EVAL_BASE_HH_

#define TEST_USE_MPI
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>

#include "fem/eval_points.hh"
#include "fem/eval_subset.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_values.hh"
#include "fields/field_set.hh"
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
#include "la/vector_mpi.hh"


/**
 * Main class of field evaluation unit tests.
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
class FieldEvalBaseTest : public testing::Test {
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

    template<class RefVal>
    class VecRef {
    public:
        VecRef(std::vector<RefVal> vector) : _vec(vector) {};
        RefVal value(int i, BulkPoint p) {
            return _vec[i];
        }

        std::vector<RefVal> _vec;
    };

    template<class RefVal>
    class SingleValRef {
    public:
    	SingleValRef(RefVal val) : _val(val) {};
    	RefVal value(int i, BulkPoint p) {
            return _val;
        }

    	RefVal _val;
    };

    template<class RefField>
    class FieldRef {
    public:
        FieldRef(RefField &field) : _field(field) {};
        typename RefField::ValueType::return_type value(int i, BulkPoint p) {
            return _field(p);
        }

        RefField &_field;
    };

    FieldEvalBaseTest() {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        eq_data_ = std::make_shared<EqData>();
    }

    ~FieldEvalBaseTest() {
        Profiler::uninitialize();
    }

    static Input::Type::Record & get_input_type() {
        return IT::Record("SomeEquation","")
                .declare_key("data", IT::Array(
                        IT::Record("SomeEquation_Data", FieldCommon::field_descriptor_record_description("SomeEquation_Data") )
                        .copy_keys( FieldEvalBaseTest::EqData().make_field_descriptor_type("SomeEquation") )
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
    template<class FieldType, class RefType>
    bool eval_bulk_field(FieldType &eval_field, RefType &ref_obj)
    {
        bool is_passed = true;
	    for(unsigned int i=0; i < mesh_->n_elements(); i++) {
            eq_data_->computed_dh_cell_ = DHCellAccessor(dh_.get(), i);
            uint dim = eq_data_->computed_dh_cell_.dim()-1;
            eq_data_->update_cache();

            auto p = *( eq_data_->mass_integral[dim]->points(eq_data_->position_in_cache(eq_data_->computed_dh_cell_.elm_idx()), eq_data_.get()).begin() );
            is_passed = is_passed & check_point_value(eval_field, p, ref_obj.value(i, p) );
        }
	    return is_passed;
    }


    /// Evaluates and compare boundary field, uses FieldFormula as reference
    template<class FieldType, class RefType>
    bool eval_boundary_field(FieldType &eval_field, RefType &ref_obj, unsigned int i_bulk_elem, unsigned int i_bdr_elem)
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
                    return check_point_value(eval_field, p_bdr, ref_obj.value(0, p_bdr)); // boundary compares only single value, we don't need 'i'
                }
            }
        }
        // Boundary not found - test doesn' pass
        std::cout << "Wrong indices pair of bulk and boundary element!" << std::endl;
        std::cout << "Boundary element of idx: " << i_bdr_elem << " is not relevant to bulk element of idx: " << i_bulk_elem << std::endl;
	    return false;
    }


    std::shared_ptr<EqData> eq_data_;
    Mesh * mesh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_;
};


/// General function compares value of Vector or Tensor fields with ref_val
template<class Val>
bool FieldEvalBaseTest::compare_vals(const Val &field_val, const Val &ref_val)
{
    return _expect_arma_eqal(field_val, ref_val, std::cout);
}


/// Template specialization for Enum fields.
template<>
bool FieldEvalBaseTest::compare_vals(const unsigned int &field_val, const unsigned int &ref_val)
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
bool FieldEvalBaseTest::compare_vals(const double &field_val, const double &ref_val)
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

#endif /* FIELD_EVAL_BASE_HH_ */
