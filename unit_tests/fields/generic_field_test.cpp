
#define TEST_USE_PETSC
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>


#include "fields/generic_field.hh"
#include "fem/eval_points.hh"
#include "fem/integral_acc.hh"
#include "fem/element_cache_map.hh"
#include "fields/field_set.hh"
#include "mesh/mesh.h"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"
#include "io/msh_gmshreader.h"
#include "system/sys_profiler.hh"
#include "fields/field_flag.hh"
#include "tools/unit_si.hh"
#include "quadrature/quadrature.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/dh_cell_accessor.hh"
#include "system/sys_profiler.hh"
#include "coupling/generic_assembly.hh"
#include "coupling/assembly_base.hh"


class GenericFieldTest : public testing::Test {

public:
    class EqOutput : public FieldSet, public ElementCacheMap {
    public:
        EqOutput() {
            *this += region_id.name("region_id")
            	        .units( UnitSI::dimensionless())
            	        .flags(FieldFlag::input_copy);

            //*this += subdomain.name("subdomain")
            //  .units( UnitSI::dimensionless() )
            //  .flags(FieldFlag::input_copy);


            // Asumme following types:
            eval_points_ = std::make_shared<EvalPoints>();
            Quadrature *q_bulk_1d = new QGauss(1, 0);
            Quadrature *q_bulk_2d = new QGauss(2, 0);
            Quadrature *q_bulk_3d = new QGauss(3, 0);
            bulk_int[0] = std::make_shared<BulkIntegral>(eval_points_, q_bulk_1d, 1);
            bulk_int[1] = std::make_shared<BulkIntegral>(eval_points_, q_bulk_2d, 2);
            bulk_int[2] = std::make_shared<BulkIntegral>(eval_points_, q_bulk_3d, 3);
            this->init(eval_points_);
        }

        void register_eval_points() {
            for(auto dh_cell : dh_->local_range() ) {
                uint subset_idx = bulk_int[dh_cell.dim()-1]->get_subset_idx();
                unsigned int reg_idx = dh_cell.elm().region_idx().idx();
                for (uint i=uint( eval_points_->subset_begin(dh_cell.dim(), subset_idx) );
                          i<uint( eval_points_->subset_end(dh_cell.dim(), subset_idx) ); ++i) {
                    this->add_eval_point(reg_idx, dh_cell.elm_idx(), i, dh_cell.local_idx());
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


        // fields
        Field<3, FieldValue<3>::Scalar> region_id;
        //Field<3, FieldValue<3>::Scalar> subdomain; // test of subdomain is not solved now
        std::shared_ptr<EvalPoints> eval_points_;
        std::array<std::shared_ptr<BulkIntegral>, 3> bulk_int;  // dim 1,2,3
        std::shared_ptr<DOFHandlerMultiDim> dh_;
    };

    GenericFieldTest()
    {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        Profiler::instance();
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        eq_output_ = std::make_shared<EqOutput>();
        eq_output_->add_coords_field();
        mesh_ = mesh_full_constructor("{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }");
        eq_output_->dh_ = std::make_shared<DOFHandlerMultiDim>(*mesh_);
    }

    ~GenericFieldTest() {
        Profiler::uninitialize();
    }

    void initialize() {
        TimeGovernor tg(0.0, 0.5);
        eq_output_->set_mesh(*mesh_);
        eq_output_->region_id = GenericField<3>::region_id(*mesh_);
        //eq_output_->subdomain = GenericField<3>::subdomain(*mesh_);
        eq_output_->set_time(tg.step(), LimitSide::right);
        eq_output_->cache_reallocate( *(eq_output_.get()), *(eq_output_.get()) );
    }


    std::shared_ptr<EqOutput> eq_output_;
    Mesh * mesh_;
};


TEST_F(GenericFieldTest, all) {
    this->initialize();
    eq_output_->start_elements_update();
    eq_output_->update_cache();

    for(auto dh_cell : eq_output_->dh_->local_range() ) {
        auto p = *( eq_output_->bulk_int[dh_cell.dim()-1]->points(eq_output_->position_in_cache(dh_cell.elm_idx()), eq_output_.get()).begin() );
        EXPECT_EQ( dh_cell.elm().region().id(), eq_output_->region_id(p) );
    }
}
