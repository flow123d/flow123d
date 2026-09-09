/*
 * fe_values_test.cpp
 *
 *  Created on: Sep 9, 2012
 *      Author: jb
 */

#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>

#include <cmath>
#include <mesh_constructor.hh>
#include "arma_expect.hh"
#include "armadillo"
#include "system/armadillo_tools.hh"
#include "system/sys_profiler.hh"
#include "quadrature/quadrature_lib.hh"
#include "fem/fe_p.hh"
#include "fem/fe_rt.hh"
#include "fem/integral_acc.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fem/patch_fe_values.hh"
#include "fem/patch_op_impl.hh"
#include "tools/revertable_list.hh"
#include "fem/mapping_p1.hh"
#include "mesh/mesh.h"
#include "mesh/elements.h"
#include "mesh/region.hh"
#include "mesh/accessors.hh"
#include "fem/fe_system.hh"



/*NodeVector make_nodes(const std::vector<string> &nodes_str)
{
  std::vector<arma::vec3> nodes;
  for(auto str : nodes_str) nodes.push_back( arma::vec3(str));
  
  NodeVector node_vector(nodes.size());
  unsigned int i=0;
  for (auto node : nodes)
  {
    node_vector.add_item(i);
    node_vector[i++] = Node(node[0], node[1], node[2]);
  }
  
  return node_vector;
}


vector<Element> make_elements(NodeVector &node_vector, const std::vector<std::vector<unsigned int> > &node_idx)
{
  vector<Element> el_vec(node_idx.size());
  
  unsigned int iel = 0;
  for (auto nodes : node_idx)
  {
    el_vec[iel].init(nodes.size()-1, iel, nullptr, RegionIdx());
    unsigned int i=0;
    for(auto node : nodes)
      el_vec[iel].node[i++] = &node_vector[node];
    
    iel++;
  }
  
  return el_vec;
}*/


// Define EXPECT_<...>_NEAR with fixed abs_error 1e-10
#define EXPECT_TEST_NEAR( A, B )\
  EXPECT_NEAR(A, B, 1e-9)

#define EXPECT_TEST_ARMA_NEAR( A, B ) \
    EXPECT_ARMA_NEAR(A, B, 1e-9);



/**
 * Base test class - simulation of equation and assembly without fields
 */
class PatchFETestBase {
public:

    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmBase {
    public:
        /**
         * Constructor
         *
         * @param quad_order      Order of Quadrature (quad_, quad_low_) objects.
         * @param quad_diff_order Order of Quadrature (quad_diff_order_, quad_low_diff_order_) objects.
         */
        AsmBase(PatchFETestBase *generic, uint quad_order, uint quad_order_diff = -1)
        : generic_(generic),
          quad_( new QGauss(dim, 2*quad_order) ),
          quad_diff_( new QGauss(dim, 2*( (quad_order_diff==-1) ? quad_order : quad_order_diff ) ) ),
          bulk_integral_( create_bulk_integral(quad_) ),
          bulk_integral_diff_( create_bulk_integral(quad_diff_) )
        {}

    	/// Destructor
        virtual ~AsmBase() {
            delete quad_;
            delete quad_diff_;
        }


        std::shared_ptr<BulkIntegralAcc<dim>> create_bulk_integral(Quadrature *quad) {
            ASSERT_PERMANENT_EQ(quad->dim(), dim);
            std::tuple<uint, uint> tpl = IntegralTplHash::integral_tuple(dim, quad->size());
            auto result = integrals_.bulk_.insert({
                    tpl,
                    std::make_shared<BulkIntegralAcc<dim>>(generic_->patch_internals_, quad)
                });
            return result.first->second;
        }

        /** Declaration of data members **/
        PatchFETestBase *generic_;                                        ///< pointer to generic object
        Quadrature *quad_;                                                ///< Quadrature (of dim)
        Quadrature *quad_diff_;                                           ///< Quadrature (of dim) used only in mixed FESystem test
        DimIntegrals<dim> integrals_;                                     ///< Set of used integrals.
        std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;             ///< BulkIntegral
        std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_diff_;        ///< BulkIntegral used only in mixed FESystem test
    };


    PatchFETestBase(std::shared_ptr<DOFHandlerMultiDim> dh)
    : dh_(dh),
	  patch_internals_(dh_->ds()->fe())
    {
        used_element_idx_ = {0};

        expected_vector_shape_ = {
            {
                {0.585410196624968, 0, 0}, {0.138196601125011, 0, 0}, {0.138196601125011, 0, 0}, {0.138196601125011, 0, 0},
                {0, 0.585410196624968, 0}, {0, 0.138196601125011, 0}, {0, 0.138196601125011, 0}, {0, 0.138196601125011, 0},
                {0, 0, 0.585410196624968}, {0, 0, 0.138196601125011}, {0, 0, 0.138196601125011}, {0, 0, 0.138196601125011}
            },
    		{
    		    {0.138196601125011, 0, 0}, {0.138196601125011, 0, 0}, {0.138196601125011, 0, 0}, {0.585410196624968, 0, 0},
    		    {0, 0.138196601125011, 0}, {0, 0.138196601125011, 0}, {0, 0.138196601125011,0 }, {0, 0.585410196624968, 0},
    		    {0, 0, 0.138196601125011}, {0, 0, 0.138196601125011}, {0, 0, 0.138196601125011}, {0, 0, 0.585410196624968}
    		},
    		{
    		    {0.138196601125011, 0, 0}, {0.138196601125011, 0, 0}, {0.585410196624968, 0, 0}, {0.138196601125011, 0, 0},
    		    {0, 0.138196601125011,0 }, {0, 0.138196601125011,0 }, {0, 0.585410196624968, 0}, {0, 0.138196601125011, 0},
    		    {0, 0, 0.138196601125011}, {0, 0, 0.138196601125011}, {0, 0, 0.585410196624968}, {0, 0, 0.138196601125011}
    		},
    		{
    		    {0.138196601125011, 0, 0}, {0.585410196624968, 0, 0}, {0.138196601125011, 0, 0}, {0.138196601125011, 0, 0},
                {0, 0.138196601125011, 0}, {0, 0.585410196624968, 0}, {0, 0.138196601125011, 0}, {0, 0.138196601125011, 0},
                {0, 0, 0.138196601125011}, {0, 0, 0.585410196624968}, {0, 0, 0.138196601125011}, {0, 0, 0.138196601125011}
    		}
        };
        expected_tensor_mats_ = {
            { 1, 0, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 1, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 1, 0, 0 },
            { 0, 1, 0, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 1, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 1, 0 },
            { 0, 0, 1, 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 1, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0, 0, 0, 1 }
        };
        expected_tensor_coefs_ = {
            {0.585410196624968, 0.138196601125011, 0.138196601125011, 0.138196601125011},
            {0.138196601125011, 0.138196601125011, 0.138196601125011, 0.585410196624968},
            {0.138196601125011, 0.138196601125011, 0.585410196624968, 0.138196601125011},
            {0.138196601125011, 0.585410196624968, 0.138196601125011, 0.138196601125011}
        };
        expected_rt_shape_ = {
            { {1.170820393249937, 0.276393202250021, 0.276393202250021}, {1.17082039324994, 0.276393202250021, -1.72360679774998},
              {1.170820393249937, -1.72360679774998, 0.276393202250021}, {-0.829179606750063, 0.276393202250021, 0.276393202250021} },
            { {0.276393202250021, 0.276393202250021, 0.276393202250021}, {0.276393202250021, 0.276393202250021, -1.72360679774998},
              {0.276393202250021, -1.72360679774998, 0.276393202250021}, {-1.72360679774998, 0.276393202250021, 0.276393202250021} },
            { {0.276393202250021, 0.276393202250021, 1.170820393249937}, {0.276393202250021, 0.276393202250021, -0.829179606750063},
              {0.276393202250021, -1.72360679774998, 1.170820393249937}, {-1.72360679774998, 0.276393202250021, 1.17082039324994} },
            { {0.276393202250021, 1.170820393249937, 0.276393202250021}, {0.276393202250021, 1.170820393249937, -1.72360679774998},
              {0.276393202250021, -0.829179606750063, 0.276393202250021}, {-1.72360679774998, 1.170820393249937, 0.276393202250021} }
    	};

    }

    ~PatchFETestBase() {}

    void add_bulk_integral(DHCellAccessor cell, std::shared_ptr<BulkIntegral> bulk_integral) {
        uint subset_idx = bulk_integral->get_subset_idx();
        bulk_integral->patch_data().emplace_back(cell);
        uint dim = cell.dim();

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( patch_internals_.eval_points_->subset_begin(dim, subset_idx) );
                  i<uint( patch_internals_.eval_points_->subset_end(dim, subset_idx) ); ++i) {
            patch_internals_.element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }
    }

    void initialize(MixedPtr<AsmBase, 1> multidim_asm) {
        this->bulk_integrals_[0] = multidim_asm[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm[3_d]->bulk_integral_;
        this->bulk_integrals_diff_[0] = multidim_asm[1_d]->bulk_integral_diff_;
        this->bulk_integrals_diff_[1] = multidim_asm[2_d]->bulk_integral_diff_;
        this->bulk_integrals_diff_[2] = multidim_asm[3_d]->bulk_integral_diff_;

        this->patch_internals_.fe_values_.init_finalize();
    }

    virtual void test_evaluation() =0;

    void create_patch(MixedPtr<AsmBase, 1> multidim_asm) {
        //for (auto elm_idx : used_element_idx_)
        {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element( used_element_idx_[0] );
            auto &ppv_bulk = patch_internals_.fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
            this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
            this->patch_internals_.fe_values_.make_permanent_ppv_data();
        }
        multidim_asm[1_d]->integrals_.make_permanent();
        multidim_asm[2_d]->integrals_.make_permanent();
        multidim_asm[3_d]->integrals_.make_permanent();

        patch_internals_.element_cache_map_.make_paermanent_eval_points();
        patch_internals_.element_cache_map_.create_patch();
    }

    void update_patch(MixedPtr<AsmBase, 1> multidim_asm) {
    	patch_internals_.fe_values_.prepare_new_patch(this->patch_internals_.eval_points_);
    	patch_internals_.fe_values_.add_patch_points<3>(multidim_asm[3_d]->integrals_, &this->patch_internals_.element_cache_map_);
    	patch_internals_.fe_values_.add_patch_points<2>(multidim_asm[2_d]->integrals_, &this->patch_internals_.element_cache_map_);
    	patch_internals_.fe_values_.add_patch_points<1>(multidim_asm[1_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.reinit_patch();
    }


    std::shared_ptr<DOFHandlerMultiDim> dh_;

    PatchInternals patch_internals_;                                          ///< Holds common patch objects (EvalPoints, ElementCacheMap ...)
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_;             ///< Bulk integrals of dim 1,2,3
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_diff_;        ///< Bulk integrals of dim 1,2,3
    std::vector<unsigned int> used_element_idx_;                              ///< List of mesh idx of elements used in tests

    // Reference values
    std::vector< std::vector<arma::vec3> > expected_vector_shape_;
    std::vector<arma::mat33> expected_tensor_mats_;
    std::vector< std::vector<double> > expected_tensor_coefs_;
    std::vector< std::vector< arma::vec3 > > expected_rt_shape_;
};


/**
 * Specialization of vector FeSystem
 */
class PatchFETestVector : public PatchFETestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmVector : public PatchFETestBase::AsmBase<dim> {
    public:
        /// Constructor
        AsmVector(PatchFETestBase *generic, uint quad_order)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order),
//          generic_inst_(generic),
          vector_shape_( this->bulk_integral_->vector_shape() ),
          grad_vector_shape_( this->bulk_integral_->grad_vector_shape() )
        {}

        /// Destructor
        virtual ~AsmVector() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int n_dofs) {
        	unsigned int k=0;
            for ( auto p : this->bulk_integral_->points( this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx()) ) ) {
                for (unsigned int i=0; i<n_dofs; i++) {
                    // check values
                    EXPECT_TEST_ARMA_NEAR( this->generic_->expected_vector_shape_[k][i], vector_shape_.shape(i)(p) );
                    for (unsigned int c=0; c<3; c++)
                    {
                        //check gradients
                        arma::rowvec gr = grad_vector_shape_.shape(i)(p).row(c);
                        if (i / 4 == c) { // gradient of nonzero component
                            switch (i%4) {
                            case 0:
                                EXPECT_ARMA_EQ( arma::rowvec("1 0 0"), gr );
                                break;
                            case 1:
                                EXPECT_ARMA_EQ( arma::rowvec("0 1 0"), gr );
                                break;
                            case 2:
                                EXPECT_ARMA_EQ( arma::rowvec("0 0 1"), gr );
                                break;
                            case 3:
                                EXPECT_ARMA_EQ( arma::rowvec("-1 -1 -1"), gr );
                                break;
                            }
                        } else {
                            EXPECT_ARMA_EQ( arma::rowvec("0 0 0"), gr );
                        }
                    }
                }
                ++k;
            }
        }

    	/** Declaration of data members **/
//        PatchFETestVector *generic_inst_;                                    ///< pointer to generic object
        FeQArray<Vector> vector_shape_;
        FeQArray<Tensor> grad_vector_shape_;
    };

	PatchFETestVector(std::shared_ptr<DOFHandlerMultiDim> dh, unsigned int quad_order)
    : PatchFETestBase(dh),
      multidim_asm_(this, quad_order)
    {
	    patch_internals_.element_cache_map_.init(patch_internals_.eval_points_);
	    this->initialize(multidim_asm_);
    }

    ~PatchFETestVector() {}

    void test_evaluation() override {
        create_patch(multidim_asm_);
        update_patch(multidim_asm_);

        DHCellAccessor dh_cell = dh_->cell_accessor_from_element( used_element_idx_[0] );
        uint n_dofs = dh_->ds()->fe()[Dim<3>{}]->n_dofs();
        multidim_asm_[3_d]->test_bulk_values(dh_cell, n_dofs);
    }

    MixedPtr<AsmVector, 1> multidim_asm_;                                       ///< Assembly object
};


/**
 * Specialization of tensor FeSystem
 */
class PatchFETestTensor : public PatchFETestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmTensor : public PatchFETestBase::AsmBase<dim> {
    public:
        /// Constructor
    	AsmTensor(PatchFETestBase *generic, uint quad_order)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order),
//          generic_inst_(generic),
          tensor_shape_( this->bulk_integral_->tensor_shape() )
        {}

        /// Destructor
        virtual ~AsmTensor() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int n_dofs) {
            unsigned int k=0;
            for ( auto p : this->bulk_integral_->points( this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx()) ) ) {
                for (unsigned int i=0; i<n_dofs; i++) {
                    // check values
                    arma::mat33 expected_val = this->generic_->expected_tensor_coefs_[k][i%4] * this->generic_->expected_tensor_mats_[i/4];
                    EXPECT_TEST_ARMA_NEAR( expected_val, tensor_shape_.shape(i)(p) );
                }
                ++k;
            }
        }

    	/** Declaration of data members **/
//        PatchFETestVector *generic_inst_;                                    ///< pointer to generic object
        FeQArray<Tensor> tensor_shape_;
    };

	PatchFETestTensor(std::shared_ptr<DOFHandlerMultiDim> dh, unsigned int quad_order)
    : PatchFETestBase(dh),
      multidim_asm_(this, quad_order)
    {
	    patch_internals_.element_cache_map_.init(patch_internals_.eval_points_);
	    this->initialize(multidim_asm_);
    }

    ~PatchFETestTensor() {}

    void test_evaluation() override {
        create_patch(multidim_asm_);
        update_patch(multidim_asm_);

        DHCellAccessor dh_cell = dh_->cell_accessor_from_element( used_element_idx_[0] );
        uint n_dofs = dh_->ds()->fe()[Dim<3>{}]->n_dofs();
        multidim_asm_[3_d]->test_bulk_values(dh_cell, n_dofs);
    }

    MixedPtr<AsmTensor, 1> multidim_asm_;                                       ///< Assembly object
};


/**
 * Specialization of mixed FeSystem
 */
class PatchFETestMixed : public PatchFETestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmMixed : public PatchFETestBase::AsmBase<dim> {
    public:
        /// Constructor
        AsmMixed(PatchFETestBase *generic, uint quad_order, uint quad_order_vec)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order, quad_order_vec),
//          generic_inst_(generic),
          scalar_shape_( this->bulk_integral_->scalar_shape(0) ),
          grad_scalar_shape_( this->bulk_integral_->grad_scalar_shape(0) ),
          vector_shape_( this->bulk_integral_diff_->vector_shape(1) ),
          grad_vector_shape_( this->bulk_integral_diff_->grad_vector_shape(1) ),
          rt_vector_shape_( this->bulk_integral_diff_->vector_shape(2) )
        {}

        /// Destructor
        virtual ~AsmMixed() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int n_dofs_vec, unsigned int n_dofs_rt) {
            for ( auto p : this->bulk_integral_->points( this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx()) ) ) {
                // check value and gradient of P0 function
                EXPECT_TEST_NEAR( 1, scalar_shape_.shape(0)(p) );
                EXPECT_TEST_ARMA_NEAR( arma::vec("0 0 0"), grad_scalar_shape_.shape(0)(p) );
            }

            unsigned int k=0;
            for ( auto p : this->bulk_integral_diff_->points( this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx()) ) ) {
                // check values and gradients of P1^3 function
                for (unsigned int i=0; i<n_dofs_vec; i++) {
                    // check values
                    EXPECT_TEST_ARMA_NEAR( this->generic_->expected_vector_shape_[k][i], vector_shape_.shape(i)(p) );
                    for (unsigned int c=0; c<3; c++)
                    {
                        //check gradients
                        arma::rowvec gr = grad_vector_shape_.shape(i)(p).row(c);
                        if (i / 4 == c) { // gradient of nonzero component
                            switch (i%4) {
                            case 0:
                                EXPECT_ARMA_EQ( arma::rowvec("1 0 0"), gr );
                                break;
                            case 1:
                                EXPECT_ARMA_EQ( arma::rowvec("0 1 0"), gr );
                                break;
                            case 2:
                                EXPECT_ARMA_EQ( arma::rowvec("0 0 1"), gr );
                                break;
                            case 3:
                                EXPECT_ARMA_EQ( arma::rowvec("-1 -1 -1"), gr );
                                break;
                            }
                        } else {
                            EXPECT_ARMA_EQ( arma::rowvec("0 0 0"), gr );
                        }
                    }
                }

                // check RT0 function
                for (unsigned int i=0; i<n_dofs_rt; i++) {
                    EXPECT_TEST_ARMA_NEAR( this->generic_->expected_rt_shape_[k][i], rt_vector_shape_.shape(i)(p) );
                }
                ++k;
            }
        }

    	/** Declaration of data members **/
//        PatchFETestVector *generic_inst_;                                    ///< pointer to generic object
        FeQArray<Scalar> scalar_shape_;
        FeQArray<Vector> grad_scalar_shape_;
        FeQArray<Vector> vector_shape_;
        FeQArray<Tensor> grad_vector_shape_;
        FeQArray<Vector> rt_vector_shape_;
    };

    PatchFETestMixed(std::shared_ptr<DOFHandlerMultiDim> dh, unsigned int quad_order, unsigned int quad_order_vec)
    : PatchFETestBase(dh),
      multidim_asm_(this, quad_order, quad_order_vec)
    {
	    patch_internals_.element_cache_map_.init(patch_internals_.eval_points_);
	    this->initialize(multidim_asm_);
    }

    ~PatchFETestMixed() {}

    void test_evaluation() override {
        //for (auto elm_idx : used_element_idx_)
        {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element( used_element_idx_[0] );
            auto &ppv_bulk = patch_internals_.fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
            this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
            this->add_bulk_integral(dh_cell, this->bulk_integrals_diff_[dh_cell.dim()-1]);
            this->patch_internals_.fe_values_.make_permanent_ppv_data();
        }
        multidim_asm_[1_d]->integrals_.make_permanent();
        multidim_asm_[2_d]->integrals_.make_permanent();
        multidim_asm_[3_d]->integrals_.make_permanent();

        patch_internals_.element_cache_map_.make_paermanent_eval_points();
        patch_internals_.element_cache_map_.create_patch();
        update_patch(multidim_asm_);

        DHCellAccessor dh_cell = dh_->cell_accessor_from_element( used_element_idx_[0] );
        FESystem<3> *fe_sys = dynamic_cast<FESystem<3>*>( dh_->ds()->fe()[Dim<3>{}].get() );
        uint n_dofs_vec_comp = fe_sys->fe()[1]->n_dofs();
        uint n_dofs_rt_comp = fe_sys->fe()[2]->n_dofs();
        multidim_asm_[3_d]->test_bulk_values(dh_cell, n_dofs_vec_comp, n_dofs_rt_comp);
    }

    MixedPtr<AsmMixed, 1> multidim_asm_;                                       ///< Assembly object
};




class FESystemTest : public testing::Test {
public:
    FESystemTest()
    {
        FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
        PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

        std::string input_str = "{ mesh_file=\"fem/one_element_mesh.msh\", optimize_mesh=false }";
        mesh_ = mesh_full_constructor(input_str);
    }

    std::shared_ptr<DOFHandlerMultiDim> create_dh(MixedPtr<FiniteElement> fe) {
        std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe);
        std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh_);
        dh->distribute_dofs(ds);
        return dh;
    }
  
protected:
    Mesh *mesh_;
};






TEST_F(FESystemTest, test_vector) {
    // Test vector-valued FESystem using P1 element on tetrahedron.
    uint quad_order = 1;
    MixedPtr<FE_P> fe_p( quad_order );
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);

    PatchFETestVector patch_fe_system(create_dh(fe), quad_order);
    patch_fe_system.test_evaluation();
}


TEST_F(FESystemTest, test_tensor) {
    // Test vector-valued FESystem using P1 element on tetrahedron.
    uint quad_order = 1;
    MixedPtr<FE_P> fe_p( quad_order );
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FETensor, 9);

    PatchFETestTensor patch_fe_system(create_dh(fe), quad_order);
    patch_fe_system.test_evaluation();
}


TEST_F(FESystemTest, test_mixed_system) {
    // Test mixed-system FE using P0, P1^3 and RT0 elements on tetrahedron.
    // The basis functions are ordered first nodal and then element-supported,
    // hence the scalar constant function from P0 comes after the linear
    // functions from P1^3 and the RT0 functions are at the end.
    uint quad_order = 0, quad_order_vec = 1;
// 	  std::shared_ptr< FiniteElement<0> > fe0_rt = std::make_shared<FE_RT0_disc<0>>();
    std::shared_ptr< FiniteElement<1> > fe1_rt = std::make_shared<FE_RT0_disc<1>>();
    std::shared_ptr< FiniteElement<2> > fe2_rt = std::make_shared<FE_RT0_disc<2>>();
    std::shared_ptr< FiniteElement<3> > fe3_rt = std::make_shared<FE_RT0_disc<3>>();
    std::shared_ptr< FiniteElement<0> > fe0_p = std::make_shared<FE_P<0>>(quad_order);
    std::shared_ptr< FiniteElement<1> > fe1_p = std::make_shared<FE_P<1>>(quad_order);
    std::shared_ptr< FiniteElement<2> > fe2_p = std::make_shared<FE_P<2>>(quad_order);
    std::shared_ptr< FiniteElement<3> > fe3_p = std::make_shared<FE_P<3>>(quad_order);
    std::shared_ptr< FiniteElement<0> > fe0_vec = std::make_shared<FESystem<0>>( std::make_shared<FE_P<0>>(quad_order_vec), FEVector, 3 );
    std::shared_ptr< FiniteElement<1> > fe1_vec = std::make_shared<FESystem<1>>( std::make_shared<FE_P<1>>(quad_order_vec), FEVector, 3 );
    std::shared_ptr< FiniteElement<2> > fe2_vec = std::make_shared<FESystem<2>>( std::make_shared<FE_P<2>>(quad_order_vec), FEVector, 3 );
    std::shared_ptr< FiniteElement<3> > fe3_vec = std::make_shared<FESystem<3>>( std::make_shared<FE_P<3>>(quad_order_vec), FEVector, 3 );
    FESystem<0> fe0_sys( {fe0_p, fe0_vec, fe0_p} );
    FESystem<1> fe1_sys( {fe1_p, fe1_vec, fe1_rt} );
    FESystem<2> fe2_sys( {fe2_p, fe2_vec, fe2_rt} );
    FESystem<3> fe3_sys( {fe3_p, fe3_vec, fe3_rt} );
    MixedPtr<FESystem> fe_sys( std::make_shared<FESystem<0>>(fe0_sys), std::make_shared<FESystem<1>>(fe1_sys),
                               std::make_shared<FESystem<2>>(fe2_sys), std::make_shared<FESystem<3>>(fe3_sys) );

    std::vector<std::vector<unsigned int> > fe_dof_indices = { fe_sys[3_d]->fe_dofs(0), fe_sys[3_d]->fe_dofs(1), fe_sys[3_d]->fe_dofs(2) };
    std::vector<std::vector<unsigned int> > ref_indices = { {0}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, {13, 14, 15, 16} };
    EXPECT_EQ( ref_indices, fe_dof_indices );

    PatchFETestMixed patch_fe_system(create_dh(fe_sys), quad_order, quad_order_vec);
    patch_fe_system.test_evaluation();
}



