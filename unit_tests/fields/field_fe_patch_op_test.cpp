#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "quadrature/quadrature_lib.hh"
#include "fem/integral_acc.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fem/patch_fe_values.hh"
#include "fem/patch_op_impl.hh"
#include "fem/fe_p.hh"
#include "tools/revertable_list.hh"
#include "fields/field_fe.hh"
#include "system/sys_profiler.hh"



// Define EXPECT_<...>_NEAR with fixed abs_error 1e-5
#define EXPECT_TEST_NEAR( A, B )\
  EXPECT_NEAR(A, B, 1e-4)

#define EXPECT_TEST_ARMA_NEAR( A, B ) \
    EXPECT_ARMA_NEAR(A, B, 1e-4);



class FieldFePatchOpTestBase {
public:

    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmBase {
    public:
        /**
         * Constructor
         *
         * @param quad_order      Order of Quadrature (quad_, quad_low_) objects.
         */
        AsmBase(FieldFePatchOpTestBase *generic, uint quad_order)
        : generic_(generic),
          quad_( new QGauss(dim, 2*quad_order) ),
          quad_low_( new QGauss(dim-1, 2*quad_order) ),
          bulk_integral_( create_bulk_integral(quad_) ),
  	      boundary_integral_( create_boundary_integral(quad_low_) )
        {}

    	/// Destructor
        virtual ~AsmBase() {
            delete quad_;
            delete quad_low_;
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

        std::shared_ptr<BoundaryIntegralAcc<dim>> create_boundary_integral(Quadrature *quad) {
            ASSERT_PERMANENT_EQ(quad->dim()+1, dim);
            std::tuple<uint, uint> tpl = IntegralTplHash::integral_tuple(dim, quad->size());
            auto result = integrals_.boundary_.insert({
                    tpl,
                    std::make_shared<BoundaryIntegralAcc<dim>>(generic_->patch_internals_, quad)
                });
            return result.first->second;
        }

        /** Declaration of data members **/
        FieldFePatchOpTestBase *generic_;                                 ///< pointer to generic object
        Quadrature *quad_;                                                ///< Quadrature (of dim)
        Quadrature *quad_low_;                                            ///< Quadrature (of dim-1).
        DimIntegrals<dim> integrals_;                                     ///< Set of used integrals.
        std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;             ///< BulkIntegral
        std::shared_ptr<BoundaryIntegralAcc<dim>> boundary_integral_;     ///< BoundaryIntegral between dim-1 and dim elements

    };


    FieldFePatchOpTestBase(std::shared_ptr<DOFHandlerMultiDim> dh)
    : dh_(dh),
	  patch_internals_(dh_->ds()->fe())
    {
        used_element_idx_ = {0, 1, 2, 3, 8}; // dimension of used elements: 1D, 2D, 2D, 3D, 3D

        ref_scalar_fe_op_ = {
            { 1.31132, 2.6, 2.6, 2.68541, 2.1 },
            { 1.3254, 1.89548, 1.89548, 1.2032, 2.54359 }
        };
        ref_vector_fe_op_ = {
            {1.73397, 1.88868, 2.88868}, {2.6, 3.6, 1.93333}, {3.6, 1.93333, 2.26667}, {1.92918, 2.37639, 2.82361}, {2.65279, 3.1, 2.99443}
        };
        ref_tensor_fe_op_ = {
            { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.232457 },
            { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.683333 },
            { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.183333 },
            { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.428409 },
            { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.290213 }
        };
    }

    ~FieldFePatchOpTestBase() {}

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

    void add_boundary_integral(DHCellAccessor cell, std::shared_ptr<BoundaryIntegral> bdr_integral) {
        if (!cell.is_own()) return; // ghost element

        uint dim = cell.dim();
        auto &ppv_side = patch_internals_.fe_values_.ppv(side_domain, dim);
        auto &ppv_bdr = patch_internals_.fe_values_.ppv(bulk_domain, dim-1);
        for( DHCellSide bdr_side : cell.side_range() ) {
            if ( (bdr_side.side().edge().n_sides() == 1) && (bdr_side.side().is_boundary()) ) { // tests if side is really boundary
                bdr_integral->patch_data().emplace_back(bdr_side);

                unsigned int reg_idx = bdr_side.element().region_idx().idx();
                ++ppv_side.n_mesh_items_;
                ++ppv_bdr.n_mesh_items_;
                for (auto p : bdr_integral->points(bdr_side, &patch_internals_.element_cache_map_) ) {
                    patch_internals_.element_cache_map_.add_eval_point(reg_idx, bdr_side.elem_idx(), p.eval_point_idx(), bdr_side.cell().local_idx());

                    BulkPoint p_bdr = p.point_bdr(bdr_side.cond().element_accessor()); // equivalent point on boundary element
                    unsigned int bdr_reg = bdr_side.cond().element_accessor().region_idx().idx();
                    // invalid local_idx value, DHCellAccessor of boundary element doesn't exist
                    patch_internals_.element_cache_map_.add_eval_point(bdr_reg, bdr_side.cond().bc_ele_idx(), p_bdr.eval_point_idx(), -1);
                }
            }
        }
    }

    virtual void update_patch() =0;

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}

    /// reset patch data
    virtual void reset() =0;

    virtual void set_integrals_arrays() =0;

    void initialize() {
        set_integrals_arrays();
        this->patch_internals_.fe_values_.init_finalize();
    }


    std::shared_ptr<DOFHandlerMultiDim> dh_;

    PatchInternals patch_internals_;                                          ///< Holds common patch objects (EvalPoints, ElementCacheMap ...)
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_;             ///< Bulk integrals of dim 1,2,3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_integrals_;             ///< Edge integrals of dim 1,2,3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_integrals_;     ///< Coupling integrals of dim 1-2,2-3
    std::array<std::shared_ptr<BoundaryIntegral>, 3> boundary_integrals_;     ///< Coupling integrals of dim 1-2,2-3
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_diff_order_;  ///< Bulk integrals of dim 1,2,3 of high order

    std::vector<unsigned int> used_element_idx_;                              ///< List of mesh idx of elements used in tests

    /* Reference values */
    std::vector< std::vector<double> > ref_scalar_fe_op_;
    std::vector<arma::vec3> ref_vector_fe_op_;
    std::vector<arma::mat33> ref_tensor_fe_op_;
};


/**
 * Specialization defining FE scalar operations
 */
class FieldFePatchOpTestScalar : public FieldFePatchOpTestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmScalar : public FieldFePatchOpTestBase::AsmBase<dim> {
    public:
        /// Constructor
        AsmScalar(FieldFePatchOpTestScalar *generic, uint quad_order)
        : FieldFePatchOpTestBase::AsmBase<dim>(generic, quad_order),
          //generic_inst_(generic),
		  scalar_field_fe_op_( this->field_fe_scalar_op() )
        {}

        /// Destructor
        virtual ~AsmScalar() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int i_run, unsigned int i_test_elem) {
            auto p = *( this->bulk_integral_->points(this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx())).begin() );
            double fe_op = this->scalar_field_fe_op_(p);

            EXPECT_TEST_NEAR( fe_op, this->generic_->ref_scalar_fe_op_[i_run][i_test_elem] );
        }

        inline FeQ<Scalar> field_fe_scalar_op()
        {
        	using ShapeSelector = internal::InputOpType<1, 1>;

            VectorMPI data_vec = this->generic_->dh_->create_vector();
            for (uint i=0; i<data_vec.size(); ++i)
                data_vec.set(i, (1.1 + i%3) );
            FieldFeOpData field_fe_op_data(this->generic_->dh_, data_vec);

            std::shared_ptr<FiniteElement<dim>> fe_component = this->generic_->patch_internals_.fe_values_.fe_comp(this->generic_->patch_internals_.fe_[Dim<dim>{}], 0);
            return FeQ<Scalar>(
                this->generic_->patch_internals_.fe_values_.template get<
                    Op::FieldFeOp<dim, Op::BulkDomain, typename ShapeSelector::type<dim, Op::BulkDomain, 3>, 3>,
                    dim
                >(*this->quad_, fe_component, field_fe_op_data)
            );
        }


        /** Declaration of data members **/
        //FieldFePatchOpTestScalar *generic_inst_;     ///< pointer to generic object
        FeQ<Scalar> scalar_field_fe_op_;
    };


    FieldFePatchOpTestScalar(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : FieldFePatchOpTestBase(dh),
      multidim_asm_(this, quad_order)
    {
        patch_internals_.element_cache_map_.init(patch_internals_.eval_points_);
    }

    ~FieldFePatchOpTestScalar() {}

    void set_integrals_arrays() override {
        this->bulk_integrals_[0] = multidim_asm_[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm_[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm_[3_d]->bulk_integral_;
        this->boundary_integrals_[0] = multidim_asm_[1_d]->boundary_integral_;
        this->boundary_integrals_[1] = multidim_asm_[2_d]->boundary_integral_;
        this->boundary_integrals_[2] = multidim_asm_[3_d]->boundary_integral_;
    }

    void update_patch() override {
        patch_internals_.fe_values_.prepare_new_patch(this->patch_internals_.eval_points_);
        patch_internals_.fe_values_.add_patch_points<3>(multidim_asm_[3_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.add_patch_points<2>(multidim_asm_[2_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.add_patch_points<1>(multidim_asm_[1_d]->integrals_, &this->patch_internals_.element_cache_map_);

        START_TIMER("reinit_patch");
        patch_internals_.fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(unsigned int i_run, bool print_tables=false) {
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            auto &ppv_bulk = patch_internals_.fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
            //this->add_boundary_integral(dh_cell, this->boundary_integrals_[dh_cell.dim()-1]);
        	this->patch_internals_.fe_values_.make_permanent_ppv_data();
        }
        multidim_asm_[1_d]->integrals_.make_permanent();
        multidim_asm_[2_d]->integrals_.make_permanent();
        multidim_asm_[3_d]->integrals_.make_permanent();

        patch_internals_.element_cache_map_.make_paermanent_eval_points();
        patch_internals_.element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        update_patch();

        if (print_tables) {
            std::stringstream ss;
            patch_internals_.fe_values_.print_operations(ss);
            WarningOut() << ss.str();
        }

        unsigned int i_test_elem = 0;
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            switch (dh_cell.dim()) {
            case 1:
                multidim_asm_[1_d]->test_bulk_values(dh_cell, i_run, i_test_elem);
                break;
            case 2:
                multidim_asm_[2_d]->test_bulk_values(dh_cell, i_run, i_test_elem);
                break;
            case 3:
                multidim_asm_[3_d]->test_bulk_values(dh_cell, i_run, i_test_elem);
                break;
            }
            ++i_test_elem;
        }
    }

    void reset() override {
        multidim_asm_[1_d]->integrals_.reset();
        multidim_asm_[2_d]->integrals_.reset();
        multidim_asm_[3_d]->integrals_.reset();
        this->patch_internals_.element_cache_map_.clear_element_eval_points_map();
        this->patch_internals_.fe_values_.reset();
    }


    MixedPtr<AsmScalar, 1> multidim_asm_;  ///< Assembly object
};


/**
 * Specialization defining FE vector operations
 */
class FieldFePatchOpTestVector : public FieldFePatchOpTestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmVector : public FieldFePatchOpTestBase::AsmBase<dim> {
    public:
        /// Constructor
    	AsmVector(FieldFePatchOpTestVector *generic, uint quad_order)
        : FieldFePatchOpTestBase::AsmBase<dim>(generic, quad_order),
          //generic_inst_(generic),
		  vector_field_fe_op_( this->field_fe_vector_op() )
        {}

        /// Destructor
        virtual ~AsmVector() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int i_test_elem) {
            auto p = *( this->bulk_integral_->points(this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx())).begin() );
            arma::vec3 fe_op = vector_field_fe_op_(p);

            EXPECT_TEST_ARMA_NEAR( fe_op, this->generic_->ref_vector_fe_op_[i_test_elem] );
        }

        inline FeQ<Vector> field_fe_vector_op()
        {
        	using ShapeSelector = internal::InputOpType<3, 1>;

            VectorMPI data_vec = this->generic_->dh_->create_vector();
            for (uint i=0; i<data_vec.size(); ++i)
                data_vec.set(i, (1.1 + i%4) );
            FieldFeOpData field_fe_op_data(this->generic_->dh_, data_vec);

            std::shared_ptr<FiniteElement<dim>> fe_component = this->generic_->patch_internals_.fe_values_.fe_comp(this->generic_->patch_internals_.fe_[Dim<dim>{}], 0);
            return FeQ<Vector>(
                this->generic_->patch_internals_.fe_values_.template get<
                    Op::FieldFeOp<dim, Op::BulkDomain, typename ShapeSelector::type<dim, Op::BulkDomain, 3>, 3>,
                    dim
                >(*this->quad_, fe_component, field_fe_op_data)
            );
        }


        /** Declaration of data members **/
        //FieldFePatchOpTestVector *generic_inst_;     ///< pointer to generic object
        FeQ<Vector> vector_field_fe_op_;
    };


    FieldFePatchOpTestVector(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : FieldFePatchOpTestBase(dh),
      multidim_asm_(this, quad_order)
    {
        patch_internals_.element_cache_map_.init(patch_internals_.eval_points_);
    }

    ~FieldFePatchOpTestVector() {}

    void set_integrals_arrays() override {
        this->bulk_integrals_[0] = multidim_asm_[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm_[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm_[3_d]->bulk_integral_;
        this->boundary_integrals_[0] = multidim_asm_[1_d]->boundary_integral_;
        this->boundary_integrals_[1] = multidim_asm_[2_d]->boundary_integral_;
        this->boundary_integrals_[2] = multidim_asm_[3_d]->boundary_integral_;
    }

    void update_patch() override {
        patch_internals_.fe_values_.prepare_new_patch(this->patch_internals_.eval_points_);
        patch_internals_.fe_values_.add_patch_points<3>(multidim_asm_[3_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.add_patch_points<2>(multidim_asm_[2_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.add_patch_points<1>(multidim_asm_[1_d]->integrals_, &this->patch_internals_.element_cache_map_);

        START_TIMER("reinit_patch");
        patch_internals_.fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(bool print_tables=false) {
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            auto &ppv_bulk = patch_internals_.fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
        	this->patch_internals_.fe_values_.make_permanent_ppv_data();
        }
        multidim_asm_[1_d]->integrals_.make_permanent();
        multidim_asm_[2_d]->integrals_.make_permanent();
        multidim_asm_[3_d]->integrals_.make_permanent();
        patch_internals_.element_cache_map_.make_paermanent_eval_points();
        patch_internals_.element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        update_patch();

        if (print_tables) {
            std::stringstream ss;
            patch_internals_.fe_values_.print_operations(ss);
            WarningOut() << ss.str();
        }

        unsigned int i_test_elem = 0;
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            switch (dh_cell.dim()) {
            case 1:
                multidim_asm_[1_d]->test_bulk_values(dh_cell, i_test_elem);
                break;
            case 2:
                multidim_asm_[2_d]->test_bulk_values(dh_cell, i_test_elem);
                break;
            case 3:
                multidim_asm_[3_d]->test_bulk_values(dh_cell, i_test_elem);
                break;
            }
            ++i_test_elem;
        }
    }

    void reset() override {
        multidim_asm_[1_d]->integrals_.reset();
        multidim_asm_[2_d]->integrals_.reset();
        multidim_asm_[3_d]->integrals_.reset();
        this->patch_internals_.element_cache_map_.clear_element_eval_points_map();
        this->patch_internals_.fe_values_.reset();
    }


    MixedPtr<AsmVector, 1> multidim_asm_;  ///< Assembly object
};


/**
 * Specialization defining FE tensor operations
 */
class FieldFePatchOpTestTensor : public FieldFePatchOpTestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmTensor : public FieldFePatchOpTestBase::AsmBase<dim> {
    public:
        /// Constructor
    	AsmTensor(FieldFePatchOpTestTensor *generic, uint quad_order)
        : FieldFePatchOpTestBase::AsmBase<dim>(generic, quad_order),
          //generic_inst_(generic),
		  tensor_field_fe_op_( this->field_fe_tensor_op() )
        {}

        /// Destructor
        virtual ~AsmTensor() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int i_test_elem) {
            auto p = *( this->bulk_integral_->points(this->generic_->patch_internals_.element_cache_map_.position_in_cache(dh_cell.elm_idx())).begin() );
            arma::mat33 fe_op = tensor_field_fe_op_(p);

            EXPECT_TEST_ARMA_NEAR( fe_op, this->generic_->ref_tensor_fe_op_[i_test_elem] );
        }

        inline FeQ<Tensor> field_fe_tensor_op()
        {
        	using ShapeSelector = internal::InputOpType<3, 3>;

            VectorMPI data_vec = this->generic_->dh_->create_vector();
            for (uint i=0; i<data_vec.size(); ++i)
                data_vec.set(i, (1.1 + i%5) );
            FieldFeOpData field_fe_op_data(this->generic_->dh_, data_vec);

            std::shared_ptr<FiniteElement<dim>> fe_component = this->generic_->patch_internals_.fe_values_.fe_comp(this->generic_->patch_internals_.fe_[Dim<dim>{}], 0);
            return FeQ<Tensor>(
                this->generic_->patch_internals_.fe_values_.template get<
                    Op::FieldFeOp<dim, Op::BulkDomain, typename ShapeSelector::type<dim, Op::BulkDomain, 3>, 3>,
                    dim
                >(*this->quad_, fe_component, field_fe_op_data)
            );
        }


        /** Declaration of data members **/
        //FieldFePatchOpTestVector *generic_inst_;     ///< pointer to generic object
        FeQ<Tensor> tensor_field_fe_op_;
    };


    FieldFePatchOpTestTensor(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : FieldFePatchOpTestBase(dh),
      multidim_asm_(this, quad_order)
    {
        patch_internals_.element_cache_map_.init(patch_internals_.eval_points_);
    }

    ~FieldFePatchOpTestTensor() {}

    void set_integrals_arrays() override {
        this->bulk_integrals_[0] = multidim_asm_[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm_[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm_[3_d]->bulk_integral_;
        this->boundary_integrals_[0] = multidim_asm_[1_d]->boundary_integral_;
        this->boundary_integrals_[1] = multidim_asm_[2_d]->boundary_integral_;
        this->boundary_integrals_[2] = multidim_asm_[3_d]->boundary_integral_;
    }

    void update_patch() override {
        patch_internals_.fe_values_.prepare_new_patch(this->patch_internals_.eval_points_);
        patch_internals_.fe_values_.add_patch_points<3>(multidim_asm_[3_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.add_patch_points<2>(multidim_asm_[2_d]->integrals_, &this->patch_internals_.element_cache_map_);
        patch_internals_.fe_values_.add_patch_points<1>(multidim_asm_[1_d]->integrals_, &this->patch_internals_.element_cache_map_);

        START_TIMER("reinit_patch");
        patch_internals_.fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(bool print_tables=false) {
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            auto &ppv_bulk = patch_internals_.fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
        	this->patch_internals_.fe_values_.make_permanent_ppv_data();
        }
        multidim_asm_[1_d]->integrals_.make_permanent();
        multidim_asm_[2_d]->integrals_.make_permanent();
        multidim_asm_[3_d]->integrals_.make_permanent();
        patch_internals_.element_cache_map_.make_paermanent_eval_points();
        patch_internals_.element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        update_patch();

        if (print_tables) {
            std::stringstream ss;
            patch_internals_.fe_values_.print_operations(ss);
            WarningOut() << ss.str();
        }

        unsigned int i_test_elem = 0;
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            switch (dh_cell.dim()) {
            case 1:
                multidim_asm_[1_d]->test_bulk_values(dh_cell, i_test_elem);
                break;
            case 2:
                multidim_asm_[2_d]->test_bulk_values(dh_cell, i_test_elem);
                break;
            case 3:
                multidim_asm_[3_d]->test_bulk_values(dh_cell, i_test_elem);
                break;
            }
            ++i_test_elem;
        }
    }

    void reset() override {
        multidim_asm_[1_d]->integrals_.reset();
        multidim_asm_[2_d]->integrals_.reset();
        multidim_asm_[3_d]->integrals_.reset();
        this->patch_internals_.element_cache_map_.clear_element_eval_points_map();
        this->patch_internals_.fe_values_.reset();
    }


    MixedPtr<AsmTensor, 1> multidim_asm_;  ///< Assembly object
};


/// Complete test with FE scalar operations
void compare_evaluation_func_scalar(Mesh* mesh, unsigned int i_run, bool print_fa_data = false) {
    std::vector<uint> quad_orders = {1, 2};
    MixedPtr<FE_P_disc> fe(quad_orders[i_run]);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    FieldFePatchOpTestScalar patch_fe(quad_orders[i_run], dh);
    patch_fe.initialize();
    patch_fe.test_evaluation(i_run, print_fa_data);
    patch_fe.reset();
    patch_fe.test_evaluation(i_run);
}


/// Complete test with FE vector operations
void compare_evaluation_func_vector(Mesh* mesh, unsigned int quad_order, bool print_fa_data = false) {
    MixedPtr<FE_P> fe_p( quad_order );
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    FieldFePatchOpTestVector patch_fe(quad_order, dh);
    patch_fe.initialize();
    patch_fe.test_evaluation(print_fa_data);
    patch_fe.reset();
    patch_fe.test_evaluation();
}


/// Complete test with FE vector operations
void compare_evaluation_func_tensor(Mesh* mesh, unsigned int quad_order, bool print_fa_data = false) {
    MixedPtr<FE_P> fe_p( quad_order );
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FETensor, 9);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    FieldFePatchOpTestTensor patch_fe(quad_order, dh);
    patch_fe.initialize();
    patch_fe.test_evaluation(print_fa_data);
    patch_fe.reset();
    patch_fe.test_evaluation();
}


TEST(FieldFePatchOpTest, complete_evaluation) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    std::string input_str = "{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }";
    Mesh* mesh = mesh_full_constructor(input_str);

    // two tests with different quad_order and Scalar / Vector FE operations
    compare_evaluation_func_scalar(mesh, 0, true);
    compare_evaluation_func_scalar(mesh, 1);
    compare_evaluation_func_vector(mesh, 1, true);
    compare_evaluation_func_tensor(mesh, 1);
}

