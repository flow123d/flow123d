#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "quadrature/quadrature_lib.hh"
#include "fem/integral_acc.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fem/fe_values.hh"
#include "fem/patch_fe_values.hh"
#include "fem/patch_op_impl.hh"
#include "fem/fe_p.hh"
#include "tools/revertable_list.hh"
#include "system/sys_profiler.hh"



class PatchFETestBase {
public:

	/**
	 * Helper structzre holds data of cell (bulk) integral
	 *
	 * Data is specified by cell and subset index in EvalPoint object
	 */
    struct BulkIntegralData {
    	/// Default constructor
        BulkIntegralData() {}

        /// Constructor with data mebers initialization
        BulkIntegralData(DHCellAccessor dhcell, unsigned int subset_idx)
        : cell(dhcell), subset_index(subset_idx) {}

        /// Copy constructor
        BulkIntegralData(const BulkIntegralData &other)
        : cell(other.cell), subset_index(other.subset_index) {}

        DHCellAccessor cell;          ///< Specified cell (element)
        unsigned int subset_index;    ///< Index (order) of subset in EvalPoints object
    };

	/**
	 * Helper structzre holds data of edge integral
	 *
	 * Data is specified by side and subset index in EvalPoint object
	 */
    struct EdgeIntegralData {
    	/// Default constructor
    	EdgeIntegralData()
    	: edge_side_range(make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() ), make_iter<DHEdgeSide, DHCellSide>( DHEdgeSide() )) {}

        /// Copy constructor
    	EdgeIntegralData(const EdgeIntegralData &other)
        : edge_side_range(other.edge_side_range), subset_index(other.subset_index) {}

        /// Constructor with data mebers initialization
    	EdgeIntegralData(RangeConvert<DHEdgeSide, DHCellSide> range, unsigned int subset_idx)
        : edge_side_range(range), subset_index(subset_idx) {}

    	RangeConvert<DHEdgeSide, DHCellSide> edge_side_range;   ///< Specified cell side (element)
        unsigned int subset_index;                              ///< Index (order) of subset in EvalPoints object
	};

	/**
	 * Helper structzre holds data of neighbour (coupling) integral
	 *
	 * Data is specified by cell, side and their subset indices in EvalPoint object
	 */
    struct CouplingIntegralData {
    	/// Default constructor
       	CouplingIntegralData() {}

        /// Constructor with data mebers initialization
       	CouplingIntegralData(DHCellAccessor dhcell, unsigned int bulk_idx, DHCellSide dhside, unsigned int side_idx)
        : cell(dhcell), bulk_subset_index(bulk_idx), side(dhside), side_subset_index(side_idx) {}

        /// Copy constructor
       	CouplingIntegralData(const CouplingIntegralData &other)
        : cell(other.cell), bulk_subset_index(other.bulk_subset_index), side(other.side), side_subset_index(other.side_subset_index) {}

        DHCellAccessor cell;
	    unsigned int bulk_subset_index;    ///< Index (order) of lower dim subset in EvalPoints object
        DHCellSide side;                   ///< Specified cell side (higher dim element)
	    unsigned int side_subset_index;    ///< Index (order) of higher dim subset in EvalPoints object
    };

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
        AsmBase(PatchFETestBase *generic, uint quad_order, uint quad_diff_order)
        : generic_(generic),
          quad_( new QGauss(dim, 2*quad_order) ),
          quad_low_( new QGauss(dim-1, 2*quad_order) ),
          quad_diff_order_( new QGauss(dim, 2*quad_diff_order) ),
          quad_low_diff_order_( new QGauss(dim-1, 2*quad_diff_order) ),
          bulk_integral_( std::make_shared<BulkIntegralAcc<dim>>(generic_->eval_points_, quad_, &generic_->patch_fe_values_, &generic_->element_cache_map_) ),
          bulk_integral_diff_order_( std::make_shared<BulkIntegralAcc<dim>>(generic_->eval_points_, quad_diff_order_, &generic_->patch_fe_values_, &generic_->element_cache_map_) ),
  	      edge_integral_( std::make_shared<EdgeIntegralAcc<dim>>(generic_->eval_points_, quad_low_, &generic_->patch_fe_values_, &generic_->element_cache_map_) ),
  	      coupling_integral_( std::make_shared<CouplingIntegralAcc<dim>>(generic_->eval_points_, quad_, &generic_->patch_fe_values_, &generic_->element_cache_map_) ),
	      det_( bulk_integral_->determinant() ),
	      jxw_( bulk_integral_->JxW() ),
	      jxw_side_( edge_integral_->JxW() ),
	      normal_vec_( edge_integral_->normal_vector() )
        {}

    	/// Destructor
        virtual ~AsmBase() {
            delete quad_;
            delete quad_low_;
            delete quad_diff_order_;
            delete quad_low_diff_order_;
        }


    	/** Declaration of data members **/
        PatchFETestBase *generic_;                                        ///< pointer to generic object
        Quadrature *quad_;                                                ///< Quadrature (of dim)
        Quadrature *quad_low_;                                            ///< Quadrature (of dim-1).
        Quadrature *quad_diff_order_;                                     ///< Quadrature of different size than previous (of dim)
        Quadrature *quad_low_diff_order_;                                 ///< Quadrature of different size than previous (of dim-1).
        std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_;             ///< BulkIntegral
        std::shared_ptr<BulkIntegralAcc<dim>> bulk_integral_diff_order_;  ///< BulkIntegralof high order
        std::shared_ptr<EdgeIntegralAcc<dim>> edge_integral_;             ///< EdgeIntegral
        std::shared_ptr<CouplingIntegralAcc<dim>> coupling_integral_;     ///< CouplingIntegral between dim and dim+1 elements
        ElQ<Scalar> det_;
        FeQ<Scalar> jxw_;
        FeQ<Scalar> jxw_side_;
        ElQ<Vector> normal_vec_;

    };


    PatchFETestBase(std::shared_ptr<DOFHandlerMultiDim> dh)
    : dh_(dh), patch_fe_values_(dh_->ds()->fe()),
      fe_(dh_->ds()->fe()), fe_values_(3), fe_values_side_(3),
	  eval_points_( std::make_shared<EvalPoints>() ),
      bulk_integral_data_(20, 10),
      edge_integral_data_(12, 6),
      coupling_integral_data_(12, 6)
    {}

    ~PatchFETestBase() {}

    /// Return BulkPoint range of appropriate dimension
    inline Range< BulkPoint > bulk_points(unsigned int dim, unsigned int element_patch_idx) const {
        return bulk_integrals_[dim-1]->points(element_patch_idx, &element_cache_map_);
    }

    /// Return EdgePoint range of appropriate dimension
    inline Range< EdgePoint > edge_points(unsigned int dim, const DHCellSide &cell_side) const {
	    return edge_integrals_[dim-1]->points(cell_side, &element_cache_map_);
    }

    /// Return CouplingPoint range of appropriate dimension
    inline Range< CouplingPoint > coupling_points(unsigned int dim, const DHCellSide &cell_side) const {
        ASSERT( cell_side.dim() > 1 ).error("Invalid cell dimension, must be 2 or 3!\n");
	    return coupling_integrals_[dim-2]->points(cell_side, &element_cache_map_);
    }

    void add_bulk_integral(DHCellAccessor cell, std::shared_ptr<BulkIntegral> bulk_integral) {
        uint subset_idx = bulk_integral->get_subset_idx();
        bulk_integral_data_.emplace_back(cell, subset_idx);
        uint dim = cell.dim();
        auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, dim);

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( eval_points_->subset_begin(dim, subset_idx) );
                  i<uint( eval_points_->subset_end(dim, subset_idx) ); ++i) {
            element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }
    }

    void add_edge_integral(DHCellAccessor cell, std::shared_ptr<EdgeIntegral> edge_integral) {
        uint dim = cell.dim();
        auto &ppv_edge = patch_fe_values_.ppv(side_domain, dim);
        for( DHCellSide cell_side : cell.side_range() ) {
            if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                auto range = cell_side.edge_sides();
                uint subset_idx = edge_integral->get_subset_idx();
                edge_integral_data_.emplace_back(range, subset_idx);

                for( DHCellSide edge_side : range ) {
                    uint dim = edge_side.dim();
                    ++ppv_edge.n_mesh_items_;
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : edge_integral->points(edge_side, &element_cache_map_) ) {
                        element_cache_map_.add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    }
                }
            }
        }
    }

//    void add_coupling_integral(DHCellAccessor cell, std::shared_ptr<CouplingIntegral> coupling_integral) {
//        bool add_low = true;
//        uint dim = cell.dim();
//        for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
//            if (cell.dim() != neighb_side.dim()-1) continue;
//            coupling_integral_data_.emplace_back(cell, coupling_integral->get_subset_low_idx(), neighb_side,
//                    coupling_integral>get_subset_high_idx());
//            table_sizes_.elem_sizes_[1][cell.dim()]++;
//            if (add_low) table_sizes_.elem_sizes_[0][cell.dim()-1]++;
//
//            unsigned int reg_idx_low = cell.elm().region_idx().idx();
//            unsigned int reg_idx_high = neighb_side.element().region_idx().idx();
//            for (auto p : coupling_integral->points(neighb_side, &element_cache_map_) ) {
//                element_cache_map_.add_eval_point(reg_idx_high, neighb_side.elem_idx(), p.eval_point_idx(), neighb_side.cell().local_idx());
//                table_sizes_.point_sizes_[1][cell.dim()]++;
//
//                if (add_low) {
//                    auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
//                    element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
//                    table_sizes_.point_sizes_[0][cell.dim()-1]++;
//                }
//            }
//            add_low = false;
//        }
//    }

    void update_patch() {
        patch_fe_values_.resize_tables(eval_points_);
        for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            uint dim = bulk_integral_data_[i].cell.dim();
            uint element_patch_idx = element_cache_map_.position_in_cache(bulk_integral_data_[i].cell.elm_idx());
            uint elm_pos = patch_fe_values_.register_element(bulk_integral_data_[i].cell, element_patch_idx);
            uint i_point = 0;
            for (auto p : this->bulk_points(dim, element_patch_idx) ) {
                patch_fe_values_.register_bulk_point(bulk_integral_data_[i].cell, elm_pos, p.value_cache_idx(), i_point++);
            }
        }
        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            auto range = edge_integral_data_[i].edge_side_range;
            uint dim = range.begin()->dim();
            for( DHCellSide edge_side : range )
            {
                uint side_pos = patch_fe_values_.register_side(edge_side, &element_cache_map_);
                uint i_point = 0;
                for (auto p : this->edge_points(dim, edge_side) ) {
                    patch_fe_values_.register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
        uint element_patch_idx, elm_pos=0;
        uint last_element_idx = -1;
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            uint dim = coupling_integral_data_[i].side.dim();
            uint side_pos = patch_fe_values_.register_side(coupling_integral_data_[i].side, &element_cache_map_);
            if (coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                element_patch_idx = this->element_cache_map_.position_in_cache(coupling_integral_data_[i].cell.elm_idx());
                elm_pos = patch_fe_values_.register_element(coupling_integral_data_[i].cell, element_patch_idx);
            }

            uint i_bulk_point = 0, i_side_point = 0;
            for (auto p_high : this->coupling_points(dim, coupling_integral_data_[i].side) )
            {
                patch_fe_values_.register_side_point(coupling_integral_data_[i].side, side_pos, p_high.value_cache_idx(), i_side_point++);
                if (coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                    auto p_low = p_high.lower_dim(coupling_integral_data_[i].cell);
                    patch_fe_values_.register_bulk_point(coupling_integral_data_[i].cell, elm_pos, p_low.value_cache_idx(), i_bulk_point++);
                }
            }
            last_element_idx = coupling_integral_data_[i].cell.elm_idx();
        }

        this->reinit_patch_fe();
    }

    virtual void reinit_patch_fe() =0;

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}

    /// reset patch data
    void reset() {
        bulk_integral_data_.reset();
        edge_integral_data_.reset();
        coupling_integral_data_.reset();
        element_cache_map_.clear_element_eval_points_map();
        patch_fe_values_.reset();
    }


    std::shared_ptr<DOFHandlerMultiDim> dh_;
    PatchFEValues<3> patch_fe_values_;                                        ///< Common FEValues object over all dimensions

    MixedPtr<FiniteElement> fe_;
    std::vector<FEValues<3>> fe_values_;                                      ///< FeValues object of elements of dim 1,2,3
    std::vector<FEValues<3>> fe_values_side_;                                 ///< FeValues object of sides of dim 0,1,2

    std::shared_ptr<EvalPoints> eval_points_;                                 ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                                       ///< ElementCacheMap according to EvalPoints
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_;             ///< Bulk integrals of dim 1,2,3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_integrals_;             ///< Edge integrals of dim 1,2,3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_integrals_;     ///< Coupling integrals of dim 1-2,2-3
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_diff_order_;  ///< Bulk integrals of dim 1,2,3 of high order
    RevertableList<BulkIntegralData> bulk_integral_data_;                     ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData> edge_integral_data_;                     ///< Holds data for computing edge integrals.
    RevertableList<CouplingIntegralData> coupling_integral_data_;             ///< Holds data for computing couplings integrals.
};


/**
 * Specialization defining FE scalar operations
 */
class PatchFETestScalar : public PatchFETestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmScalar : public PatchFETestBase::AsmBase<dim> {
    public:
        /// Constructor
        AsmScalar(PatchFETestBase *generic, uint quad_order, uint quad_diff_order)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order, quad_diff_order),
          //generic_inst_(generic),
	      scalar_shape_( this->bulk_integral_->scalar_shape() ),
	      scalar_shape_side_( this->edge_integral_->scalar_shape() ),
	      grad_scalar_shape_( this->bulk_integral_->grad_scalar_shape() ),
	      grad_scalar_shape_side_( this->edge_integral_->grad_scalar_shape() ),
	      conc_join_shape_( FeQJoin<Scalar>( this->coupling_integral_->scalar_join_shape() ) )
        {}

        /// Destructor
        virtual ~AsmScalar() {}


    	/** Declaration of data members **/
//        PatchFETestScalar *generic_inst_;                                    ///< pointer to generic object
        FeQArray<Scalar> scalar_shape_;
        FeQArray<Scalar> scalar_shape_side_;
        FeQArray<Vector> grad_scalar_shape_;
        FeQArray<Vector> grad_scalar_shape_side_;
        FeQJoin<Scalar> conc_join_shape_;
    };


    PatchFETestScalar(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : PatchFETestBase(dh),
      multidim_asm_(this, quad_order, 0)
    {
        element_cache_map_.init(eval_points_);
    }

    ~PatchFETestScalar() {}

    void set_integrals_arrays() {
        this->bulk_integrals_[0] = multidim_asm_[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm_[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm_[3_d]->bulk_integral_;
        this->edge_integrals_[0] = multidim_asm_[1_d]->edge_integral_;
        this->edge_integrals_[1] = multidim_asm_[2_d]->edge_integral_;
        this->edge_integrals_[2] = multidim_asm_[3_d]->edge_integral_;
        this->coupling_integrals_[0] = multidim_asm_[1_d]->coupling_integral_;
        this->coupling_integrals_[1] = multidim_asm_[2_d]->coupling_integral_;
        this->bulk_integrals_diff_order_[0] = multidim_asm_[1_d]->bulk_integral_diff_order_;
        this->bulk_integrals_diff_order_[1] = multidim_asm_[2_d]->bulk_integral_diff_order_;
        this->bulk_integrals_diff_order_[2] = multidim_asm_[3_d]->bulk_integral_diff_order_;
    }

    void initialize() {
        set_integrals_arrays();

        UpdateFlags u = update_values | update_inverse_jacobians | update_JxW_values | update_quadrature_points | update_volume_elements | update_gradients;
        UpdateFlags u_side = update_values | update_inverse_jacobians | update_side_JxW_values | update_normal_vectors | update_quadrature_points | update_gradients;
        fe_values_[0].initialize(*multidim_asm_[1_d]->quad_, *fe_[Dim<1>{}], u);
        fe_values_[1].initialize(*multidim_asm_[2_d]->quad_, *fe_[Dim<2>{}], u);
        fe_values_[2].initialize(*multidim_asm_[3_d]->quad_, *fe_[Dim<3>{}], u);
        fe_values_side_[0].initialize(*multidim_asm_[1_d]->quad_low_, *fe_[Dim<1>{}], u_side);
        fe_values_side_[1].initialize(*multidim_asm_[2_d]->quad_low_, *fe_[Dim<2>{}], u_side);
        fe_values_side_[2].initialize(*multidim_asm_[3_d]->quad_low_, *fe_[Dim<3>{}], u_side);

        this->patch_fe_values_.init_finalize();
    }

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(unsigned int i_run, bool print_tables=false) {
        std::vector< std::vector<double> > ref_jxw = {
            { 1.73205080756887720, 0.94280904158206336, 0.94280904158206336, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331 },
	        { 0.96225044864937626, 0.63181854741421128, 0.63181854741421128, 0.09799072415514927, 0.09799072415514927, 0.09799072415514927, 0.09799072415514927, 0.09799072415514927, 0.09799072415514927 }
        };

        for(auto cell_it = dh_->local_range().begin(); cell_it != dh_->local_range().end(); ++cell_it) {
            auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, cell_it->dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(*cell_it, this->bulk_integrals_[cell_it->dim()-1]);
        	this->add_edge_integral(*cell_it, this->edge_integrals_[cell_it->dim()-1]);
//            this->add_coupling_integral(*cell_it, this->coupling_integrals_[cell_it->dim()-1]);
        	this->patch_fe_values_.make_permanent_ppv_data();
        }
        bulk_integral_data_.make_permanent();
        edge_integral_data_.make_permanent();
        coupling_integral_data_.make_permanent();
        element_cache_map_.make_paermanent_eval_points();
        element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        update_patch();

        if (print_tables) {
            std::stringstream ss;
            patch_fe_values_.print_operations(ss);
            WarningOut() << ss.str();
        }

        for(auto dh_cell : dh_->local_range() ) {
            ElementAccessor<3> elm = dh_cell.elm();
            auto p = *( bulk_integrals_[dh_cell.dim()-1]->points(element_cache_map_.position_in_cache(dh_cell.elm_idx()), &element_cache_map_).begin() );
            double jxw = 0.0, det = 0.0, det_ref = 0.0;
            double scalar_shape_dof0 = 0.0, scalar_shape_dof0_ref = 0.0;
            double scalar_shape_dof1 = 0.0, scalar_shape_dof1_ref = 0.0;
            arma::vec3 grad_scalar_dof0("0 0 0");
            arma::vec3 grad_scalar_dof0_ref("0 0 0");
            arma::vec3 grad_scalar_dof1("0 0 0");
            arma::vec3 grad_scalar_dof1_ref("0 0 0");
            switch (dh_cell.dim()) {
            case 1:
                fe_values_[0].reinit(elm);
                jxw = multidim_asm_[1_d]->jxw_(p);
                det = multidim_asm_[1_d]->det_(p);
                scalar_shape_dof0 = multidim_asm_[1_d]->scalar_shape_.shape(0)(p);
                grad_scalar_dof0 = multidim_asm_[1_d]->grad_scalar_shape_.shape(0)(p);
                det_ref = fe_values_[0].determinant(0);
                scalar_shape_dof0_ref = fe_values_[0].shape_value(0, 0);
                grad_scalar_dof0_ref = fe_values_[0].shape_grad(0, 0);
                break;
            case 2:
                fe_values_[1].reinit(elm);
                jxw = multidim_asm_[2_d]->jxw_(p);
                det = multidim_asm_[2_d]->det_(p);
                scalar_shape_dof0 = multidim_asm_[2_d]->scalar_shape_.shape(0)(p);
                scalar_shape_dof1 = multidim_asm_[2_d]->scalar_shape_.shape(1)(p);
                grad_scalar_dof0 = multidim_asm_[2_d]->grad_scalar_shape_.shape(0)(p);
                grad_scalar_dof1 = multidim_asm_[2_d]->grad_scalar_shape_.shape(1)(p);
                det_ref = fe_values_[1].determinant(0);
                scalar_shape_dof0_ref = fe_values_[1].shape_value(0, 0);
                scalar_shape_dof1_ref = fe_values_[1].shape_value(1, 0);
                grad_scalar_dof0_ref = fe_values_[1].shape_grad(0, 0);
                grad_scalar_dof1_ref = fe_values_[1].shape_grad(1, 0);
                break;
            case 3:
                fe_values_[2].reinit(elm);
                jxw = multidim_asm_[3_d]->jxw_(p);
                det = multidim_asm_[3_d]->det_(p);
                scalar_shape_dof0 = multidim_asm_[3_d]->scalar_shape_.shape(0)(p);
                scalar_shape_dof1 = multidim_asm_[3_d]->scalar_shape_.shape(1)(p);
                grad_scalar_dof0 = multidim_asm_[3_d]->grad_scalar_shape_.shape(0)(p);
                grad_scalar_dof1 = multidim_asm_[3_d]->grad_scalar_shape_.shape(1)(p);
                det_ref = fe_values_[2].determinant(0);
                scalar_shape_dof0_ref = fe_values_[2].shape_value(0, 0);
                scalar_shape_dof1_ref = fe_values_[2].shape_value(1, 0);
                grad_scalar_dof0_ref = fe_values_[2].shape_grad(0, 0);
                grad_scalar_dof1_ref = fe_values_[2].shape_grad(1, 0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, ref_jxw[i_run][dh_cell.elm_idx()] );
//            std::cout << "Elem: " << dh_cell.elm_idx() << ", dim " << dh_cell.dim() << ", det " << det << ", ref " << det_ref << std::endl;
            EXPECT_DOUBLE_EQ( det, det_ref );
            EXPECT_DOUBLE_EQ( scalar_shape_dof0, scalar_shape_dof0_ref );
            EXPECT_DOUBLE_EQ( scalar_shape_dof1, scalar_shape_dof1_ref );
            EXPECT_ARMA_EQ( grad_scalar_dof0, grad_scalar_dof0_ref );
            EXPECT_ARMA_EQ( grad_scalar_dof1, grad_scalar_dof1_ref );
        }

        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            auto range = edge_integral_data_[i].edge_side_range;

            auto zero_edge_side = *range.begin();
            auto p = *( edge_integrals_[zero_edge_side.dim()-1]->points(zero_edge_side, &element_cache_map_).begin() );

            double jxw = 0.0, jxw_ref = 0.0;
            double scalar_shape = 0.0, scalar_shape_ref = 0.0;
            arma::vec3 normal_vec = {0.0, 0.0, 0.0};
            arma::vec3 normal_vec_ref = {0.0, 0.0, 0.0};
            arma::vec3 grad_scalar("0 0 0");
            arma::vec3 grad_scalar_ref("0 0 0");
            switch (zero_edge_side.dim()) {
            case 1:
                jxw = multidim_asm_[1_d]->jxw_side_(p);
                normal_vec = multidim_asm_[1_d]->normal_vec_(p);
                scalar_shape = multidim_asm_[1_d]->scalar_shape_side_.shape(0)(p);
                grad_scalar = multidim_asm_[1_d]->grad_scalar_shape_side_.shape(0)(p);
                fe_values_side_[0].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[0].JxW(0);
                normal_vec_ref = fe_values_side_[0].normal_vector(0);
                scalar_shape_ref = fe_values_side_[0].shape_value(0, 0);
                grad_scalar_ref = fe_values_side_[0].shape_grad(0, 0);
                break;
            case 2:
                jxw = multidim_asm_[2_d]->jxw_side_(p);
                normal_vec = multidim_asm_[2_d]->normal_vec_(p);
                scalar_shape = multidim_asm_[2_d]->scalar_shape_side_.shape(0)(p);
                grad_scalar = multidim_asm_[2_d]->grad_scalar_shape_side_.shape(0)(p);
                fe_values_side_[1].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[1].JxW(0);
                normal_vec_ref = fe_values_side_[1].normal_vector(0);
                scalar_shape_ref = fe_values_side_[1].shape_value(0, 0);
                grad_scalar_ref = fe_values_side_[1].shape_grad(0, 0);
                break;
            case 3:
                jxw = multidim_asm_[3_d]->jxw_side_(p);
                normal_vec = multidim_asm_[3_d]->normal_vec_(p);
                scalar_shape = multidim_asm_[3_d]->scalar_shape_side_.shape(0)(p);
                grad_scalar = multidim_asm_[3_d]->grad_scalar_shape_side_.shape(0)(p);
                fe_values_side_[2].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[2].JxW(0);
                normal_vec_ref = fe_values_side_[2].normal_vector(0);
                scalar_shape_ref = fe_values_side_[2].shape_value(0, 0);
                grad_scalar_ref = fe_values_side_[2].shape_grad(0, 0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, jxw_ref );
            EXPECT_ARMA_EQ( normal_vec, normal_vec_ref );
            EXPECT_DOUBLE_EQ( scalar_shape, scalar_shape_ref );
            EXPECT_ARMA_EQ( grad_scalar, grad_scalar_ref );
        }

        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            DHCellAccessor cell_lower_dim = coupling_integral_data_[i].cell;
            DHCellSide neighb_side = coupling_integral_data_[i].side;;
            //std::cout << " el high " << neighb_side.elem_idx() << ", el low: " << cell_lower_dim.elm_idx() << std::endl;

            auto p_high = *( coupling_points(neighb_side.dim(), neighb_side).begin() );
            auto p_low = p_high.lower_dim(cell_lower_dim);

            uint i_dof_high=0, i_dof_low=0;
            switch (neighb_side.dim()) {
            case 2:
                fe_values_[0].reinit(cell_lower_dim.elm());
                fe_values_side_[1].reinit(neighb_side.side());
                for (uint i_dof=0; i_dof<multidim_asm_[1_d]->conc_join_shape_.n_dofs_both(); ++i_dof) {
                    if (multidim_asm_[1_d]->conc_join_shape_.is_high_dim(i_dof)) {
                        auto result = multidim_asm_[1_d]->conc_join_shape_.shape(i_dof)(p_high);
                        auto ref = fe_values_side_[1].shape_value(i_dof_high, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_low = multidim_asm_[1_d]->conc_join_shape_.shape(i_dof)(p_low);
                        EXPECT_DOUBLE_EQ( result_low, 0.0 );
                    	i_dof_high++;
                    }
                    else {
                        auto result = multidim_asm_[1_d]->conc_join_shape_.shape(i_dof)(p_low);
                        auto ref = fe_values_[0].shape_value(i_dof_low, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_high = multidim_asm_[1_d]->conc_join_shape_.shape(i_dof)(p_high);
                        EXPECT_DOUBLE_EQ( result_high, 0.0 );
                        i_dof_low++;
                    }
                }
                break;
            case 3:
                fe_values_[1].reinit(cell_lower_dim.elm());
                fe_values_side_[2].reinit(neighb_side.side());
                for (uint i_dof=0; i_dof<multidim_asm_[2_d]->conc_join_shape_.n_dofs_both(); ++i_dof) {
                    if (multidim_asm_[2_d]->conc_join_shape_.is_high_dim(i_dof)) {
                        auto result = multidim_asm_[2_d]->conc_join_shape_.shape(i_dof)(p_high);
                        auto ref = fe_values_side_[2].shape_value(i_dof_high, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_low = multidim_asm_[2_d]->conc_join_shape_.shape(i_dof)(p_low);
                        EXPECT_DOUBLE_EQ( result_low, 0.0 );
                        i_dof_high++;
                    }
                    else {
                   	    auto result = multidim_asm_[2_d]->conc_join_shape_.shape(i_dof)(p_low);
                        auto ref = fe_values_[1].shape_value(i_dof_low, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_high = multidim_asm_[2_d]->conc_join_shape_.shape(i_dof)(p_high);
                        EXPECT_DOUBLE_EQ( result_high, 0.0 );
                        i_dof_low++;
                    }
                }
                break;
            }
        }

    }

    MixedPtr<AsmScalar, 1> multidim_asm_;  ///< Assembly object
};


/**
 * Specialization defining FE vector operations
 */
class PatchFETestVector : public PatchFETestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmVector : public PatchFETestBase::AsmBase<dim> {
    public:
        /// Constructor
        AsmVector(PatchFETestBase *generic, uint quad_order, uint quad_diff_order)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order, quad_diff_order),
//          generic_inst_(generic),
          vector_shape_( this->bulk_integral_->vector_shape() ),
          vector_shape_side_( this->edge_integral_->vector_shape() ),
          grad_vector_shape_( this->bulk_integral_->grad_vector_shape() ),
          grad_vector_shape_side_( this->edge_integral_->grad_vector_shape() ),
          sym_grad_( this->bulk_integral_->vector_sym_grad() ),
          sym_grad_side_( this->edge_integral_->vector_sym_grad() ),
          divergence_( this->bulk_integral_->vector_divergence() ),
          divergence_side_( this->edge_integral_->vector_divergence() ),
          vector_join_( this->coupling_integral_->vector_join_shape() ),
          vector_join_grad_( this->coupling_integral_->gradient_vector_join_shape() )
        {}

        /// Destructor
        virtual ~AsmVector() {}


    	/** Declaration of data members **/
//        PatchFETestVector *generic_inst_;                                    ///< pointer to generic object
        FeQArray<Vector> vector_shape_;
        FeQArray<Vector> vector_shape_side_;
        FeQArray<Tensor> grad_vector_shape_;
        FeQArray<Tensor> grad_vector_shape_side_;
        FeQArray<Tensor> sym_grad_;
        FeQArray<Tensor> sym_grad_side_;
        FeQArray<Scalar> divergence_;
        FeQArray<Scalar> divergence_side_;
        FeQJoin<Vector> vector_join_;
        FeQJoin<Tensor> vector_join_grad_;
    };


	PatchFETestVector(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : PatchFETestBase(dh),
      multidim_asm_(this, quad_order, 0)
    {
		element_cache_map_.init(eval_points_);
    }

    ~PatchFETestVector() {}

    void set_integrals_arrays() {
        this->bulk_integrals_[0] = multidim_asm_[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm_[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm_[3_d]->bulk_integral_;
        this->edge_integrals_[0] = multidim_asm_[1_d]->edge_integral_;
        this->edge_integrals_[1] = multidim_asm_[2_d]->edge_integral_;
        this->edge_integrals_[2] = multidim_asm_[3_d]->edge_integral_;
        this->coupling_integrals_[0] = multidim_asm_[1_d]->coupling_integral_;
        this->coupling_integrals_[1] = multidim_asm_[2_d]->coupling_integral_;
        this->bulk_integrals_diff_order_[0] = multidim_asm_[1_d]->bulk_integral_diff_order_;
        this->bulk_integrals_diff_order_[1] = multidim_asm_[2_d]->bulk_integral_diff_order_;
        this->bulk_integrals_diff_order_[2] = multidim_asm_[3_d]->bulk_integral_diff_order_;
    }

    void initialize() {
        set_integrals_arrays();

        UpdateFlags u = update_values | update_inverse_jacobians | update_JxW_values | update_quadrature_points | update_volume_elements | update_gradients;
        UpdateFlags u_side = update_values | update_inverse_jacobians | update_side_JxW_values | update_normal_vectors | update_quadrature_points | update_gradients;
        fe_values_[0].initialize(*multidim_asm_[1_d]->quad_, *fe_[Dim<1>{}], u);
        fe_values_[1].initialize(*multidim_asm_[2_d]->quad_, *fe_[Dim<2>{}], u);
        fe_values_[2].initialize(*multidim_asm_[3_d]->quad_, *fe_[Dim<3>{}], u);
        fe_values_side_[0].initialize(*multidim_asm_[1_d]->quad_low_, *fe_[Dim<1>{}], u_side);
        fe_values_side_[1].initialize(*multidim_asm_[2_d]->quad_low_, *fe_[Dim<2>{}], u_side);
        fe_values_side_[2].initialize(*multidim_asm_[3_d]->quad_low_, *fe_[Dim<3>{}], u_side);
	    vec_view_1d_ = &fe_values_[0].vector_view(0);
	    vec_view_2d_ = &fe_values_[1].vector_view(0);
	    vec_view_3d_ = &fe_values_[2].vector_view(0);
	    vec_view_side_1d_ = &fe_values_side_[0].vector_view(0);
	    vec_view_side_2d_ = &fe_values_side_[1].vector_view(0);
	    vec_view_side_3d_ = &fe_values_side_[2].vector_view(0);
        this->patch_fe_values_.init_finalize();
    }

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(bool print_tables=false) {
        for(auto cell_it = dh_->local_range().begin(); cell_it != dh_->local_range().end(); ++cell_it) {
            auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, cell_it->dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(*cell_it, this->bulk_integrals_[cell_it->dim()-1]);
        	this->add_edge_integral(*cell_it, this->edge_integrals_[cell_it->dim()-1]);
//            this->add_coupling_integral(*cell_it, this->coupling_integrals_[cell_it->dim()-1]);
        	this->patch_fe_values_.make_permanent_ppv_data();
        }
        bulk_integral_data_.make_permanent();
        edge_integral_data_.make_permanent();
        coupling_integral_data_.make_permanent();
        element_cache_map_.make_paermanent_eval_points();
        element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        update_patch();

        if (print_tables) {
            std::stringstream ss;
            patch_fe_values_.print_operations(ss);
            WarningOut() << ss.str();
        }

        for(auto dh_cell : dh_->local_range() ) {
            ElementAccessor<3> elm = dh_cell.elm();
            auto p = *( bulk_integrals_[dh_cell.dim()-1]->points(element_cache_map_.position_in_cache(dh_cell.elm_idx()), &element_cache_map_).begin() );
            double jxw = 0.0, jxw_ref = 0.0;
            arma::vec3 vector_shape_dof0 = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_dof0_ref = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_dof1 = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_dof1_ref = {0.0, 0.0, 0.0};
            arma::mat33 grad_vector_dof0 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 grad_vector_dof0_ref = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 grad_vector_dof1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 grad_vector_dof1_ref = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 sym_grad_dof0 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 sym_grad_dof0_ref = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 sym_grad_dof1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 sym_grad_dof1_ref = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            double div_dof0 = 0.0, div_dof0_ref = 0.0;
            double div_dof1 = 0.0, div_dof1_ref = 0.0;
            switch (dh_cell.dim()) {
            case 1:
                fe_values_[0].reinit(elm);
                jxw = multidim_asm_[1_d]->jxw_(p);
                vector_shape_dof0 = multidim_asm_[1_d]->vector_shape_.shape(0)(p);
                grad_vector_dof0 = multidim_asm_[1_d]->grad_vector_shape_.shape(0)(p);
                sym_grad_dof0 = multidim_asm_[1_d]->sym_grad_.shape(0)(p);
                div_dof0 = multidim_asm_[1_d]->divergence_.shape(0)(p);
                jxw_ref = fe_values_[0].JxW(0);
                vector_shape_dof0_ref = vec_view_1d_->value(0, 0);
                grad_vector_dof0_ref = vec_view_1d_->grad(0, 0);
                sym_grad_dof0_ref = vec_view_1d_->sym_grad(0, 0);
                div_dof0_ref = vec_view_1d_->divergence(0, 0);
                break;
            case 2:
                fe_values_[1].reinit(elm);
                jxw = multidim_asm_[2_d]->jxw_(p);
                vector_shape_dof0 = multidim_asm_[2_d]->vector_shape_.shape(0)(p);
                vector_shape_dof1 = multidim_asm_[2_d]->vector_shape_.shape(1)(p);
                grad_vector_dof0 = multidim_asm_[2_d]->grad_vector_shape_.shape(0)(p);
                grad_vector_dof1 = multidim_asm_[2_d]->grad_vector_shape_.shape(1)(p);
                sym_grad_dof0 = multidim_asm_[2_d]->sym_grad_.shape(0)(p);
                sym_grad_dof1 = multidim_asm_[2_d]->sym_grad_.shape(1)(p);
                div_dof0 = multidim_asm_[2_d]->divergence_.shape(0)(p);
                div_dof1 = multidim_asm_[2_d]->divergence_.shape(1)(p);
                jxw_ref = fe_values_[1].JxW(0);
                vector_shape_dof0_ref = vec_view_2d_->value(0, 0);
                vector_shape_dof1_ref = vec_view_2d_->value(1, 0);
                grad_vector_dof0_ref = vec_view_2d_->grad(0, 0);
                grad_vector_dof1_ref = vec_view_2d_->grad(1, 0);
                sym_grad_dof0_ref = vec_view_2d_->sym_grad(0, 0);
                sym_grad_dof1_ref = vec_view_2d_->sym_grad(1, 0);
                div_dof0_ref = vec_view_2d_->divergence(0, 0);
                div_dof1_ref = vec_view_2d_->divergence(1, 0);
                break;
            case 3:
                fe_values_[2].reinit(elm);
                jxw = multidim_asm_[3_d]->jxw_(p);
                vector_shape_dof0 = multidim_asm_[3_d]->vector_shape_.shape(0)(p);
                vector_shape_dof1 = multidim_asm_[3_d]->vector_shape_.shape(1)(p);
                grad_vector_dof0 = multidim_asm_[3_d]->grad_vector_shape_.shape(0)(p);
                grad_vector_dof1 = multidim_asm_[3_d]->grad_vector_shape_.shape(1)(p);
                sym_grad_dof0 = multidim_asm_[3_d]->sym_grad_.shape(0)(p);
                sym_grad_dof1 = multidim_asm_[3_d]->sym_grad_.shape(1)(p);
                div_dof0 = multidim_asm_[3_d]->divergence_.shape(0)(p);
                div_dof1 = multidim_asm_[3_d]->divergence_.shape(1)(p);
                jxw_ref = fe_values_[2].JxW(0);
                vector_shape_dof0_ref = vec_view_3d_->value(0, 0);
                vector_shape_dof1_ref = vec_view_3d_->value(1, 0);
                grad_vector_dof0_ref = vec_view_3d_->grad(0, 0);
                grad_vector_dof1_ref = vec_view_3d_->grad(1, 0);
                sym_grad_dof0_ref = vec_view_3d_->sym_grad(0, 0);
                sym_grad_dof1_ref = vec_view_3d_->sym_grad(1, 0);
                div_dof0_ref = vec_view_3d_->divergence(0, 0);
                div_dof1_ref = vec_view_3d_->divergence(1, 0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, jxw_ref );
            EXPECT_ARMA_EQ( vector_shape_dof0, vector_shape_dof0_ref );
            EXPECT_ARMA_EQ( vector_shape_dof1, vector_shape_dof1_ref );
            EXPECT_ARMA_EQ( grad_vector_dof0, grad_vector_dof0_ref );
            EXPECT_ARMA_EQ( grad_vector_dof1, grad_vector_dof1_ref );
            EXPECT_ARMA_EQ( sym_grad_dof0, sym_grad_dof0_ref );
            EXPECT_ARMA_EQ( sym_grad_dof1, sym_grad_dof1_ref );
            EXPECT_DOUBLE_EQ( div_dof0, div_dof0_ref );
            EXPECT_DOUBLE_EQ( div_dof1, div_dof1_ref );
        }

        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            auto range = edge_integral_data_[i].edge_side_range;

            auto zero_edge_side = *range.begin();
            auto p = *( edge_integrals_[zero_edge_side.dim()-1]->points(zero_edge_side, &element_cache_map_).begin() );

            double jxw = 0.0, jxw_ref = 0.0;
            arma::vec3 vector_shape = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_ref = {0.0, 0.0, 0.0};
            arma::mat33 grad_vector = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 grad_vector_ref = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 sym_grad = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            arma::mat33 sym_grad_ref = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            double div = 0.0, div_ref = 0.0;
            switch (zero_edge_side.dim()) {
            case 1:
                jxw = multidim_asm_[1_d]->jxw_side_(p);
                vector_shape = multidim_asm_[1_d]->vector_shape_side_.shape(0)(p);
                grad_vector = multidim_asm_[1_d]->grad_vector_shape_side_.shape(0)(p);
                sym_grad = multidim_asm_[1_d]->sym_grad_side_.shape(0)(p);
                div = multidim_asm_[1_d]->divergence_side_.shape(0)(p);
                fe_values_side_[0].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[0].JxW(0);
                vector_shape_ref = vec_view_side_1d_->value(0, 0);
                grad_vector_ref = vec_view_side_1d_->grad(0, 0);
                sym_grad_ref = vec_view_side_1d_->sym_grad(0, 0);
                div_ref = vec_view_side_1d_->divergence(0, 0);
                break;
            case 2:
                jxw = multidim_asm_[2_d]->jxw_side_(p);
                vector_shape = multidim_asm_[2_d]->vector_shape_side_.shape(0)(p);
                grad_vector = multidim_asm_[2_d]->grad_vector_shape_side_.shape(0)(p);
                sym_grad = multidim_asm_[2_d]->sym_grad_side_.shape(0)(p);
                div = multidim_asm_[2_d]->divergence_side_.shape(0)(p);
                fe_values_side_[1].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[1].JxW(0);
                vector_shape_ref = vec_view_side_2d_->value(0, 0);
                grad_vector_ref = vec_view_side_2d_->grad(0, 0);
                sym_grad_ref = vec_view_side_2d_->sym_grad(0, 0);
                div_ref = vec_view_side_2d_->divergence(0, 0);
                break;
            case 3:
                jxw = multidim_asm_[3_d]->jxw_side_(p);
                vector_shape = multidim_asm_[3_d]->vector_shape_side_.shape(0)(p);
                grad_vector = multidim_asm_[3_d]->grad_vector_shape_side_.shape(0)(p);
                sym_grad = multidim_asm_[3_d]->sym_grad_side_.shape(0)(p);
                div = multidim_asm_[3_d]->divergence_side_.shape(0)(p);
                fe_values_side_[2].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[2].JxW(0);
                vector_shape_ref = vec_view_side_3d_->value(0, 0);
                grad_vector_ref = vec_view_side_3d_->grad(0, 0);
                sym_grad_ref = vec_view_side_3d_->sym_grad(0, 0);
                div_ref = vec_view_side_3d_->divergence(0, 0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, jxw_ref );
            EXPECT_ARMA_EQ( vector_shape, vector_shape_ref );
            EXPECT_ARMA_EQ( grad_vector, grad_vector_ref );
            EXPECT_ARMA_EQ( sym_grad, sym_grad_ref );
            EXPECT_DOUBLE_EQ( div, div_ref );
        }

        arma::vec3 arma_zero_vec = arma::zeros(3);
        arma::mat33 arma_zero_mat = arma::zeros(3,3);
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            DHCellAccessor cell_lower_dim = coupling_integral_data_[i].cell;
            DHCellSide neighb_side = coupling_integral_data_[i].side;;
            //std::cout << " el high " << neighb_side.elem_idx() << ", el low: " << cell_lower_dim.elm_idx() << std::endl;

            auto p_high = *( coupling_points(neighb_side.dim(), neighb_side).begin() );
            p_high.inc();
            auto p_low = p_high.lower_dim(cell_lower_dim);

            uint i_dof_high=0, i_dof_low=0;
            switch (neighb_side.dim()) {
            case 2:
                fe_values_[0].reinit(cell_lower_dim.elm());
                fe_values_side_[1].reinit(neighb_side.side());
                for (uint i_dof=0; i_dof<multidim_asm_[1_d]->vector_join_.n_dofs_both(); ++i_dof) {
                    if (multidim_asm_[1_d]->vector_join_.is_high_dim(i_dof)) {
                        auto result = multidim_asm_[1_d]->vector_join_.shape(i_dof)(p_high);
                        auto ref = vec_view_side_2d_->value(i_dof_high, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_low = multidim_asm_[1_d]->vector_join_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( result_low, arma_zero_vec );

                        auto grad_result = multidim_asm_[1_d]->vector_join_grad_.shape(i_dof)(p_high);
                        auto grad_ref = vec_view_side_2d_->grad(i_dof_high, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_low = multidim_asm_[1_d]->vector_join_grad_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( grad_result_low, arma_zero_mat );
                        i_dof_high++;
                    }
                    else {
                        auto result = multidim_asm_[1_d]->vector_join_.shape(i_dof)(p_low);
                        auto ref = vec_view_1d_->value(i_dof_low, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_high = multidim_asm_[1_d]->vector_join_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( result_high, arma_zero_vec );

                        auto grad_result = multidim_asm_[1_d]->vector_join_grad_.shape(i_dof)(p_low);
                        auto grad_ref = vec_view_1d_->grad(i_dof_low, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_high = multidim_asm_[1_d]->vector_join_grad_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( grad_result_high, arma_zero_mat );
                        i_dof_low++;
                    }
                }
                break;
            case 3:
                fe_values_[1].reinit(cell_lower_dim.elm());
                fe_values_side_[2].reinit(neighb_side.side());
                for (uint i_dof=0; i_dof<multidim_asm_[2_d]->vector_join_.n_dofs_both(); ++i_dof) {
                    if (multidim_asm_[2_d]->vector_join_.is_high_dim(i_dof)) {
                        auto result = multidim_asm_[2_d]->vector_join_.shape(i_dof)(p_high);
                        auto ref = vec_view_side_3d_->value(i_dof_high, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_low = multidim_asm_[2_d]->vector_join_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( result_low, arma_zero_vec );

                        auto grad_result = multidim_asm_[2_d]->vector_join_grad_.shape(i_dof)(p_high);
                        auto grad_ref = vec_view_side_3d_->grad(i_dof_high, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_low = multidim_asm_[2_d]->vector_join_grad_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( grad_result_low, arma_zero_mat );
                        i_dof_high++;
                    }
                    else {
                        auto result = multidim_asm_[2_d]->vector_join_.shape(i_dof)(p_low);
                        auto ref = vec_view_2d_->value(i_dof_low, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_high = multidim_asm_[2_d]->vector_join_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( result_high, arma_zero_vec );

                        auto grad_result = multidim_asm_[2_d]->vector_join_grad_.shape(i_dof)(p_low);
                        auto grad_ref = vec_view_2d_->grad(i_dof_low, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_high = multidim_asm_[2_d]->vector_join_grad_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( grad_result_high, arma_zero_mat );
                        i_dof_low++;
                    }
                }
                break;
            }
        }

    }

    MixedPtr<AsmVector, 1> multidim_asm_;  ///< Assembly object

    ///< Vector view in cell calculation.
    const FEValuesViews::Vector<3> * vec_view_1d_;
    const FEValuesViews::Vector<3> * vec_view_2d_;
    const FEValuesViews::Vector<3> * vec_view_3d_;
    ///< Vector view in SIDE calculation.
    const FEValuesViews::Vector<3> * vec_view_side_1d_;
    const FEValuesViews::Vector<3> * vec_view_side_2d_;
    const FEValuesViews::Vector<3> * vec_view_side_3d_;
};


/**
 * Specialization defining FE vector operations
 */
class PatchFETestQuadOrders : public PatchFETestBase {
public:
    /// Represent assembly class similar to assembly objects in equations
    template <unsigned int dim>
    class AsmQuadOrders : public PatchFETestBase::AsmBase<dim> {
    public:
        /// Constructor
        AsmQuadOrders(PatchFETestBase *generic, uint quad_order, uint quad_diff_order)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order, quad_diff_order),
          //generic_inst_(generic),
	      jxw_diff_order_( this->bulk_integral_diff_order_->JxW() ),
	      vector_shape_diff_order_( this->bulk_integral_diff_order_->vector_shape() )
        {}

        /// Destructor
        virtual ~AsmQuadOrders() {}


    	/** Declaration of data members **/
//        PatchFETestQuadOrders *generic_inst_;                                    ///< pointer to generic object
        FeQ<Scalar> jxw_diff_order_;
        FeQArray<Vector> vector_shape_diff_order_;
    };


    PatchFETestQuadOrders(unsigned int quad_order_1, unsigned int quad_order_2, std::shared_ptr<DOFHandlerMultiDim> dh_1, std::shared_ptr<DOFHandlerMultiDim> dh_2)
    : PatchFETestBase(dh_1), fe_values_diff_order_(3),
      multidim_asm_(this, quad_order_1, quad_order_2)
    {
        element_cache_map_.init(eval_points_);
    }

    void set_integrals_arrays() {
        this->bulk_integrals_[0] = multidim_asm_[1_d]->bulk_integral_;
        this->bulk_integrals_[1] = multidim_asm_[2_d]->bulk_integral_;
        this->bulk_integrals_[2] = multidim_asm_[3_d]->bulk_integral_;
        this->edge_integrals_[0] = multidim_asm_[1_d]->edge_integral_;
        this->edge_integrals_[1] = multidim_asm_[2_d]->edge_integral_;
        this->edge_integrals_[2] = multidim_asm_[3_d]->edge_integral_;
        this->coupling_integrals_[0] = multidim_asm_[1_d]->coupling_integral_;
        this->coupling_integrals_[1] = multidim_asm_[2_d]->coupling_integral_;
        this->bulk_integrals_diff_order_[0] = multidim_asm_[1_d]->bulk_integral_diff_order_;
        this->bulk_integrals_diff_order_[1] = multidim_asm_[2_d]->bulk_integral_diff_order_;
        this->bulk_integrals_diff_order_[2] = multidim_asm_[3_d]->bulk_integral_diff_order_;
    }

    void initialize() {
        set_integrals_arrays();

        UpdateFlags u = update_values | update_inverse_jacobians | update_JxW_values | update_quadrature_points | update_volume_elements | update_gradients;
        UpdateFlags u_side = update_values | update_inverse_jacobians | update_side_JxW_values | update_normal_vectors | update_quadrature_points | update_gradients;
        fe_values_[0].initialize(*multidim_asm_[1_d]->quad_, *fe_[Dim<1>{}], u);
        fe_values_[1].initialize(*multidim_asm_[2_d]->quad_, *fe_[Dim<2>{}], u);
        fe_values_[2].initialize(*multidim_asm_[3_d]->quad_, *fe_[Dim<3>{}], u);
        fe_values_side_[0].initialize(*multidim_asm_[1_d]->quad_low_, *fe_[Dim<1>{}], u_side);
        fe_values_side_[1].initialize(*multidim_asm_[2_d]->quad_low_, *fe_[Dim<2>{}], u_side);
        fe_values_side_[2].initialize(*multidim_asm_[3_d]->quad_low_, *fe_[Dim<3>{}], u_side);

        fe_values_diff_order_[0].initialize(*multidim_asm_[1_d]->quad_diff_order_, *fe_[Dim<1>{}], u);
        fe_values_diff_order_[1].initialize(*multidim_asm_[2_d]->quad_diff_order_, *fe_[Dim<2>{}], u);
        fe_values_diff_order_[2].initialize(*multidim_asm_[3_d]->quad_diff_order_, *fe_[Dim<3>{}], u);
        vec_view_1d_ = &fe_values_diff_order_[0].vector_view(0);
	    vec_view_2d_ = &fe_values_diff_order_[1].vector_view(0);
	    vec_view_3d_ = &fe_values_diff_order_[2].vector_view(0);

	    this->patch_fe_values_.init_finalize();
    }

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void update_patch() {
        patch_fe_values_.resize_tables(this->eval_points_);
        for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            uint dim = bulk_integral_data_[i].cell.dim();
            uint element_patch_idx = element_cache_map_.position_in_cache(bulk_integral_data_[i].cell.elm_idx());
            uint elm_pos = patch_fe_values_.register_element(bulk_integral_data_[i].cell, element_patch_idx);
            if ( bulk_integral_data_[i].subset_index == (unsigned int)(bulk_integrals_[dim-1]->get_subset_idx()) ) {
                uint i_point = 0;
                for (auto p : this->bulk_integrals_[dim-1]->points(element_patch_idx, &element_cache_map_)) {
                    patch_fe_values_.register_bulk_point(bulk_integral_data_[i].cell, elm_pos, p.value_cache_idx(), i_point++);
                }
            } else if ( bulk_integral_data_[i].subset_index == (unsigned int)(bulk_integrals_diff_order_[dim-1]->get_subset_idx()) ) {
                uint i_point = 0;
                for (auto p : this->bulk_integrals_diff_order_[dim-1]->points(element_patch_idx, &element_cache_map_)) {
                    patch_fe_values_.register_bulk_point(bulk_integral_data_[i].cell, elm_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            auto range = edge_integral_data_[i].edge_side_range;
            uint dim = range.begin()->dim();
            for( DHCellSide edge_side : range )
            {
                uint side_pos = patch_fe_values_.register_side(edge_side, &element_cache_map_);
                uint i_point = 0;
                for (auto p : this->edge_points(dim, edge_side) ) {
                    patch_fe_values_.register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
        this->reinit_patch_fe();
    }

    void test_evaluation(bool print_tables=false) {
        for(auto cell_it = dh_->local_range().begin(); cell_it != dh_->local_range().end(); ++cell_it) {
            auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, cell_it->dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(*cell_it, this->bulk_integrals_[cell_it->dim()-1]);
        	this->add_bulk_integral(*cell_it, this->bulk_integrals_diff_order_[cell_it->dim()-1]);
        	this->add_edge_integral(*cell_it, this->edge_integrals_[cell_it->dim()-1]);
        	this->patch_fe_values_.make_permanent_ppv_data();
        }
        bulk_integral_data_.make_permanent();
        edge_integral_data_.make_permanent();
        coupling_integral_data_.make_permanent();
        element_cache_map_.make_paermanent_eval_points();
        element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        update_patch();

        if (print_tables) {
            std::stringstream ss;
            patch_fe_values_.print_operations(ss);
            WarningOut() << ss.str();
        }

        for(auto dh_cell : dh_->local_range() ) {
            ElementAccessor<3> elm = dh_cell.elm();
            double jxw = 0.0, jxw_ref = 0.0;
            double det = 0.0, det_ref = 0.0;
            arma::vec3 vector_shape_dof0 = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_dof0_ref = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_dof1 = {0.0, 0.0, 0.0};
            arma::vec3 vector_shape_dof1_ref = {0.0, 0.0, 0.0};
            fe_values_[dh_cell.dim()-1].reinit(elm);
            fe_values_diff_order_[dh_cell.dim()-1].reinit(elm);

            {
                auto p = *( bulk_integrals_[dh_cell.dim()-1]->points(element_cache_map_.position_in_cache(dh_cell.elm_idx()), &element_cache_map_).begin() );
                switch (dh_cell.dim()) {
                case 1:
                    jxw = multidim_asm_[1_d]->jxw_(p);
                    det = multidim_asm_[1_d]->det_(p);
                    jxw_ref = fe_values_[0].JxW(0);
                    det_ref = fe_values_[0].determinant(0);
                    break;
                case 2:
                    jxw = multidim_asm_[2_d]->jxw_(p);
                    det = multidim_asm_[2_d]->det_(p);
                    jxw_ref = fe_values_[1].JxW(0);
                    det_ref = fe_values_[1].determinant(0);
                    break;
                case 3:
                    jxw = multidim_asm_[3_d]->jxw_(p);
                    det = multidim_asm_[3_d]->det_(p);
                    jxw_ref = fe_values_[2].JxW(0);
                    det_ref = fe_values_[2].determinant(0);
                    break;
                }
                EXPECT_DOUBLE_EQ( jxw, jxw_ref );
                EXPECT_DOUBLE_EQ( det, det_ref );
            }

            uint k=0;
            for (auto pt : bulk_integrals_diff_order_[dh_cell.dim()-1]->points(element_cache_map_.position_in_cache(dh_cell.elm_idx()), &element_cache_map_)) {
                switch (dh_cell.dim()) {
                case 1:
                    jxw = multidim_asm_[1_d]->jxw_diff_order_(pt);
                    vector_shape_dof0 = multidim_asm_[1_d]->vector_shape_diff_order_.shape(0)(pt);
                    jxw_ref = fe_values_diff_order_[0].JxW(k);
                    vector_shape_dof0_ref = vec_view_1d_->value(0, 0);
                    break;
                case 2:
                    jxw = multidim_asm_[2_d]->jxw_diff_order_(pt);
                    vector_shape_dof0 = multidim_asm_[2_d]->vector_shape_diff_order_.shape(0)(pt);
                    vector_shape_dof1 = multidim_asm_[2_d]->vector_shape_diff_order_.shape(1)(pt);
                    jxw_ref = fe_values_diff_order_[1].JxW(k);
                    vector_shape_dof0_ref = vec_view_2d_->value(0, 0);
                    vector_shape_dof1_ref = vec_view_2d_->value(1, 0);
                    break;
                case 3:
                    jxw = multidim_asm_[3_d]->jxw_diff_order_(pt);
                    vector_shape_dof0 = multidim_asm_[3_d]->vector_shape_diff_order_.shape(0)(pt);
                    vector_shape_dof1 = multidim_asm_[3_d]->vector_shape_diff_order_.shape(1)(pt);
                    jxw_ref = fe_values_diff_order_[2].JxW(k);
                    vector_shape_dof0_ref = vec_view_3d_->value(0, 0);
                    vector_shape_dof1_ref = vec_view_3d_->value(1, 0);
                    break;
                }
                EXPECT_DOUBLE_EQ( jxw, jxw_ref );
                EXPECT_ARMA_EQ( vector_shape_dof0, vector_shape_dof0_ref );
                EXPECT_ARMA_EQ( vector_shape_dof1, vector_shape_dof1_ref );
                ++k;
            }
        }

        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            auto range = edge_integral_data_[i].edge_side_range;

            auto zero_edge_side = *range.begin();
            auto p = *( edge_integrals_[zero_edge_side.dim()-1]->points(zero_edge_side, &element_cache_map_).begin() );

            double jxw = 0.0, jxw_ref = 0.0;
            arma::vec3 normal_vec = {0.0, 0.0, 0.0};
            arma::vec3 normal_vec_ref = {0.0, 0.0, 0.0};
            switch (zero_edge_side.dim()) {
            case 1:
                jxw = multidim_asm_[1_d]->jxw_side_(p);
                normal_vec = multidim_asm_[1_d]->normal_vec_(p);
                fe_values_side_[0].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[0].JxW(0);
                normal_vec_ref = fe_values_side_[0].normal_vector(0);
                break;
            case 2:
                jxw = multidim_asm_[2_d]->jxw_side_(p);
                normal_vec = multidim_asm_[2_d]->normal_vec_(p);
                fe_values_side_[1].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[1].JxW(0);
                normal_vec_ref = fe_values_side_[1].normal_vector(0);
                break;
            case 3:
                jxw = multidim_asm_[3_d]->jxw_side_(p);
                normal_vec = multidim_asm_[3_d]->normal_vec_(p);
                fe_values_side_[2].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[2].JxW(0);
                normal_vec_ref = fe_values_side_[2].normal_vector(0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, jxw_ref );
            EXPECT_ARMA_EQ( normal_vec, normal_vec_ref );
        }
    }

    MixedPtr<AsmQuadOrders, 1> multidim_asm_;  ///< Assembly object

    std::vector<FEValues<3>> fe_values_diff_order_;                           ///< FeValues object of elements of dim 1,2,3

    ///< Vector view in cell calculation.
    const FEValuesViews::Vector<3> * vec_view_1d_;
    const FEValuesViews::Vector<3> * vec_view_2d_;
    const FEValuesViews::Vector<3> * vec_view_3d_;
};


/**
 * Used in speed_comparation test
 */
//class PatchFETestCompare : public PatchFETestBase {
//public:
//	/// Number of repeats
//	static const unsigned int n_loops = 1e05;
//
//	PatchFETestCompare(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
//    : PatchFETestBase(quad_order, dh)
//    {}
//
//    ~PatchFETestCompare() {}
//
//    void reinit_patch_fe() override {
//        START_TIMER("reinit_patch_fe");
//        for (uint i=0; i<PatchFETestCompare::n_loops; ++i)
//            patch_fe_values_.reinit_patch();
//        END_TIMER("reinit_patch_fe");
//    }
//
//    void reinit_fe_values() {
//        START_TIMER("reinit_fe_values");
//        for (uint i=0; i<PatchFETestCompare::n_loops; ++i) {
//
//            for(auto cell_it = this->dh_->local_range().begin(); cell_it != this->dh_->local_range().end(); ++cell_it) {
//                ElementAccessor<3> elm = cell_it->elm();
//                switch (elm.dim()) {
//                case 1:
//                    fe_values_[0].reinit(elm);
//                    for (uint si=0; si<2; ++si)
//                        fe_values_side_[0].reinit(*elm.side(si));
//                    break;
//                case 2:
//                    fe_values_[1].reinit(elm);
//                    for (uint si=0; si<3; ++si)
//                        fe_values_side_[1].reinit(*elm.side(si));
//                    break;
//                case 3:
//                    fe_values_[2].reinit(elm);
//                    for (uint si=0; si<4; ++si)
//                        fe_values_side_[2].reinit(*elm.side(si));
//                    break;
//                }
//            }
//
//        }
//        END_TIMER("reinit_fe_values");
//    }
//};


/// Complete test with FE scalar operations
void compare_evaluation_func_scalar(Mesh* mesh, unsigned int i_run, bool print_fa_data = false) {
    std::vector<uint> quad_orders = {1, 2};
    MixedPtr<FE_P_disc> fe(quad_orders[i_run]);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    PatchFETestScalar patch_fe(quad_orders[i_run], dh);
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

    PatchFETestVector patch_fe(quad_order, dh);
    patch_fe.initialize();
    patch_fe.test_evaluation(print_fa_data);
    patch_fe.reset();
    patch_fe.test_evaluation();
}


/// Complete test with FE vector operations
void compare_evaluation_diff_orders(Mesh* mesh, unsigned int quad_order_1, unsigned int quad_order_2, bool print_fa_data = false) {
    MixedPtr<FE_P> fe_p( quad_order_1 );
    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    PatchFETestQuadOrders patch_fe(quad_order_1, quad_order_2, dh, dh); // fix passing of dh
    patch_fe.initialize();
    patch_fe.test_evaluation(print_fa_data);
    patch_fe.reset();
    patch_fe.test_evaluation();
}



TEST(PatchFeTest, complete_evaluation) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    std::string input_str = "{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }";
    Mesh* mesh = mesh_full_constructor(input_str);

    // two tests with different quad_order and Scalar / Vector FE operations
    compare_evaluation_func_scalar(mesh, 0, true);
    compare_evaluation_func_scalar(mesh, 1);
    compare_evaluation_func_vector(mesh, 1, true);
    std::cout << " - Different quad orders test -----------------------------" << std::endl;
    compare_evaluation_diff_orders(mesh, 0, 2, true);
}


TEST(PatchFeTest, arena_alloc) {
    AssemblyArena asm_arena(1024*1024, 256);

    char *c = (char *)asm_arena.allocate_8<char>(4);
    c[0] = 'a';
    ArenaVec<double> vec1(10, asm_arena);
    ArenaVec<double> vec2(16, asm_arena);
    ArenaVec<double> vec3(12, asm_arena);
    PatchArena *patch_arena = asm_arena.get_child_arena();

    ArenaVec<double> vec4(10, *patch_arena);

    std::cout << "Sizes: " << vec1.data_size() << ", " << vec2.data_size() << ", " << vec3.data_size() << ", " << vec4.data_size() << ", " << c[0] << std::endl;

    patch_arena->reset();
    std::cout << "Reset: " << vec1.data_size() << ", " << vec2.data_size() << ", " << vec3.data_size() << ", " << c[0] << std::endl;
}
