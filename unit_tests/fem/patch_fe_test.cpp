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
#include "system/sys_profiler.hh"


// Define EXPECT_<...>_NEAR with fixed abs_error 1e-5
#define EXPECT_TEST_NEAR( A, B )\
  EXPECT_NEAR(A, B, 1e-4)

#define EXPECT_TEST_ARMA_NEAR( A, B ) \
    EXPECT_ARMA_NEAR(A, B, 1e-4);



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
      fe_(dh_->ds()->fe()),
	  eval_points_( std::make_shared<EvalPoints>() ),
      bulk_integral_data_(20, 10),
      edge_integral_data_(12, 6),
      coupling_integral_data_(12, 6)
    {
        used_element_idx_ = {0, 1, 2, 3, 8}; // dimension of used elements: 1D, 2D, 2D, 3D, 3D

        ref_bulk_jxw_ = {
            { 1.73205, 0.94281, 0.94281, 0.33333, 0.33333 },
            { 0.96225, 0.63182, 0.63182, 0.09799, 0.09799 }
        };
        ref_bulk_det_ = { 3.46410, 5.65685, 5.65685, 8.00000, 8.00000 };
        ref_bulk_scalar_shape_dof0_ = {
            { 0.78868, 0.66667, 0.66667, 0.58541, 0.58541 },
            { 0.68730, -0.08473, -0.08473, 0.32018, 0.32018 }
        };
        ref_bulk_scalar_shape_dof1_ = {
            { 0.21132, 0.16667, 0.16667, 0.13820, 0.13820 },
            { 0.40000, 0.19283, 0.19283, 0.26774, 0.26774 }
        };
        ref_bulk_grad_scalar_dof0_ = {
            { {-0.16667, -0.16667, 0.16667}, {0.00000, 0.00000, 0.50000}, {0.25000, 0.25000, 0.50000}, {0.50000, 0.00000, 0.50000}, {0.00000, 0.50000, 0.50000} },
            { {-0.42487, -0.42487, 0.42487}, {0.00000, 0.00000, -0.28379}, {-0.14190, -0.14190, -0.28379}, {0.94359, 0.00000, 0.94359}, {0.00000, 0.94359, 0.94359} }
        };
        ref_bulk_grad_scalar_dof1_ = {
            { {0.16667, 0.16667, -0.16667}, {-0.25000, -0.25000, -0.50000}, {-0.25000, -0.25000, 0.00000}, {-0.50000, 0.50000, -0.00000}, {-0.50000, -0.00000, -0.00000} },
            { {0.51640, 0.51640, -0.51640}, {-0.10810, -0.10810, 0.67569}, {0.33785, 0.33785, 0.89190}, {-1.25812, 1.44359, 0.18547}, {-1.44359, 0.18547, 0.18547} }
        };
        ref_bulk_vector_shape_dof0_ = {
            {0.78868, 0.00000, 0.00000}, {0.66667, 0.00000, 0.00000}, {0.66667, 0.00000, 0.00000}, {0.58541, 0.00000, 0.00000}, {0.58541, 0.00000, 0.00000}
        };
        ref_bulk_vector_shape_dof1_ = {
            {0.21132, 0.00000, 0.00000}, {0.16667, 0.00000, 0.00000}, {0.16667, 0.00000, 0.00000}, {0.13820, 0.00000, 0.00000}, {0.13820, 0.00000, 0.00000}
        };
        ref_bulk_grad_vector_dof0_ = {
            { -0.16667, 0.00000, 0.00000, -0.16667, 0.00000, 0.00000, 0.16667, 0.00000, 0.00000 },
            { 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.50000, 0.00000, 0.00000 },
            { 0.25000, 0.00000, 0.00000, 0.25000, 0.00000, 0.00000, 0.50000, 0.00000, 0.00000 },
            { 0.50000, -0.00000, -0.00000, 0.00000, 0.00000, 0.00000, 0.50000, -0.00000, -0.00000 },
            { 0.00000, 0.00000, 0.00000, 0.50000, -0.00000, -0.00000, 0.50000, -0.00000, -0.00000 }
        };
        ref_bulk_grad_vector_dof1_ = {
            { 0.16667, 0.00000, 0.00000, 0.16667, 0.00000, 0.00000, -0.16667, 0.00000, 0.00000 },
            { -0.25000, 0.00000, 0.00000, -0.25000, 0.00000, 0.00000, -0.50000, 0.00000, 0.00000 },
            { -0.25000, 0.00000, 0.00000, -0.25000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 },
            { -0.50000, -0.00000, -0.00000, 0.50000, 0.00000, 0.00000, -0.00000, -0.00000, -0.00000 },
            { -0.50000, 0.00000, 0.00000, -0.00000, -0.00000, -0.00000, -0.00000, -0.00000, -0.00000 }
        };
        ref_bulk_sym_grad_dof0_ = {
            { -0.16667, -0.08333, 0.08333, -0.08333, 0.00000, 0.00000, 0.08333, 0.00000, 0.00000 },
            { 0.00000, 0.00000, 0.25000, 0.00000, 0.00000, 0.00000, 0.25000, 0.00000, 0.00000 },
            { 0.25000, 0.12500, 0.25000, 0.12500, 0.00000, 0.00000, 0.25000, 0.00000, 0.00000 },
            { 0.50000, 0.00000, 0.25000, 0.00000, 0.00000, 0.00000, 0.25000, 0.00000, -0.00000 },
            { 0.00000, 0.25000, 0.25000, 0.25000, -0.00000, -0.00000, 0.25000, -0.00000, -0.00000 }
        };
        ref_bulk_sym_grad_dof1_ = {
            { 0.16667, 0.08333, -0.08333, 0.08333, 0.00000, 0.00000, -0.08333, 0.00000, 0.00000 },
            { -0.25000, -0.12500, -0.25000, -0.12500, 0.00000, 0.00000, -0.25000, 0.00000, 0.00000 },
            { -0.25000, -0.12500, 0.00000, -0.12500, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000 },
            { -0.50000, 0.25000, -0.00000, 0.25000, 0.00000, 0.00000, -0.00000, 0.00000, -0.00000 },
            { -0.50000, 0.00000, 0.00000, 0.00000, -0.00000, -0.00000, 0.00000, -0.00000, -0.00000 }
        };
        ref_bulk_div_dof0_ = { -0.16667, 0.00000, 0.25000, 0.50000, 0.00000 };
        ref_bulk_div_dof1_ = { 0.16667, -0.25000, -0.25000, -0.50000, -0.50000 };
        ref_side_jxw_ = {
            { 0.94281, 0.94281 },
            { 0.63182, 0.63182 }
        };
        ref_normal_vec_ = { {-0.70711, 0.00000, -0.70711}, {0.70711, 0.00000, 0.70711} };
        ref_side_scalar_shape_ = {
            { -0.00000, 0.66667 },
            { -0.00000, -0.08473 }
        };
        ref_side_grad_scalar_ = {
            { {0.50000, 0.00000, 0.50000}, {0.00000, 0.50000, 0.50000} },
            { {-0.50000, 0.00000, -0.50000}, {0.00000, -0.28379, -0.28379} }
        };
        ref_side_vector_shape_ = { {-0.00000, 0.00000, 0.00000}, {0.66667, 0.00000, 0.00000} };
        ref_side_grad_vector_ = { { 0.50000, -0.00000, -0.00000, 0.00000, 0.00000, 0.00000, 0.50000, -0.00000, -0.00000 },
                                  { 0.00000, 0.00000, 0.00000, 0.50000, 0.00000, 0.00000, 0.50000, 0.00000, 0.00000 } };
        ref_side_sym_grad_ = { { 0.50000, 0.00000, 0.25000, 0.00000, 0.00000, 0.00000, 0.25000, 0.00000, -0.00000 },
                               { 0.00000, 0.25000, 0.25000, 0.25000, 0.00000, 0.00000, 0.25000, 0.00000, 0.00000 } };
        ref_side_vector_div_ = { 0.50000, 0.00000 };
        ref_join_scalar_shape_ = { { 0.21132, 0.21132, 0.00000, 0.16667, 0.16667, 0.16667 },
                                   { -0.08730, -0.08730, 0.00000, -0.04821, -0.04821, -0.04821 } };
        ref_join_vector_shape_ = { {0.00000, 0.00000, 0.78868}, {0.00000, 0.00000, 0.78868}, {0.00000, 0.00000, 0.00000},
                                   {0.00000, 0.00000, 0.66667}, {0.00000, 0.00000, 0.66667}, {0.00000, 0.00000, 0.66667} };
        ref_join_grad_vector_ = {
                { 0.00000, 0.00000, 0.25000, 0.00000, 0.00000, 0.25000, 0.00000, 0.00000, 0.00000 },
                { 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, -0.50000 },
                { 0.00000, 0.00000, -0.50000, 0.00000, 0.00000, 0.50000, -0.00000, -0.00000, -0.00000 },
                { 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.50000, 0.00000, 0.00000, 0.00000 },
                { -0.00000, -0.00000, -0.00000, 0.00000, 0.00000, 0.00000, -0.00000, -0.00000, -0.50000 },
                { 0.00000, 0.00000, 0.00000, -0.00000, -0.00000, -0.00000, -0.00000, -0.00000, -0.50000 }
        };
    }

    ~PatchFETestBase() {}

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

    void add_coupling_integral(DHCellAccessor cell, std::shared_ptr<CouplingIntegral> coupling_integral) {
        bool add_low = true;
        uint dim = cell.dim();
        for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
            if (cell.dim() != neighb_side.dim()-1) continue;
            coupling_integral_data_.emplace_back(cell, coupling_integral->get_subset_low_idx(), neighb_side,
                    coupling_integral->get_subset_high_idx());
            auto &ppv_high = patch_fe_values_.ppv(side_domain, dim+1);
            ++ppv_high.n_mesh_items_;

            unsigned int reg_idx_low = cell.elm().region_idx().idx();
            unsigned int reg_idx_high = neighb_side.element().region_idx().idx();
            for (auto p : coupling_integral->points(neighb_side, &element_cache_map_) ) {
                element_cache_map_.add_eval_point(reg_idx_high, neighb_side.elem_idx(), p.eval_point_idx(), neighb_side.cell().local_idx());

                if (add_low) {
                    auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                    element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
                }
            }
            add_low = false;
        }
    }

    void resize_tables() {
        for (uint i=1; i<=3; ++i) {
            patch_fe_values_.ppv(bulk_domain, i).resize_tables(eval_points_->get_max_bulk_quad_size(i), patch_fe_values_.patch_arena());
            patch_fe_values_.ppv(side_domain, i).resize_tables(eval_points_->get_max_side_quad_size(i), patch_fe_values_.patch_arena());
        }
        patch_fe_values_.clean_elements_map();
    }

    void update_patch() {
        this->resize_tables();
        for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            uint dim = bulk_integral_data_[i].cell.dim();
            uint element_patch_idx = element_cache_map_.position_in_cache(bulk_integral_data_[i].cell.elm_idx());
            uint elm_pos = patch_fe_values_.register_element(bulk_integral_data_[i].cell, element_patch_idx);
            uint i_point = 0;
            for (auto p : bulk_integrals_[dim-1]->points(element_patch_idx, &element_cache_map_) ) {
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
                for (auto p : edge_integrals_[dim-1]->points(edge_side, &element_cache_map_) ) {
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
            for (auto p_high : coupling_integrals_[dim-2]->points(coupling_integral_data_[i].side, &element_cache_map_) )
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

    std::shared_ptr<EvalPoints> eval_points_;                                 ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                                       ///< ElementCacheMap according to EvalPoints
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_;             ///< Bulk integrals of dim 1,2,3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_integrals_;             ///< Edge integrals of dim 1,2,3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_integrals_;     ///< Coupling integrals of dim 1-2,2-3
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_diff_order_;  ///< Bulk integrals of dim 1,2,3 of high order
    RevertableList<BulkIntegralData> bulk_integral_data_;                     ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData> edge_integral_data_;                     ///< Holds data for computing edge integrals.
    RevertableList<CouplingIntegralData> coupling_integral_data_;             ///< Holds data for computing couplings integrals.

    std::vector<unsigned int> used_element_idx_;                              ///< List of mesh idx of elements used in tests

    /* Reference values */
    std::vector< std::vector<double> > ref_bulk_jxw_;
    std::vector<double> ref_bulk_det_;
    std::vector< std::vector<double> > ref_bulk_scalar_shape_dof0_;
    std::vector< std::vector<double> > ref_bulk_scalar_shape_dof1_;
    std::vector< std::vector<arma::vec3> > ref_bulk_grad_scalar_dof0_;
    std::vector< std::vector<arma::vec3> > ref_bulk_grad_scalar_dof1_;
    std::vector<arma::vec3> ref_bulk_vector_shape_dof0_;
    std::vector<arma::vec3> ref_bulk_vector_shape_dof1_;
    std::vector<arma::mat33> ref_bulk_grad_vector_dof0_;
    std::vector<arma::mat33> ref_bulk_grad_vector_dof1_;
    std::vector<arma::mat33> ref_bulk_sym_grad_dof0_;
    std::vector<arma::mat33> ref_bulk_sym_grad_dof1_;
    std::vector<double> ref_bulk_div_dof0_;
    std::vector<double> ref_bulk_div_dof1_;
    std::vector< std::vector<double> > ref_side_jxw_;
    std::vector<arma::vec3> ref_normal_vec_;
    std::vector< std::vector<double> > ref_side_scalar_shape_;
    std::vector< std::vector<arma::vec3> > ref_side_grad_scalar_;
    std::vector<arma::vec3> ref_side_vector_shape_;
    std::vector<arma::mat33> ref_side_grad_vector_;
    std::vector<arma::mat33> ref_side_sym_grad_;
    std::vector<double> ref_side_vector_div_;
    std::vector< std::vector<double> > ref_join_scalar_shape_;
    std::vector<arma::vec3> ref_join_vector_shape_;
    std::vector<arma::mat33> ref_join_grad_vector_;
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
        AsmScalar(PatchFETestScalar *generic, uint quad_order, uint quad_diff_order)
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

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int i_run, unsigned int i_test_elem) {
            auto p = *( this->bulk_integral_->points(this->generic_->element_cache_map_.position_in_cache(dh_cell.elm_idx())).begin() );
            double jxw = this->jxw_(p);
            double det = this->det_(p);
            double scalar_shape_dof0 = scalar_shape_.shape(0)(p);
            double scalar_shape_dof1 = scalar_shape_.shape(1)(p);
            arma::vec3 grad_scalar_dof0 = grad_scalar_shape_.shape(0)(p);
            arma::vec3 grad_scalar_dof1 = grad_scalar_shape_.shape(1)(p);

            EXPECT_TEST_NEAR( jxw, this->generic_->ref_bulk_jxw_[i_run][i_test_elem] );
            EXPECT_TEST_NEAR( det, this->generic_->ref_bulk_det_[i_test_elem] );
            EXPECT_TEST_NEAR( scalar_shape_dof0, this->generic_->ref_bulk_scalar_shape_dof0_[i_run][i_test_elem] );
            EXPECT_TEST_NEAR( scalar_shape_dof1, this->generic_->ref_bulk_scalar_shape_dof1_[i_run][i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( grad_scalar_dof0, this->generic_->ref_bulk_grad_scalar_dof0_[i_run][i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( grad_scalar_dof1, this->generic_->ref_bulk_grad_scalar_dof1_[i_run][i_test_elem] );
        }

        void test_side_values(RangeConvert<DHEdgeSide, DHCellSide> range, unsigned int i_run) {
            uint k=0;
            for (DHCellSide cell_side : range) {
                auto p = *( this->edge_integral_->points(cell_side).begin() );

                double jxw = this->jxw_side_(p);
                arma::vec3 normal_vec = this->normal_vec_(p);
                double scalar_shape = scalar_shape_side_.shape(0)(p);
                arma::vec3 grad_scalar = grad_scalar_shape_side_.shape(0)(p);

                EXPECT_TEST_NEAR( jxw, this->generic_->ref_side_jxw_[i_run][k] );
                EXPECT_TEST_ARMA_NEAR( normal_vec, this->generic_->ref_normal_vec_[k] );
                EXPECT_TEST_NEAR( scalar_shape, this->generic_->ref_side_scalar_shape_[i_run][k] );
                EXPECT_TEST_ARMA_NEAR( grad_scalar, this->generic_->ref_side_grad_scalar_[i_run][k] );
                ++k;
            }
        }

        void test_join_values(DHCellAccessor cell_lower_dim, DHCellSide neighb_side, unsigned int i_run, unsigned int i_test_join) {
            auto p_high = *( this->coupling_integral_->points(neighb_side).begin() );
            auto p_low = p_high.lower_dim(cell_lower_dim);

            double result = 0.0, result_zero = 0.0;
            for (uint i_dof=0; i_dof<conc_join_shape_.n_dofs_both(); ++i_dof) {
                if (conc_join_shape_.is_high_dim(i_dof)) {
                    result = conc_join_shape_.shape(i_dof)(p_high);
                    result_zero = conc_join_shape_.shape(i_dof)(p_low);
                }
                else {
                    result = conc_join_shape_.shape(i_dof)(p_low);
                    result_zero = conc_join_shape_.shape(i_dof)(p_high);
                }
            }
            EXPECT_TEST_NEAR( result, this->generic_->ref_join_scalar_shape_[i_run][i_test_join] );
            EXPECT_TEST_NEAR( result_zero, 0.0 );
        }


    	/** Declaration of data members **/
        //PatchFETestScalar *generic_inst_;                                    ///< pointer to generic object
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
        this->patch_fe_values_.init_finalize();
    }

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(unsigned int i_run, bool print_tables=false) {
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
        	this->add_edge_integral(dh_cell, this->edge_integrals_[dh_cell.dim()-1]);
            this->add_coupling_integral(dh_cell, this->coupling_integrals_[dh_cell.dim()-1]);
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

        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            RangeConvert<DHEdgeSide, DHCellSide> range = edge_integral_data_[i].edge_side_range;
            switch (range.begin()->dim()) {
            case 1:
                multidim_asm_[1_d]->test_side_values(range, i_run);
                break;
            case 2:
                multidim_asm_[2_d]->test_side_values(range, i_run);
                break;
            case 3:
                multidim_asm_[3_d]->test_side_values(range, i_run);
                break;
            }
        }

        unsigned int i_test_join=0;
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            DHCellAccessor cell_lower_dim = coupling_integral_data_[i].cell;
            DHCellSide neighb_side = coupling_integral_data_[i].side;;
            switch (cell_lower_dim.dim()) {
            case 1:
                multidim_asm_[1_d]->test_join_values(cell_lower_dim, neighb_side, i_run, i_test_join);
                break;
            case 2:
                multidim_asm_[2_d]->test_join_values(cell_lower_dim, neighb_side, i_run, i_test_join);
                break;
            }
            ++i_test_join;
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

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int i_test_elem) {
            auto p = *( this->bulk_integral_->points(this->generic_->element_cache_map_.position_in_cache(dh_cell.elm_idx())).begin() );
            double jxw = this->jxw_(p);
            arma::vec3 vector_shape_dof0 = vector_shape_.shape(0)(p);
            arma::vec3 vector_shape_dof1 = vector_shape_.shape(1)(p);
            arma::mat33 grad_vector_dof0 = grad_vector_shape_.shape(0)(p);
            arma::mat33 grad_vector_dof1 = grad_vector_shape_.shape(1)(p);
            arma::mat33 sym_grad_dof0 = sym_grad_.shape(0)(p);
            arma::mat33 sym_grad_dof1 = sym_grad_.shape(1)(p);
            double div_dof0 = divergence_.shape(0)(p);
            double div_dof1 = divergence_.shape(1)(p);

            EXPECT_TEST_NEAR( jxw, this->generic_->ref_bulk_jxw_[0][i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( vector_shape_dof0, this->generic_->ref_bulk_vector_shape_dof0_[i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( vector_shape_dof1, this->generic_->ref_bulk_vector_shape_dof1_[i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( grad_vector_dof0, this->generic_->ref_bulk_grad_vector_dof0_[i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( grad_vector_dof1, this->generic_->ref_bulk_grad_vector_dof1_[i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( sym_grad_dof0, this->generic_->ref_bulk_sym_grad_dof0_[i_test_elem] );
            EXPECT_TEST_ARMA_NEAR( sym_grad_dof1, this->generic_->ref_bulk_sym_grad_dof1_[i_test_elem] );
            EXPECT_TEST_NEAR( div_dof0, this->generic_->ref_bulk_div_dof0_[i_test_elem] );
            EXPECT_TEST_NEAR( div_dof1, this->generic_->ref_bulk_div_dof1_[i_test_elem] );
        }

        void test_side_values(RangeConvert<DHEdgeSide, DHCellSide> range) {
            uint k=0;
            for (DHCellSide cell_side : range) {
                auto p = *( this->edge_integral_->points(cell_side).begin() );

                double jxw = this->jxw_side_(p);
                arma::vec3 vector_shape = vector_shape_side_.shape(0)(p);
                arma::mat33 grad_vector = grad_vector_shape_side_.shape(0)(p);
                arma::mat33 sym_grad = sym_grad_side_.shape(0)(p);
                double div = divergence_side_.shape(0)(p);
                EXPECT_TEST_NEAR( jxw, this->generic_->ref_side_jxw_[0][k] );
                EXPECT_TEST_ARMA_NEAR( vector_shape, this->generic_->ref_side_vector_shape_[k] );
                EXPECT_TEST_ARMA_NEAR( grad_vector, this->generic_->ref_side_grad_vector_[k] );
                EXPECT_TEST_ARMA_NEAR( sym_grad, this->generic_->ref_side_sym_grad_[k] );
                EXPECT_TEST_NEAR( div, this->generic_->ref_side_vector_div_[k] );
                ++k;
            }
        }

        void test_join_values(DHCellAccessor cell_lower_dim, DHCellSide neighb_side, unsigned int i_test_join) {
            auto p_high = *( this->coupling_integral_->points(neighb_side).begin() );
            p_high.inc();
            auto p_low = p_high.lower_dim(cell_lower_dim);

            arma::vec3 arma_zero_vec = arma::zeros(3);
            arma::mat33 arma_zero_mat = arma::zeros(3,3);
            arma::vec3 shape_result = arma::zeros(3);
            arma::vec3 shape_zero = arma::zeros(3);
            arma::mat33 grad_result = arma::zeros(3,3);
            arma::mat33 grad_zero = arma::zeros(3,3);
            for (uint i_dof=0; i_dof<vector_join_.n_dofs_both(); ++i_dof) {
                if (vector_join_.is_high_dim(i_dof)) {
                    shape_result = vector_join_.shape(i_dof)(p_high);
                    shape_zero = vector_join_.shape(i_dof)(p_low);

                    grad_result = vector_join_grad_.shape(i_dof)(p_high);
                    grad_zero = vector_join_grad_.shape(i_dof)(p_low);
                }
                else {
                    shape_result = vector_join_.shape(i_dof)(p_low);
                    shape_zero = vector_join_.shape(i_dof)(p_high);

                    grad_result = vector_join_grad_.shape(i_dof)(p_low);
                    grad_zero = vector_join_grad_.shape(i_dof)(p_high);
                }
            }
            EXPECT_TEST_ARMA_NEAR( shape_result, this->generic_->ref_join_vector_shape_[i_test_join] );
            EXPECT_TEST_ARMA_NEAR( shape_zero, arma_zero_vec );
            EXPECT_TEST_ARMA_NEAR( grad_result, this->generic_->ref_join_grad_vector_[i_test_join] );
            EXPECT_TEST_ARMA_NEAR( grad_zero, arma_zero_mat );
        }

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
        this->patch_fe_values_.init_finalize();
    }

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(bool print_tables=false) {
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
        	this->add_edge_integral(dh_cell, this->edge_integrals_[dh_cell.dim()-1]);
            this->add_coupling_integral(dh_cell, this->coupling_integrals_[dh_cell.dim()-1]);
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

        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            RangeConvert<DHEdgeSide, DHCellSide> range = edge_integral_data_[i].edge_side_range;
            switch (range.begin()->dim()) {
            case 1:
                multidim_asm_[1_d]->test_side_values(range);
                break;
            case 2:
                multidim_asm_[2_d]->test_side_values(range);
                break;
            case 3:
                multidim_asm_[3_d]->test_side_values(range);
                break;
            }
        }

        unsigned int i_test_join=0;
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            DHCellAccessor cell_lower_dim = coupling_integral_data_[i].cell;
            DHCellSide neighb_side = coupling_integral_data_[i].side;;
            switch (cell_lower_dim.dim()) {
            case 1:
                multidim_asm_[1_d]->test_join_values(cell_lower_dim, neighb_side, i_test_join);
                break;
            case 2:
                multidim_asm_[2_d]->test_join_values(cell_lower_dim, neighb_side, i_test_join);
                break;
            }
            ++i_test_join;
        }

    }

    MixedPtr<AsmVector, 1> multidim_asm_;  ///< Assembly object
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
        AsmQuadOrders(PatchFETestQuadOrders *generic, uint quad_order, uint quad_diff_order)
        : PatchFETestBase::AsmBase<dim>(generic, quad_order, quad_diff_order),
          generic_inst_(generic),
	      jxw_diff_order_( this->bulk_integral_diff_order_->JxW() ),
	      vector_shape_diff_order_( this->bulk_integral_diff_order_->vector_shape() )
        {}

        /// Destructor
        virtual ~AsmQuadOrders() {}

        void test_bulk_values(DHCellAccessor dh_cell, unsigned int i_test_elem) {
            // low order quadrature
            {
                auto p = *( this->bulk_integral_->points(generic_inst_->element_cache_map_.position_in_cache(dh_cell.elm_idx())).begin() );
                double jxw = this->jxw_(p);
                double det = this->det_(p);
                EXPECT_TEST_NEAR( jxw, generic_inst_->ref_bulk_jxw_low_order_[i_test_elem] );
                EXPECT_TEST_NEAR( det, generic_inst_->ref_bulk_det_[i_test_elem] );
            }

            // high order quadrature
            uint k=0;
            for (auto pt : this->bulk_integral_diff_order_->points(generic_inst_->element_cache_map_.position_in_cache(dh_cell.elm_idx()))) {
                double jxw = jxw_diff_order_(pt);
                arma::vec3 vector_shape_dof0 = vector_shape_diff_order_.shape(0)(pt);
                arma::vec3 vector_shape_dof1 = vector_shape_diff_order_.shape(1)(pt);
                EXPECT_TEST_NEAR( jxw, generic_inst_->ref_bulk_jxw_high_order_[i_test_elem][k] );
                EXPECT_TEST_ARMA_NEAR( vector_shape_dof0, generic_inst_->ref_vector_shape_dof0_[i_test_elem][k] );
                EXPECT_TEST_ARMA_NEAR( vector_shape_dof1, generic_inst_->ref_vector_shape_dof1_[i_test_elem][k] );
                ++k;
            }
        }

        void test_side_values(RangeConvert<DHEdgeSide, DHCellSide> range) {
            unsigned int k=0;
            for (DHCellSide cell_side : range) {
                auto p = *( this->edge_integral_->points(cell_side).begin() );

                double jxw = this->jxw_side_(p);
                arma::vec3 normal_vec = this->normal_vec_(p);
                EXPECT_TEST_NEAR( jxw, generic_inst_->ref_side_jxw_low_order_[k] );
                EXPECT_TEST_ARMA_NEAR( normal_vec, generic_inst_->ref_normal_vec_low_order_[k] );
                ++k;
            }
        }


    	/** Declaration of data members **/
        PatchFETestQuadOrders *generic_inst_;                                    ///< pointer to generic object
        FeQ<Scalar> jxw_diff_order_;
        FeQArray<Vector> vector_shape_diff_order_;
    };


    PatchFETestQuadOrders(unsigned int quad_order_1, unsigned int quad_order_2, std::shared_ptr<DOFHandlerMultiDim> dh_1, std::shared_ptr<DOFHandlerMultiDim> dh_2)
    : PatchFETestBase(dh_1),
      multidim_asm_(this, quad_order_1, quad_order_2)
    {
        element_cache_map_.init(eval_points_);

        ref_bulk_jxw_low_order_ = { 3.46410, 2.82843, 2.82843, 1.33333, 1.33333 };
        ref_bulk_jxw_high_order_ = {
            { 0.96225, 0.96225, 1.53960 },
            { 0.63182, 0.63182, 0.63182, 0.31099, 0.31099, 0.31099 },
            { 0.63182, 0.63182, 0.63182, 0.31099, 0.31099, 0.31099 },
            { 0.09799, 0.09799, 0.09799, 0.09799, 0.15025, 0.15025, 0.15025, 0.15025, 0.05673, 0.05673, 0.05673, 0.05673, 0.05673, 0.05673 },
            { 0.09799, 0.09799, 0.09799, 0.09799, 0.15025, 0.15025, 0.15025, 0.15025, 0.05673, 0.05673, 0.05673, 0.05673, 0.05673, 0.05673 }
        };
        ref_vector_shape_dof0_ = {
            { {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000} },
            { {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000},
              {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000} },
            { {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000},
              {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000} },
            { {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000},
              {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000},
              {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000} },
            { {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000},
              {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000},
              {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000}, {1.00000, 0.00000, 0.00000} }
        };
        ref_vector_shape_dof1_ = {
            { {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000} },
            { {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000},
              {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000} },
            { {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000},
              {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000} },
            { {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000},
              {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000},
              {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000} },
            { {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000},
              {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000},
              {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000}, {0.00000, 1.00000, 0.00000} }
        };
        ref_side_jxw_low_order_ = { 2.82843, 2.82843 };
        ref_normal_vec_low_order_ = { {-0.70711, 0.00000, -0.70711}, {0.70711, 0.00000, 0.70711} };
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
	    this->patch_fe_values_.init_finalize();
    }

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void update_patch() {
        this->resize_tables();
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
                for (auto p : edge_integrals_[dim-1]->points(edge_side, &element_cache_map_) ) {
                    patch_fe_values_.register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
        this->reinit_patch_fe();
    }

    void test_evaluation(bool print_tables=false) {
        for (auto elm_idx : used_element_idx_) {
            DHCellAccessor dh_cell = dh_->cell_accessor_from_element(elm_idx);
            auto &ppv_bulk = patch_fe_values_.ppv(bulk_domain, dh_cell.dim());
            ++ppv_bulk.n_mesh_items_;
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_[dh_cell.dim()-1]);
        	this->add_bulk_integral(dh_cell, this->bulk_integrals_diff_order_[dh_cell.dim()-1]);
        	this->add_edge_integral(dh_cell, this->edge_integrals_[dh_cell.dim()-1]);
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

        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
            RangeConvert<DHEdgeSide, DHCellSide> range = edge_integral_data_[i].edge_side_range;
            switch (range.begin()->dim()) {
            case 1:
                multidim_asm_[1_d]->test_side_values(range);
                break;
            case 2:
                multidim_asm_[2_d]->test_side_values(range);
                break;
            case 3:
                multidim_asm_[3_d]->test_side_values(range);
                break;
            }
        }
    }

    MixedPtr<AsmQuadOrders, 1> multidim_asm_;  ///< Assembly object

    std::vector<double> ref_bulk_jxw_low_order_;
    std::vector< std::vector<double> > ref_bulk_jxw_high_order_;
    std::vector< std::vector<arma::vec3> > ref_vector_shape_dof0_;
    std::vector< std::vector<arma::vec3> > ref_vector_shape_dof1_;
    std::vector<double> ref_side_jxw_low_order_;
    std::vector<arma::vec3> ref_normal_vec_low_order_;
};


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
