#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "quadrature/quadrature_lib.hh"
#include "fem/eval_subset.hh"
#include "fem/eval_points.hh"
#include "fem/element_cache_map.hh"
#include "fem/fe_values.hh"
#include "fem/patch_fe_values.hh"
#include "fem/op_factory.hh"
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


    PatchFETestBase(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : dh_(dh), patch_fe_values_(dh_->ds()->fe()),
      fe_(dh_->ds()->fe()), fe_values_(3), fe_values_side_(3),
      bulk_integral_data_(20, 10),
      edge_integral_data_(12, 6),
      coupling_integral_data_(12, 6),
      ref_quad_0_(new QGauss(0, 2*quad_order)),
      ref_quad_1_(new QGauss(1, 2*quad_order)),
      ref_quad_2_(new QGauss(2, 2*quad_order)),
      ref_quad_3_(new QGauss(3, 2*quad_order)),
      det_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).determinant() ),
      det_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).determinant() ),
      det_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).determinant() ),
      jxw_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).JxW() ),
      jxw_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).JxW() ),
      jxw_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).JxW() ),
      jxw_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).JxW() ),
      jxw_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).JxW() ),
      jxw_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).JxW() ),
      normal_vec_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).normal_vector() ),
      normal_vec_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).normal_vector() ),
      normal_vec_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).normal_vector() )
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize PatchFEValues on all dimensions
        this->create_integrals();
        element_cache_map_.init(eval_points_);

        UpdateFlags u = update_values | update_inverse_jacobians | update_JxW_values | update_quadrature_points | update_volume_elements | update_gradients;
        UpdateFlags u_side = update_values | update_inverse_jacobians | update_side_JxW_values | update_normal_vectors | update_quadrature_points | update_gradients;
        fe_values_[0].initialize(*ref_quad_1_, *fe_[Dim<1>{}], u);
        fe_values_[1].initialize(*ref_quad_2_, *fe_[Dim<2>{}], u);
        fe_values_[2].initialize(*ref_quad_3_, *fe_[Dim<3>{}], u);
        fe_values_side_[0].initialize(*ref_quad_0_, *fe_[Dim<1>{}], u_side);
        fe_values_side_[1].initialize(*ref_quad_1_, *fe_[Dim<2>{}], u_side);
        fe_values_side_[2].initialize(*ref_quad_2_, *fe_[Dim<3>{}], u_side);
    }

    ~PatchFETestBase() {}

    void create_integrals() {
        bulk_integrals_[0] = std::make_shared<BulkIntegral>(eval_points_, ref_quad_1_, 1);
        bulk_integrals_[1] = std::make_shared<BulkIntegral>(eval_points_, ref_quad_2_, 2);
        bulk_integrals_[2] = std::make_shared<BulkIntegral>(eval_points_, ref_quad_3_, 3);
        edge_integrals_[0] = std::make_shared<EdgeIntegral>(eval_points_, ref_quad_0_, 1);
        edge_integrals_[1] = std::make_shared<EdgeIntegral>(eval_points_, ref_quad_1_, 2);
        edge_integrals_[2] = std::make_shared<EdgeIntegral>(eval_points_, ref_quad_2_, 3);
        coupling_integrals_[0] = std::make_shared<CouplingIntegral>(eval_points_, ref_quad_1_, 1);
        coupling_integrals_[1] = std::make_shared<CouplingIntegral>(eval_points_, ref_quad_2_, 2);
    }

    void initialize() {
        this->patch_fe_values_.initialize<1>();
        this->patch_fe_values_.initialize<2>();
        this->patch_fe_values_.initialize<3>();
        this->patch_fe_values_.initialize<1>(false);
        this->patch_fe_values_.initialize<2>(false);
        this->patch_fe_values_.initialize<3>(false);
        this->patch_fe_values_.init_finalize();
    }

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

    void add_integrals(DHCellAccessor cell) {
        // Bulk integral
        uint subset_idx = bulk_integrals_[cell.dim()-1]->get_subset_idx();
        bulk_integral_data_.emplace_back(cell, subset_idx);
        uint dim = cell.dim();
        auto &ppv_bulk = patch_fe_values_.ppv(0, dim);
        ++ppv_bulk.n_elems_;
        ppv_bulk.n_points_ += eval_points_->subset_size(dim, subset_idx); // add rows for bulk points to table

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( eval_points_->subset_begin(dim, subset_idx) );
                  i<uint( eval_points_->subset_end(dim, subset_idx) ); ++i) {
            element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }

        // Edge integral
        auto &ppv_edge = patch_fe_values_.ppv(1, dim);
        for( DHCellSide cell_side : cell.side_range() ) {
            if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                auto range = cell_side.edge_sides();
                uint subset_idx = edge_integrals_[range.begin()->dim()-1]->get_subset_idx();
                edge_integral_data_.emplace_back(range, subset_idx);

                for( DHCellSide edge_side : range ) {
                    uint dim = edge_side.dim();
                    ++ppv_edge.n_elems_;
                    ppv_edge.n_points_ += eval_points_->subset_size(dim, subset_idx) / (dim+1); // add rows for side points to table
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : edge_integrals_[range.begin()->dim()-1]->points(edge_side, &element_cache_map_) ) {
                        element_cache_map_.add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    }
                }
            }
        }

//        // Coupling integral
//        bool add_low = true;
//    	for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
//            if (cell.dim() != neighb_side.dim()-1) continue;
//            coupling_integral_data_.emplace_back(cell, coupling_integrals_[cell.dim()-1]->get_subset_low_idx(), neighb_side,
//                    coupling_integrals_[cell.dim()-1]->get_subset_high_idx());
//            table_sizes_.elem_sizes_[1][cell.dim()]++;
//            if (add_low) table_sizes_.elem_sizes_[0][cell.dim()-1]++;
//
//            unsigned int reg_idx_low = cell.elm().region_idx().idx();
//            unsigned int reg_idx_high = neighb_side.element().region_idx().idx();
//            for (auto p : coupling_integrals_[cell.dim()-1]->points(neighb_side, &element_cache_map_) ) {
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
        patch_fe_values_.make_permanent_ppv_data();
    }

    void update_patch() {
        patch_fe_values_.resize_tables();
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
                uint side_pos = patch_fe_values_.register_side(edge_side);
                uint i_point = 0;
                for (auto p : this->edge_points(dim, edge_side) ) {
                    patch_fe_values_.register_side_point(edge_side, side_pos, p.value_cache_idx(), i_point++);
                }
            }
        }
        uint side_pos, element_patch_idx, elm_pos=0;
        uint last_element_idx = -1;
        for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
            uint dim = coupling_integral_data_[i].side.dim();
            side_pos = patch_fe_values_.register_side(coupling_integral_data_[i].side);
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
    PatchFEValues<3> patch_fe_values_;                                     ///< Common FEValues object over all dimensions

    MixedPtr<FiniteElement> fe_;
    std::vector<FEValues<3>> fe_values_;                                   ///< FeValues object of elements of dim 1,2,3
    std::vector<FEValues<3>> fe_values_side_;                              ///< FeValues object of sides of dim 0,1,2

    std::shared_ptr<EvalPoints> eval_points_;                              ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                                    ///< ElementCacheMap according to EvalPoints
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_;          ///< Bulk integrals of dim 1,2,3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_integrals_;          ///< Edge integrals of dim 1,2,3
    std::array<std::shared_ptr<CouplingIntegral>, 2> coupling_integrals_;  ///< Coupling integrals of dim 1-2,2-3
    RevertableList<BulkIntegralData> bulk_integral_data_;                  ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData> edge_integral_data_;                  ///< Holds data for computing edge integrals.
    RevertableList<CouplingIntegralData> coupling_integral_data_;          ///< Holds data for computing couplings integrals.
    Quadrature *ref_quad_0_;                                               ///< Reference quadrature used in FeValues
    Quadrature *ref_quad_1_;
    Quadrature *ref_quad_2_;
    Quadrature *ref_quad_3_;
    ElQ<Scalar> det_1d_;
    ElQ<Scalar> det_2d_;
    ElQ<Scalar> det_3d_;
    FeQ<Scalar> jxw_1d_;
    FeQ<Scalar> jxw_2d_;
    FeQ<Scalar> jxw_3d_;
    FeQ<Scalar> jxw_side_1d_;
    FeQ<Scalar> jxw_side_2d_;
    FeQ<Scalar> jxw_side_3d_;
    ElQ<Vector> normal_vec_1d_;
    ElQ<Vector> normal_vec_2d_;
    ElQ<Vector> normal_vec_3d_;
};


/**
 * Specialization defining FE scalar operations
 */
class PatchFETestScalar : public PatchFETestBase {
public:
	PatchFETestScalar(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : PatchFETestBase(quad_order, dh),
      scalar_shape_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).scalar_shape() ),
      scalar_shape_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).scalar_shape() ),
      scalar_shape_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).scalar_shape() ),
      scalar_shape_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).scalar_shape() ),
      scalar_shape_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).scalar_shape() ),
      scalar_shape_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).scalar_shape() ),
      grad_scalar_shape_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).grad_scalar_shape() ),
      grad_scalar_shape_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).grad_scalar_shape() ),
      grad_scalar_shape_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).grad_scalar_shape() ),
      grad_scalar_shape_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).grad_scalar_shape() ),
      grad_scalar_shape_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).grad_scalar_shape() ),
      grad_scalar_shape_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).grad_scalar_shape() ),
      conc_join_shape_2d_( FeQJoin<Scalar>( this->patch_fe_values_.template join_values<1>(ref_quad_1_).scalar_join_shape() ) ),
      conc_join_shape_3d_( FeQJoin<Scalar>( this->patch_fe_values_.template join_values<2>(ref_quad_2_).scalar_join_shape() ) )
    {}

    ~PatchFETestScalar() {}

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(bool print_tables=false) {
        for(auto cell_it = dh_->local_range().begin(); cell_it != dh_->local_range().end(); ++cell_it) {
            add_integrals(*cell_it);
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
            double det = 0.0, det_ref = 0.0;
            double scalar_shape_dof0 = 0.0, scalar_shape_dof0_ref = 0.0;
            double scalar_shape_dof1 = 0.0, scalar_shape_dof1_ref = 0.0;
            arma::vec3 grad_scalar_dof0("0 0 0");
            arma::vec3 grad_scalar_dof0_ref("0 0 0");
            arma::vec3 grad_scalar_dof1("0 0 0");
            arma::vec3 grad_scalar_dof1_ref("0 0 0");
            switch (dh_cell.dim()) {
            case 1:
                fe_values_[0].reinit(elm);
                jxw = jxw_1d_(p);
                det = det_1d_(p);
                scalar_shape_dof0 = scalar_shape_1d_.shape(0)(p);
                grad_scalar_dof0 = grad_scalar_shape_1d_.shape(0)(p);
                jxw_ref = fe_values_[0].JxW(0);
                det_ref = fe_values_[0].determinant(0);
                scalar_shape_dof0_ref = fe_values_[0].shape_value(0, 0);
                grad_scalar_dof0_ref = fe_values_[0].shape_grad(0, 0);
                break;
            case 2:
                fe_values_[1].reinit(elm);
                jxw = jxw_2d_(p);
                det = det_2d_(p);
                scalar_shape_dof0 = scalar_shape_2d_.shape(0)(p);
                scalar_shape_dof1 = scalar_shape_2d_.shape(1)(p);
                grad_scalar_dof0 = grad_scalar_shape_2d_.shape(0)(p);
                grad_scalar_dof1 = grad_scalar_shape_2d_.shape(1)(p);
                jxw_ref = fe_values_[1].JxW(0);
                det_ref = fe_values_[1].determinant(0);
                scalar_shape_dof0_ref = fe_values_[1].shape_value(0, 0);
                scalar_shape_dof1_ref = fe_values_[1].shape_value(1, 0);
                grad_scalar_dof0_ref = fe_values_[1].shape_grad(0, 0);
                grad_scalar_dof1_ref = fe_values_[1].shape_grad(1, 0);
                break;
            case 3:
                fe_values_[2].reinit(elm);
                jxw = jxw_3d_(p);
                det = det_3d_(p);
                scalar_shape_dof0 = scalar_shape_3d_.shape(0)(p);
                scalar_shape_dof1 = scalar_shape_3d_.shape(1)(p);
                grad_scalar_dof0 = grad_scalar_shape_3d_.shape(0)(p);
                grad_scalar_dof1 = grad_scalar_shape_3d_.shape(1)(p);
                jxw_ref = fe_values_[2].JxW(0);
                det_ref = fe_values_[2].determinant(0);
                scalar_shape_dof0_ref = fe_values_[2].shape_value(0, 0);
                scalar_shape_dof1_ref = fe_values_[2].shape_value(1, 0);
                grad_scalar_dof0_ref = fe_values_[2].shape_grad(0, 0);
                grad_scalar_dof1_ref = fe_values_[2].shape_grad(1, 0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, jxw_ref );
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
                jxw = jxw_side_1d_(p);
                normal_vec = normal_vec_1d_(p);
                scalar_shape = scalar_shape_side_1d_.shape(0)(p);
                grad_scalar = grad_scalar_shape_side_1d_.shape(0)(p);
                fe_values_side_[0].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[0].JxW(0);
                normal_vec_ref = fe_values_side_[0].normal_vector(0);
                scalar_shape_ref = fe_values_side_[0].shape_value(0, 0);
                grad_scalar_ref = fe_values_side_[0].shape_grad(0, 0);
                break;
            case 2:
                jxw = jxw_side_2d_(p);
                normal_vec = normal_vec_2d_(p);
                scalar_shape = scalar_shape_side_2d_.shape(0)(p);
                grad_scalar = grad_scalar_shape_side_2d_.shape(0)(p);
                fe_values_side_[1].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[1].JxW(0);
                normal_vec_ref = fe_values_side_[1].normal_vector(0);
                scalar_shape_ref = fe_values_side_[1].shape_value(0, 0);
                grad_scalar_ref = fe_values_side_[1].shape_grad(0, 0);
                break;
            case 3:
                jxw = jxw_side_3d_(p);
                normal_vec = normal_vec_3d_(p);
                scalar_shape = scalar_shape_side_3d_.shape(0)(p);
                grad_scalar = grad_scalar_shape_side_3d_.shape(0)(p);
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
                for (uint i_dof=0; i_dof<conc_join_shape_2d_.n_dofs_both(); ++i_dof) {
                    if (conc_join_shape_2d_.is_high_dim(i_dof)) {
                        auto result = conc_join_shape_2d_.shape(i_dof)(p_high);
                        auto ref = fe_values_side_[1].shape_value(i_dof_high, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_low = conc_join_shape_2d_.shape(i_dof)(p_low);
                        EXPECT_DOUBLE_EQ( result_low, 0.0 );
                    	i_dof_high++;
                    }
                    else {
                        auto result = conc_join_shape_2d_.shape(i_dof)(p_low);
                        auto ref = fe_values_[0].shape_value(i_dof_low, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_high = conc_join_shape_2d_.shape(i_dof)(p_high);
                        EXPECT_DOUBLE_EQ( result_high, 0.0 );
                        i_dof_low++;
                    }
                }
                break;
            case 3:
                fe_values_[1].reinit(cell_lower_dim.elm());
                fe_values_side_[2].reinit(neighb_side.side());
                for (uint i_dof=0; i_dof<conc_join_shape_3d_.n_dofs_both(); ++i_dof) {
                    if (conc_join_shape_3d_.is_high_dim(i_dof)) {
                        auto result = conc_join_shape_3d_.shape(i_dof)(p_high);
                        auto ref = fe_values_side_[2].shape_value(i_dof_high, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_low = conc_join_shape_3d_.shape(i_dof)(p_low);
                        EXPECT_DOUBLE_EQ( result_low, 0.0 );
                        i_dof_high++;
                    }
                    else {
                   	    auto result = conc_join_shape_3d_.shape(i_dof)(p_low);
                        auto ref = fe_values_[1].shape_value(i_dof_low, 0);
                        EXPECT_DOUBLE_EQ( result, ref );
                        auto result_high = conc_join_shape_3d_.shape(i_dof)(p_high);
                        EXPECT_DOUBLE_EQ( result_high, 0.0 );
                        i_dof_low++;
                    }
                }
                break;
            }
        }

    }

    FeQArray<Scalar> scalar_shape_1d_;
    FeQArray<Scalar> scalar_shape_2d_;
    FeQArray<Scalar> scalar_shape_3d_;
    FeQArray<Scalar> scalar_shape_side_1d_;
    FeQArray<Scalar> scalar_shape_side_2d_;
    FeQArray<Scalar> scalar_shape_side_3d_;
    FeQArray<Vector> grad_scalar_shape_1d_;
    FeQArray<Vector> grad_scalar_shape_2d_;
    FeQArray<Vector> grad_scalar_shape_3d_;
    FeQArray<Vector> grad_scalar_shape_side_1d_;
    FeQArray<Vector> grad_scalar_shape_side_2d_;
    FeQArray<Vector> grad_scalar_shape_side_3d_;
    FeQJoin<Scalar> conc_join_shape_2d_;
    FeQJoin<Scalar> conc_join_shape_3d_;
};


/**
 * Specialization defining FE vector operations
 */
class PatchFETestVector : public PatchFETestBase {
public:
	PatchFETestVector(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : PatchFETestBase(quad_order, dh),
      vector_shape_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).vector_shape() ),
      vector_shape_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).vector_shape() ),
      vector_shape_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).vector_shape() ),
      vector_shape_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).vector_shape() ),
      vector_shape_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).vector_shape() ),
      vector_shape_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).vector_shape() ),
      grad_vector_shape_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).grad_vector_shape() ),
      grad_vector_shape_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).grad_vector_shape() ),
      grad_vector_shape_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).grad_vector_shape() ),
      grad_vector_shape_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).grad_vector_shape() ),
      grad_vector_shape_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).grad_vector_shape() ),
      grad_vector_shape_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).grad_vector_shape() ),
      sym_grad_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).vector_sym_grad() ),
      sym_grad_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).vector_sym_grad() ),
      sym_grad_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).vector_sym_grad() ),
      sym_grad_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).vector_sym_grad() ),
      sym_grad_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).vector_sym_grad() ),
      sym_grad_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).vector_sym_grad() ),
      divergence_1d_( this->patch_fe_values_.bulk_values<1>(ref_quad_1_).vector_divergence() ),
      divergence_2d_( this->patch_fe_values_.bulk_values<2>(ref_quad_2_).vector_divergence() ),
      divergence_3d_( this->patch_fe_values_.bulk_values<3>(ref_quad_3_).vector_divergence() ),
      divergence_side_1d_( this->patch_fe_values_.side_values<1>(ref_quad_0_).vector_divergence() ),
      divergence_side_2d_( this->patch_fe_values_.side_values<2>(ref_quad_1_).vector_divergence() ),
      divergence_side_3d_( this->patch_fe_values_.side_values<3>(ref_quad_2_).vector_divergence() ),
      vector_join_2d_( this->patch_fe_values_.join_values<1>(ref_quad_1_).vector_join_shape() ),
      vector_join_3d_( this->patch_fe_values_.join_values<2>(ref_quad_2_).vector_join_shape() ),
      vector_join_grad_2d_( this->patch_fe_values_.join_values<1>(ref_quad_1_).gradient_vector_join_shape() ),
      vector_join_grad_3d_( this->patch_fe_values_.join_values<2>(ref_quad_2_).gradient_vector_join_shape() )
    {
	    vec_view_1d_ = &fe_values_[0].vector_view(0);
	    vec_view_2d_ = &fe_values_[1].vector_view(0);
	    vec_view_3d_ = &fe_values_[2].vector_view(0);
	    vec_view_side_1d_ = &fe_values_side_[0].vector_view(0);
	    vec_view_side_2d_ = &fe_values_side_[1].vector_view(0);
	    vec_view_side_3d_ = &fe_values_side_[2].vector_view(0);
    }

    ~PatchFETestVector() {}

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch");
        patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch");
    }

    void test_evaluation(bool print_tables=false) {
        for(auto cell_it = dh_->local_range().begin(); cell_it != dh_->local_range().end(); ++cell_it) {
            add_integrals(*cell_it);
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
                jxw = jxw_1d_(p);
                vector_shape_dof0 = vector_shape_1d_.shape(0)(p);
                grad_vector_dof0 = grad_vector_shape_1d_.shape(0)(p);
                sym_grad_dof0 = sym_grad_1d_.shape(0)(p);
                div_dof0 = divergence_1d_.shape(0)(p);
                jxw_ref = fe_values_[0].JxW(0);
                vector_shape_dof0_ref = vec_view_1d_->value(0, 0);
                grad_vector_dof0_ref = vec_view_1d_->grad(0, 0);
                sym_grad_dof0_ref = vec_view_1d_->sym_grad(0, 0);
                div_dof0_ref = vec_view_1d_->divergence(0, 0);
                break;
            case 2:
                fe_values_[1].reinit(elm);
                jxw = jxw_2d_(p);
                vector_shape_dof0 = vector_shape_2d_.shape(0)(p);
                vector_shape_dof1 = vector_shape_2d_.shape(1)(p);
                grad_vector_dof0 = grad_vector_shape_2d_.shape(0)(p);
                grad_vector_dof1 = grad_vector_shape_2d_.shape(1)(p);
                sym_grad_dof0 = sym_grad_2d_.shape(0)(p);
                sym_grad_dof1 = sym_grad_2d_.shape(1)(p);
                div_dof0 = divergence_2d_.shape(0)(p);
                div_dof1 = divergence_2d_.shape(1)(p);
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
                jxw = jxw_3d_(p);
                vector_shape_dof0 = vector_shape_3d_.shape(0)(p);
                vector_shape_dof1 = vector_shape_3d_.shape(1)(p);
                grad_vector_dof0 = grad_vector_shape_3d_.shape(0)(p);
                grad_vector_dof1 = grad_vector_shape_3d_.shape(1)(p);
                sym_grad_dof0 = sym_grad_3d_.shape(0)(p);
                sym_grad_dof1 = sym_grad_3d_.shape(1)(p);
                div_dof0 = divergence_3d_.shape(0)(p);
                div_dof1 = divergence_3d_.shape(1)(p);
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
                jxw = jxw_side_1d_(p);
                vector_shape = vector_shape_side_1d_.shape(0)(p);
                grad_vector = grad_vector_shape_side_1d_.shape(0)(p);
                sym_grad = sym_grad_side_1d_.shape(0)(p);
                div = divergence_side_1d_.shape(0)(p);
                fe_values_side_[0].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[0].JxW(0);
                vector_shape_ref = vec_view_side_1d_->value(0, 0);
                grad_vector_ref = vec_view_side_1d_->grad(0, 0);
                sym_grad_ref = vec_view_side_1d_->sym_grad(0, 0);
                div_ref = vec_view_side_1d_->divergence(0, 0);
                break;
            case 2:
                jxw = jxw_side_2d_(p);
                vector_shape = vector_shape_side_2d_.shape(0)(p);
                grad_vector = grad_vector_shape_side_2d_.shape(0)(p);
                sym_grad = sym_grad_side_2d_.shape(0)(p);
                div = divergence_side_2d_.shape(0)(p);
                fe_values_side_[1].reinit(zero_edge_side.side());
                jxw_ref = fe_values_side_[1].JxW(0);
                vector_shape_ref = vec_view_side_2d_->value(0, 0);
                grad_vector_ref = vec_view_side_2d_->grad(0, 0);
                sym_grad_ref = vec_view_side_2d_->sym_grad(0, 0);
                div_ref = vec_view_side_2d_->divergence(0, 0);
                break;
            case 3:
                jxw = jxw_side_3d_(p);
                vector_shape = vector_shape_side_3d_.shape(0)(p);
                grad_vector = grad_vector_shape_side_3d_.shape(0)(p);
                sym_grad = sym_grad_side_3d_.shape(0)(p);
                div = divergence_side_3d_.shape(0)(p);
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
                for (uint i_dof=0; i_dof<vector_join_2d_.n_dofs_both(); ++i_dof) {
                    if (vector_join_2d_.is_high_dim(i_dof)) {
                        auto result = vector_join_2d_.shape(i_dof)(p_high);
                        auto ref = vec_view_side_2d_->value(i_dof_high, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_low = vector_join_2d_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( result_low, arma_zero_vec );

                        auto grad_result = vector_join_grad_2d_.shape(i_dof)(p_high);
                        auto grad_ref = vec_view_side_2d_->grad(i_dof_high, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_low = vector_join_grad_2d_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( grad_result_low, arma_zero_mat );
                        i_dof_high++;
                    }
                    else {
                        auto result = vector_join_2d_.shape(i_dof)(p_low);
                        auto ref = vec_view_1d_->value(i_dof_low, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_high = vector_join_2d_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( result_high, arma_zero_vec );

                        auto grad_result = vector_join_grad_2d_.shape(i_dof)(p_low);
                        auto grad_ref = vec_view_1d_->grad(i_dof_low, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_high = vector_join_grad_2d_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( grad_result_high, arma_zero_mat );
                        i_dof_low++;
                    }
                }
                break;
            case 3:
                fe_values_[1].reinit(cell_lower_dim.elm());
                fe_values_side_[2].reinit(neighb_side.side());
                for (uint i_dof=0; i_dof<vector_join_3d_.n_dofs_both(); ++i_dof) {
                    if (vector_join_3d_.is_high_dim(i_dof)) {
                        auto result = vector_join_3d_.shape(i_dof)(p_high);
                        auto ref = vec_view_side_3d_->value(i_dof_high, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_low = vector_join_3d_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( result_low, arma_zero_vec );

                        auto grad_result = vector_join_grad_3d_.shape(i_dof)(p_high);
                        auto grad_ref = vec_view_side_3d_->grad(i_dof_high, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_low = vector_join_grad_3d_.shape(i_dof)(p_low);
                        EXPECT_ARMA_EQ( grad_result_low, arma_zero_mat );
                        i_dof_high++;
                    }
                    else {
                        auto result = vector_join_3d_.shape(i_dof)(p_low);
                        auto ref = vec_view_2d_->value(i_dof_low, 1);
                        EXPECT_ARMA_EQ( result, ref );
                        auto result_high = vector_join_3d_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( result_high, arma_zero_vec );

                        auto grad_result = vector_join_grad_3d_.shape(i_dof)(p_low);
                        auto grad_ref = vec_view_2d_->grad(i_dof_low, 1);
                        EXPECT_ARMA_EQ( grad_result, grad_ref );
                        auto grad_result_high = vector_join_grad_3d_.shape(i_dof)(p_high);
                        EXPECT_ARMA_EQ( grad_result_high, arma_zero_mat );
                        i_dof_low++;
                    }
                }
                break;
            }
        }

    }

    ///< Vector view in cell calculation.
    const FEValuesViews::Vector<3> * vec_view_1d_;
    const FEValuesViews::Vector<3> * vec_view_2d_;
    const FEValuesViews::Vector<3> * vec_view_3d_;
    ///< Vector view in SIDE calculation.
    const FEValuesViews::Vector<3> * vec_view_side_1d_;
    const FEValuesViews::Vector<3> * vec_view_side_2d_;
    const FEValuesViews::Vector<3> * vec_view_side_3d_;

    FeQArray<Vector> vector_shape_1d_;
    FeQArray<Vector> vector_shape_2d_;
    FeQArray<Vector> vector_shape_3d_;
    FeQArray<Vector> vector_shape_side_1d_;
    FeQArray<Vector> vector_shape_side_2d_;
    FeQArray<Vector> vector_shape_side_3d_;
    FeQArray<Tensor> grad_vector_shape_1d_;
    FeQArray<Tensor> grad_vector_shape_2d_;
    FeQArray<Tensor> grad_vector_shape_3d_;
    FeQArray<Tensor> grad_vector_shape_side_1d_;
    FeQArray<Tensor> grad_vector_shape_side_2d_;
    FeQArray<Tensor> grad_vector_shape_side_3d_;
    FeQArray<Tensor> sym_grad_1d_;
    FeQArray<Tensor> sym_grad_2d_;
    FeQArray<Tensor> sym_grad_3d_;
    FeQArray<Tensor> sym_grad_side_1d_;
    FeQArray<Tensor> sym_grad_side_2d_;
    FeQArray<Tensor> sym_grad_side_3d_;
    FeQArray<Scalar> divergence_1d_;
    FeQArray<Scalar> divergence_2d_;
    FeQArray<Scalar> divergence_3d_;
    FeQArray<Scalar> divergence_side_1d_;
    FeQArray<Scalar> divergence_side_2d_;
    FeQArray<Scalar> divergence_side_3d_;
    FeQJoin<Vector> vector_join_2d_;
    FeQJoin<Vector> vector_join_3d_;
    FeQJoin<Tensor> vector_join_grad_2d_;
    FeQJoin<Tensor> vector_join_grad_3d_;
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
void compare_evaluation_func_scalar(Mesh* mesh, unsigned int quad_order, bool print_fa_data = false) {
	MixedPtr<FE_P_disc> fe(quad_order);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    PatchFETestScalar patch_fe(quad_order, dh);
    patch_fe.initialize();
    patch_fe.test_evaluation(print_fa_data);
    patch_fe.reset();
    patch_fe.test_evaluation();
}

/// Complete test with FE scalar operations
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



TEST(PatchFeTest, complete_evaluation) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    std::string input_str = "{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }";
    Mesh* mesh = mesh_full_constructor(input_str);

    // two tests with different quad_order and Scalar / Vector FE operations
    compare_evaluation_func_scalar(mesh, 1, true);
    compare_evaluation_func_scalar(mesh, 2);
    compare_evaluation_func_vector(mesh, 1, true);
}

//TEST(PatchFeTest, speed_comparation) {
//    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
//    Profiler::instance();
//    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);
//
//    std::string input_str = "{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }";
//    Mesh* mesh = mesh_full_constructor(input_str);
//
//    MixedPtr<FE_P_disc> fe(1);
//    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
//    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
//    dh->distribute_dofs(ds);
//    unsigned int quad_order = 1;
//
//    PatchFETestCompare patch_fe(quad_order, dh);
//    patch_fe.initialize();
//    for(auto cell_it = dh->local_range().begin(); cell_it != dh->local_range().end(); ++cell_it) {
//    	patch_fe.add_integrals(*cell_it);
//    }
//    patch_fe.bulk_integral_data_.make_permanent();
//    patch_fe.edge_integral_data_.make_permanent();
//    patch_fe.coupling_integral_data_.make_permanent();
//    patch_fe.element_cache_map_.make_paermanent_eval_points();
//    patch_fe.element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
//    patch_fe.update_patch();
//
//    patch_fe.reinit_fe_values();
//
//    patch_fe.profiler_output("patch_fe");
//
//}


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

//TEST(PatchFeTest, multi_quad_orders_simple) {
//    MixedPtr<FE_P> fe_p(1);
//    MixedPtr<FiniteElement> fe = mixed_fe_system(fe_p, FEVector, 3);
//    PatchFEValues<3> pfev(fe);
//    Quadrature *quad = new QGauss(2, 0);
//    Quadrature *quad_high_order = new QGauss(2, 2);
//
//    auto fe_component = fe[2_d];
//    pfev.template get< Op::JxW<2, Op::BulkDomain, 3>, 2 >(quad);
//    pfev.template get< Op::DispatchGradVectorShape<2, Op::BulkDomain, 3>, 2 >(quad, fe_component);
//    pfev.template get< Op::DispatchGradVectorShape<2, Op::BulkDomain, 3>, 2 >(quad_high_order, fe_component);
//}


/*************************************************
 * New test of operation dependency
 */

class PatchFEMultiQuadTest {
public:

    PatchFEMultiQuadTest(std::shared_ptr<DOFHandlerMultiDim> dh)
    : dh_(dh), patch_fe_values_(dh_->ds()->fe()),
	  fe_(dh_->ds()->fe()),
	  //fe_values_(3), fe_values_side_(3),
	  eval_points_( std::make_shared<EvalPoints>() ),
	  quad_0d_low_order_(new QGauss(0, 0)),
	  quad_1d_low_order_(new QGauss(1, 0)),
	  quad_1d_high_order_(new QGauss(1, 2)),
	  quad_2d_low_order_(new QGauss(2, 0)),
	  bulk_int_1d_low_order_( this->create_bulk_integral(quad_1d_low_order_) ),
	  bulk_int_1d_high_order_( this->create_bulk_integral(quad_1d_high_order_) ),
	  bulk_int_2d_low_order_( this->create_bulk_integral(quad_2d_low_order_) ),
	  edge_int_1d_low_order_( this->create_edge_integral(quad_0d_low_order_) ),
	  edge_int_2d_low_order_( this->create_edge_integral(quad_1d_low_order_) ),
	  edge_int_2d_high_order_( this->create_edge_integral(quad_1d_high_order_) ),
	  coupling_int_1d_2d_low_order_( this->create_coupling_integral(quad_1d_low_order_) ),
      bulk_integral_data_(20, 10),
      edge_integral_data_(20, 10),
	  coupling_integral_data_(20, 10),
      det_1d_low_order_( this->patch_fe_values_.bulk_values<1>(quad_1d_low_order_).determinant() ),
      det_1d_high_order_( this->patch_fe_values_.bulk_values<1>(quad_1d_high_order_).determinant() ),
      det_2d_low_order_( this->patch_fe_values_.bulk_values<1>(quad_2d_low_order_).determinant() )
    {
        element_cache_map_.init(eval_points_);
        initialize();

//        UpdateFlags u = update_values | update_inverse_jacobians | update_JxW_values | update_quadrature_points | update_volume_elements | update_gradients;
//        UpdateFlags u_side = update_values | update_inverse_jacobians | update_side_JxW_values | update_normal_vectors | update_quadrature_points | update_gradients;
//        fe_values_[0].initialize(*patch_fe_values_.get_quadrature(1, true), *fe_[Dim<1>{}], u);
//        fe_values_[1].initialize(*patch_fe_values_.get_quadrature(2, true), *fe_[Dim<2>{}], u);
//        fe_values_[2].initialize(*patch_fe_values_.get_quadrature(3, true), *fe_[Dim<3>{}], u);
//        fe_values_side_[0].initialize(*patch_fe_values_.get_quadrature(1, false), *fe_[Dim<1>{}], u_side);
//        fe_values_side_[1].initialize(*patch_fe_values_.get_quadrature(2, false), *fe_[Dim<2>{}], u_side);
//        fe_values_side_[2].initialize(*patch_fe_values_.get_quadrature(3, false), *fe_[Dim<3>{}], u_side);
    }

    ~PatchFEMultiQuadTest() {}

    void initialize() {
        this->patch_fe_values_.initialize<1>(true);
        this->patch_fe_values_.initialize<2>(true);
//        this->patch_fe_values_.initialize<3>(true);
        this->patch_fe_values_.initialize<1>(false);
        this->patch_fe_values_.initialize<2>(false);
//        this->patch_fe_values_.initialize<3>(false);
        this->patch_fe_values_.init_finalize();
    }

    /// Return point range of appropriate integral
    template <class QIntegral>
    Range< typename QIntegral::PointType > points(std::shared_ptr<QIntegral> integral, typename QIntegral::MeshItem mesh_item) const {
    	return integral->points(mesh_item, &element_cache_map_);
    }

    /// Create bulk integral - used in constructor initializer list
    std::shared_ptr<BulkIntegral> create_bulk_integral(Quadrature *quad) {
        auto result = integrals_.bulk_.insert( std::make_shared<BulkIntegral>(eval_points_, quad, quad->dim()) );
	    return *result.first;
    }

    /// Create edge integral - used in constructor initializer list
    std::shared_ptr<EdgeIntegral> create_edge_integral(Quadrature *quad) {
        auto result = integrals_.edge_.insert( std::make_shared<EdgeIntegral>(eval_points_, quad, quad->dim()+1) );
	    return *result.first;
    }

    /// Create coupling integral - used in constructor initializer list
    std::shared_ptr<CouplingIntegral> create_coupling_integral(Quadrature *quad) {
        auto result = integrals_.coupling_.insert( std::make_shared<CouplingIntegral>(eval_points_, quad, quad->dim()) );
	    return *result.first;
    }

	/// Perform profiler output.
    void profiler_output(std::string file_name) {
		FilePath fp(file_name + "_profiler.json", FilePath::output_file);
		Profiler::instance()->output(MPI_COMM_WORLD, fp.filename());
	}

    /// reset patch data
    void reset() {
        element_cache_map_.clear_element_eval_points_map();
        patch_fe_values_.reset();
    }

    void test_evaluation() {
        for(auto cell_it = dh_->local_range().begin(); cell_it != dh_->local_range().end(); ++cell_it) {
            if (cell_it->dim() == 1) {
                add_bulk_integral(*cell_it, bulk_int_1d_low_order_);
                add_bulk_integral(*cell_it, bulk_int_1d_high_order_);
                add_edge_integral(*cell_it, edge_int_1d_low_order_);
                add_coupling_integral(*cell_it, coupling_int_1d_2d_low_order_);
            } else if (cell_it->dim() == 2) {
                add_bulk_integral(*cell_it, bulk_int_2d_low_order_);
                add_edge_integral(*cell_it, edge_int_2d_low_order_);
                add_edge_integral(*cell_it, edge_int_2d_high_order_);
            }
        }
        bulk_int_1d_low_order_->ppv().make_permanent_mesh_items();
        bulk_int_1d_high_order_->ppv().make_permanent_mesh_items();
        bulk_int_2d_low_order_->ppv().make_permanent_mesh_items();
        edge_int_1d_low_order_->ppv_side().make_permanent_mesh_items();
        edge_int_2d_low_order_->ppv_side().make_permanent_mesh_items();
        edge_int_2d_high_order_->ppv_side().make_permanent_mesh_items();
        coupling_int_1d_2d_low_order_->ppv().make_permanent_mesh_items();
        coupling_int_1d_2d_low_order_->ppv_side().make_permanent_mesh_items();
        element_cache_map_.make_paermanent_eval_points();
        element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
        std::cout << "Bulk 1D low   " << bulk_int_1d_low_order_->ppv().n_elems() << " - " << bulk_int_1d_low_order_->ppv().n_points() << std::endl;
        std::cout << "Bulk 1D high  " << bulk_int_1d_high_order_->ppv().n_elems() << " - " << bulk_int_1d_high_order_->ppv().n_points() << std::endl;
        std::cout << "Bulk 2D low   " << bulk_int_2d_low_order_->ppv().n_elems() << " - " << bulk_int_2d_low_order_->ppv().n_points() << std::endl;
        std::cout << "Edge 1D low   " << edge_int_1d_low_order_->ppv_side().n_elems() << " - " << edge_int_1d_low_order_->ppv_side().n_points() << std::endl;
        std::cout << "Edge 2D low   " << edge_int_2d_low_order_->ppv_side().n_elems() << " - " << edge_int_2d_low_order_->ppv_side().n_points() << std::endl;
        std::cout << "Edge 2D high  " << edge_int_2d_high_order_->ppv_side().n_elems() << " - " << edge_int_2d_high_order_->ppv_side().n_points() << std::endl;
        std::cout << "Coupling low  " << coupling_int_1d_2d_low_order_->ppv().n_elems() << " - " << coupling_int_1d_2d_low_order_->ppv().n_points() << std::endl;
        std::cout << "Coupling high " << coupling_int_1d_2d_low_order_->ppv_side().n_elems() << " - " << coupling_int_1d_2d_low_order_->ppv_side().n_points() << std::endl;

        update_patch();
    }

    void add_bulk_integral(DHCellAccessor cell, std::shared_ptr<BulkIntegral> bulk_int) {
        uint subset_idx = bulk_int->get_subset_idx();
        bulk_integral_data_.emplace_back(cell, subset_idx);
        uint dim = cell.dim();
        auto &ppv_bulk = bulk_int->ppv();
        ++ppv_bulk.n_elems_;
        ppv_bulk.n_points_ += eval_points_->subset_size(dim, subset_idx); // add rows for bulk points to table

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( eval_points_->subset_begin(dim, subset_idx) );
                  i<uint( eval_points_->subset_end(dim, subset_idx) ); ++i) {
            element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }
    }

    void add_edge_integral(DHCellAccessor cell, std::shared_ptr<EdgeIntegral> edge_int) {
        auto &ppv_edge = edge_int->ppv_side();
        for( DHCellSide cell_side : cell.side_range() ) {
            if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                auto range = cell_side.edge_sides();
                uint subset_idx = edge_int->get_subset_idx();
                edge_integral_data_.emplace_back(range, subset_idx);

                for( DHCellSide edge_side : range ) {
                    uint dim = edge_side.dim();
                    ++ppv_edge.n_elems_;
                    ppv_edge.n_points_ += eval_points_->subset_size(dim, subset_idx) / (dim+1); // add rows for side points to table
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : edge_int->points(edge_side, &element_cache_map_) ) {
                        element_cache_map_.add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    }
                }
            }
        }
    }

    void add_coupling_integral(DHCellAccessor cell, std::shared_ptr<CouplingIntegral> coupling_int) {
        auto &ppv_low = coupling_int->ppv();
        auto &ppv_high = coupling_int->ppv_side();

//        // Adds data of bulk points only if bulk point were not added during processing of bulk integral
//        bool add_bulk_points = !( (integrals_.bulk_.size() > 0) & cell.is_own() );
        // add points of low dim element only one time and only if they have not been added in BulkIntegral
        for( DHCellSide ngh_side : cell.neighb_sides() ) {
            unsigned int reg_idx_low = cell.elm().region_idx().idx();
            ++ppv_low.n_elems_;
            for (auto p : coupling_int->points(ngh_side, &element_cache_map_) ) {
//                if (add_bulk_points) {
                    auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                    element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
//                }
                ++ppv_low.n_points_;
            }
            break;
        }

        // Adds data of side points of all neighbour objects
        for( DHCellSide ngh_side : cell.neighb_sides() ) { // cell -> elm lower dim, ngh_side -> elm higher dim
            coupling_integral_data_.emplace_back(cell, coupling_int->get_subset_low_idx(), ngh_side,
                    coupling_int->get_subset_high_idx());
            ++ppv_high.n_elems_;

            unsigned int reg_idx_high = ngh_side.element().region_idx().idx();
            for (auto p : coupling_int->points(ngh_side, &element_cache_map_) ) {
                element_cache_map_.add_eval_point(reg_idx_high, ngh_side.elem_idx(), p.eval_point_idx(), ngh_side.cell().local_idx());
                ++ppv_high.n_points_;
            }
        }
    }

    void update_patch() {
        // register points of bulk integral
        for (auto integral_it : integrals_.bulk_) {
            auto &ppv = integral_it->ppv();
		    ppv.resize_tables( patch_fe_values_.patch_arena() );

            for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
                if ( bulk_integral_data_[i].subset_index != (unsigned int)(integral_it->get_subset_idx()) ) continue;
                uint element_patch_idx = element_cache_map_.position_in_cache(bulk_integral_data_[i].cell.elm_idx());
                uint elm_pos = ppv.register_element(bulk_integral_data_[i].cell, element_patch_idx);
                uint i_point = 0;
                for (auto p : integral_it->points(element_patch_idx, &element_cache_map_) ) {
                    ppv.register_bulk_point(elm_pos, p.value_cache_idx(), bulk_integral_data_[i].cell.elm_idx(), i_point++);
                }
            }
        }

    	// register points of edge integral
        for (auto integral_it : integrals_.edge_) {
            auto &ppv = integral_it->ppv_side();
		    ppv.resize_tables( patch_fe_values_.patch_arena() );

            for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
                if ( edge_integral_data_[i].subset_index != (unsigned int)(integral_it->get_subset_idx()) ) continue;
            	auto range = edge_integral_data_[i].edge_side_range;
                for( DHCellSide edge_side : range )
                {
                	uint side_pos = ppv.register_side(edge_side);
                    uint i_point = 0;
                    for (auto p : integral_it->points(edge_side, &element_cache_map_) ) {
                        ppv.register_side_point(side_pos, p.value_cache_idx(), edge_side.elem_idx(), edge_side.side_idx(), i_point++);
                    }
                }
            }
        }

        // add coupling points
        for (auto integral_it : integrals_.coupling_) {
            uint side_pos, element_patch_idx, elm_pos=0;
            uint last_element_idx = -1;
            auto &ppv_low = integral_it->ppv();
            auto &ppv_high = integral_it->ppv_side();

            for (unsigned int i=0; i<coupling_integral_data_.permanent_size(); ++i) {
                if ( coupling_integral_data_[i].bulk_subset_index != (unsigned int)(integral_it->get_subset_low_idx()) ) continue;
                side_pos = ppv_high.register_side(coupling_integral_data_[i].side);
                if (coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                    element_patch_idx = element_cache_map_.position_in_cache(coupling_integral_data_[i].cell.elm_idx());
                    elm_pos = ppv_low.register_element(coupling_integral_data_[i].cell, element_patch_idx);
                }

                uint i_bulk_point = 0, i_side_point = 0;
                for (auto p_high : integral_it->points(coupling_integral_data_[i].side, &element_cache_map_) )
                {
                    ppv_high.register_side_point(side_pos, p_high.value_cache_idx(), coupling_integral_data_[i].side.elem_idx(),
                            coupling_integral_data_[i].side.side_idx(), i_side_point++);
                    if (coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                        auto p_low = p_high.lower_dim(coupling_integral_data_[i].cell);
                        ppv_low.register_bulk_point(elm_pos, p_low.value_cache_idx(), coupling_integral_data_[i].cell.elm_idx(), i_bulk_point++);
                    }
                }
                last_element_idx = coupling_integral_data_[i].cell.elm_idx();
            }
        }
    }


    /** Data members **/

    std::shared_ptr<DOFHandlerMultiDim> dh_;
    PatchFEValues<3> patch_fe_values_;                                     ///< Common FEValues object over all dimensions

    MixedPtr<FiniteElement> fe_;
//    std::vector<FEValues<3>> fe_values_;                                   ///< FeValues object of elements of dim 1,2,3
//    std::vector<FEValues<3>> fe_values_side_;                              ///< FeValues object of sides of dim 0,1,2

    DimIntegrals integrals_;
    std::shared_ptr<EvalPoints> eval_points_;                                      ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                                            ///< ElementCacheMap according to EvalPoints
    Quadrature *quad_0d_low_order_;                                                ///< Quadrature - 0D, order = 0
    Quadrature *quad_1d_low_order_;                                                ///< Quadrature - 1D, order = 0
    Quadrature *quad_1d_high_order_;                                               ///< Quadrature - 1D, order = 2
    Quadrature *quad_2d_low_order_;                                                ///< Quadrature - 2D, order = 0
    std::shared_ptr<BulkIntegral> bulk_int_1d_low_order_;                          ///< Bulk integral - 1D, order = 0
    std::shared_ptr<BulkIntegral> bulk_int_1d_high_order_;                         ///< Bulk integral - 1D, order = 2
    std::shared_ptr<BulkIntegral> bulk_int_2d_low_order_;                          ///< Bulk integral - 2D, order = 0
    std::shared_ptr<EdgeIntegral> edge_int_1d_low_order_;                          ///< Edge integral - 1D element, order = 0
    std::shared_ptr<EdgeIntegral> edge_int_2d_low_order_;                          ///< Edge integral - 2D element, order = 0
    std::shared_ptr<EdgeIntegral> edge_int_2d_high_order_;                         ///< Edge integral - 2D element, order = 2
    std::shared_ptr<CouplingIntegral> coupling_int_1d_2d_low_order_;               ///< Coupling integral - 1D-2D, order = 0
    RevertableList<PatchFETestBase::BulkIntegralData> bulk_integral_data_;         ///< Holds data for computing bulk integrals.
    RevertableList<PatchFETestBase::EdgeIntegralData> edge_integral_data_;         ///< Holds data for computing edge integrals.
    RevertableList<PatchFETestBase::CouplingIntegralData> coupling_integral_data_; ///< Holds data for computing edge integrals.

    ElQ<Scalar> det_1d_low_order_;
    ElQ<Scalar> det_1d_high_order_;
    ElQ<Scalar> det_2d_low_order_;

//    FeQArray<Vector> grad_scalar_shape_1d_;
//    FeQArray<Vector> grad_scalar_shape_2d_;
//    FeQArray<Vector> grad_scalar_shape_3d_;
};

/// Complete test with FE scalar operations
void eval_func_scalar(Mesh* mesh, unsigned int quad_order) {
	MixedPtr<FE_P_disc> fe(quad_order);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    PatchFEMultiQuadTest patch_fe(dh);
    patch_fe.test_evaluation();
//    patch_fe.reset();
//    patch_fe.test_evaluation();
}

TEST(PatchFeTest, new_op_dependency) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    std::string input_str = "{ mesh_file=\"mesh/simplest_square.msh\", optimize_mesh=false }";
    Mesh* mesh = mesh_full_constructor(input_str);

    // two tests with different quad_order and Scalar / Vector FE operations
    eval_func_scalar(mesh, 1);
    eval_func_scalar(mesh, 2);
}

