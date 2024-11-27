#define TEST_USE_MPI
#define FEAL_OVERRIDE_ASSERTS
#include <flow_gtest_mpi.hh>
#include <mesh_constructor.hh>
#include "arma_expect.hh"

#include "quadrature/quadrature_lib.hh"
#include "fields/eval_subset.hh"
#include "fields/eval_points.hh"
#include "fields/field_value_cache.hh"
#include "fem/fe_values.hh"
#include "fem/patch_fe_values.hh"
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
    : dh_(dh), patch_fe_values_(quad_order, dh_->ds()->fe()),
	  fe_(quad_order), fe_values_(3), fe_values_side_(3),
	  bulk_integral_data_(20, 10),
	  edge_integral_data_(12, 6),
	  coupling_integral_data_(12, 6),
	  det_1d_( this->patch_fe_values_.bulk_values<1>().determinant() ),
	  det_2d_( this->patch_fe_values_.bulk_values<2>().determinant() ),
	  det_3d_( this->patch_fe_values_.bulk_values<3>().determinant() ),
	  jxw_1d_( this->patch_fe_values_.bulk_values<1>().JxW() ),
	  jxw_2d_( this->patch_fe_values_.bulk_values<2>().JxW() ),
	  jxw_3d_( this->patch_fe_values_.bulk_values<3>().JxW() ),
	  jxw_side_1d_( this->patch_fe_values_.side_values<1>().JxW() ),
	  jxw_side_2d_( this->patch_fe_values_.side_values<2>().JxW() ),
	  jxw_side_3d_( this->patch_fe_values_.side_values<3>().JxW() ),
	  normal_vec_1d_( this->patch_fe_values_.side_values<1>().normal_vector() ),
	  normal_vec_2d_( this->patch_fe_values_.side_values<2>().normal_vector() ),
	  normal_vec_3d_( this->patch_fe_values_.side_values<3>().normal_vector() ),
	  scalar_shape_1d_( this->patch_fe_values_.bulk_values<1>().scalar_shape() ),
	  scalar_shape_2d_( this->patch_fe_values_.bulk_values<2>().scalar_shape() ),
	  scalar_shape_3d_( this->patch_fe_values_.bulk_values<3>().scalar_shape() ),
	  scalar_shape_side_1d_( this->patch_fe_values_.side_values<1>().scalar_shape() ),
	  scalar_shape_side_2d_( this->patch_fe_values_.side_values<2>().scalar_shape() ),
	  scalar_shape_side_3d_( this->patch_fe_values_.side_values<3>().scalar_shape() ),
	  conc_join_shape_2d_( FeQJoin<Scalar>( this->patch_fe_values_.template join_values<2>().scalar_join_shape() ) ),
	  conc_join_shape_3d_( FeQJoin<Scalar>( this->patch_fe_values_.template join_values<3>().scalar_join_shape() ) )
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize PatchFEValues on all dimensions
        this->create_integrals();
        element_cache_map_.init(eval_points_);

        UpdateFlags u = update_values | update_inverse_jacobians | update_JxW_values | update_quadrature_points | update_volume_elements | update_gradients;
        UpdateFlags u_side = update_values | update_inverse_jacobians | update_side_JxW_values | update_normal_vectors | update_quadrature_points | update_gradients;
        fe_values_[0].initialize(*patch_fe_values_.get_quadrature(1, true), *fe_[Dim<1>{}], u);
        fe_values_[1].initialize(*patch_fe_values_.get_quadrature(2, true), *fe_[Dim<2>{}], u);
        fe_values_[2].initialize(*patch_fe_values_.get_quadrature(3, true), *fe_[Dim<3>{}], u);
        fe_values_side_[0].initialize(*patch_fe_values_.get_quadrature(1, false), *fe_[Dim<1>{}], u_side);
        fe_values_side_[1].initialize(*patch_fe_values_.get_quadrature(2, false), *fe_[Dim<2>{}], u_side);
        fe_values_side_[2].initialize(*patch_fe_values_.get_quadrature(3, false), *fe_[Dim<3>{}], u_side);
    }

    ~PatchFETestBase() {}

    void create_integrals() {
        bulk_integrals_[0] = eval_points_->add_bulk<1>(*patch_fe_values_.get_quadrature(1,true));
        bulk_integrals_[1] = eval_points_->add_bulk<2>(*patch_fe_values_.get_quadrature(2,true));
        bulk_integrals_[2] = eval_points_->add_bulk<3>(*patch_fe_values_.get_quadrature(3,true));
        edge_integrals_[0] = eval_points_->add_edge<1>(*patch_fe_values_.get_quadrature(1,false));
        edge_integrals_[1] = eval_points_->add_edge<2>(*patch_fe_values_.get_quadrature(2,false));
        edge_integrals_[2] = eval_points_->add_edge<3>(*patch_fe_values_.get_quadrature(3,false));
        coupling_integrals_[0] = eval_points_->add_coupling<2>(*patch_fe_values_.get_quadrature(1,true));
        coupling_integrals_[1] = eval_points_->add_coupling<3>(*patch_fe_values_.get_quadrature(2,true));
    }

    void initialize() {
        this->patch_fe_values_.initialize<1>(*patch_fe_values_.get_quadrature(1,true));
        this->patch_fe_values_.initialize<2>(*patch_fe_values_.get_quadrature(2,true));
        this->patch_fe_values_.initialize<3>(*patch_fe_values_.get_quadrature(3,true));
        this->patch_fe_values_.initialize<1>(*patch_fe_values_.get_quadrature(1,false));
        this->patch_fe_values_.initialize<2>(*patch_fe_values_.get_quadrature(2,false));
        this->patch_fe_values_.initialize<3>(*patch_fe_values_.get_quadrature(3,false));
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
        table_sizes_.elem_sizes_[0][dim-1]++;
        table_sizes_.point_sizes_[0][dim-1] += eval_points_->subset_size(dim, subset_idx); // add rows for bulk points to table

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( eval_points_->subset_begin(cell.dim(), subset_idx) );
                  i<uint( eval_points_->subset_end(cell.dim(), subset_idx) ); ++i) {
            element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }

        // Edge integral
        for( DHCellSide cell_side : cell.side_range() ) {
            if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                auto range = cell_side.edge_sides();
                uint subset_idx = edge_integrals_[range.begin()->dim()-1]->get_subset_idx();
                edge_integral_data_.emplace_back(range, subset_idx);

                for( DHCellSide edge_side : range ) {
                    uint dim = edge_side.dim();
                    table_sizes_.elem_sizes_[1][dim-1]++;
                    table_sizes_.point_sizes_[1][dim-1] += eval_points_->subset_size(dim, subset_idx) / (dim+1); // add rows for side points to table
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : edge_integrals_[range.begin()->dim()-1]->points(edge_side, &element_cache_map_) ) {
                        element_cache_map_.add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    }
                }
            }
        }

        // Coupling integral
        bool add_low = true;
    	for( DHCellSide neighb_side : cell.neighb_sides() ) { // cell -> elm lower dim, neighb_side -> elm higher dim
            if (cell.dim() != neighb_side.dim()-1) continue;
            coupling_integral_data_.emplace_back(cell, coupling_integrals_[cell.dim()-1]->get_subset_low_idx(), neighb_side,
                    coupling_integrals_[cell.dim()-1]->get_subset_high_idx());
            table_sizes_.elem_sizes_[1][cell.dim()]++;
            if (add_low) table_sizes_.elem_sizes_[0][cell.dim()-1]++;

            unsigned int reg_idx_low = cell.elm().region_idx().idx();
            unsigned int reg_idx_high = neighb_side.element().region_idx().idx();
            for (auto p : coupling_integrals_[cell.dim()-1]->points(neighb_side, &element_cache_map_) ) {
                element_cache_map_.add_eval_point(reg_idx_high, neighb_side.elem_idx(), p.eval_point_idx(), neighb_side.cell().local_idx());
                table_sizes_.point_sizes_[1][cell.dim()]++;

                if (add_low) {
                    auto p_low = p.lower_dim(cell); // equivalent point on low dim cell
                    element_cache_map_.add_eval_point(reg_idx_low, cell.elm_idx(), p_low.eval_point_idx(), cell.local_idx());
                    table_sizes_.point_sizes_[0][cell.dim()-1]++;
                }
            }
            add_low = false;
        }
    }

    void update_patch() {
        patch_fe_values_.resize_tables(table_sizes_);
        for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            uint dim = bulk_integral_data_[i].cell.dim();
            uint element_patch_idx = element_cache_map_.position_in_cache(bulk_integral_data_[i].cell.elm_idx());
            uint elm_pos = patch_fe_values_.register_element(bulk_integral_data_[i].cell, element_patch_idx);
            uint i_point = 0;
            for (auto p : this->bulk_points(dim, element_patch_idx) ) {
                unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                patch_fe_values_.register_bulk_point(bulk_integral_data_[i].cell, elm_pos, value_cache_idx, i_point++);
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
                    unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                    patch_fe_values_.register_side_point(edge_side, side_pos, value_cache_idx, i_point++);
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
                unsigned int value_cache_idx = p_high.elm_cache_map()->element_eval_point(p_high.elem_patch_idx(), p_high.eval_point_idx());
                patch_fe_values_.register_side_point(coupling_integral_data_[i].side, side_pos, value_cache_idx, i_side_point++);
                if (coupling_integral_data_[i].cell.elm_idx() != last_element_idx) {
                    auto p_low = p_high.lower_dim(coupling_integral_data_[i].cell);
                    value_cache_idx = p_low.elm_cache_map()->element_eval_point(p_low.elem_patch_idx(), p_low.eval_point_idx());
                    patch_fe_values_.register_bulk_point(coupling_integral_data_[i].cell, elm_pos, value_cache_idx, i_bulk_point++);
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
        table_sizes_.reset();
        patch_fe_values_.reset();
    }


    std::shared_ptr<DOFHandlerMultiDim> dh_;
    PatchFEValues<3> patch_fe_values_;                                     ///< Common FEValues object over all dimensions

    MixedPtr<FE_P_disc> fe_;
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
    FeQArray<Scalar> scalar_shape_1d_;
    FeQArray<Scalar> scalar_shape_2d_;
    FeQArray<Scalar> scalar_shape_3d_;
    FeQArray<Scalar> scalar_shape_side_1d_;
    FeQArray<Scalar> scalar_shape_side_2d_;
    FeQArray<Scalar> scalar_shape_side_3d_;
    FeQJoin<Scalar> conc_join_shape_2d_;
    FeQJoin<Scalar> conc_join_shape_3d_;

    /**
     * Struct for pre-computing number of elements, sides, bulk points and side points on each dimension.
     * Format:
     *  { {n_bulk_points_1D, 2D, 3D },
     *    {n_side_points_1D, 2D, 3D } }
     *
     * Passes its to PatchFEValues and sets size of tables in this object
     */
    PatchFEValues<3>::TableSizes table_sizes_;
};


class PatchFETestFull : public PatchFETestBase {
public:
	PatchFETestFull(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : PatchFETestBase(quad_order, dh),
	  grad_scalar_shape_1d_( this->patch_fe_values_.bulk_values<1>().grad_scalar_shape() ),
	  grad_scalar_shape_2d_( this->patch_fe_values_.bulk_values<2>().grad_scalar_shape() ),
	  grad_scalar_shape_3d_( this->patch_fe_values_.bulk_values<3>().grad_scalar_shape() ),
	  grad_scalar_shape_side_1d_( this->patch_fe_values_.side_values<1>().grad_scalar_shape() ),
	  grad_scalar_shape_side_2d_( this->patch_fe_values_.side_values<2>().grad_scalar_shape() ),
	  grad_scalar_shape_side_3d_( this->patch_fe_values_.side_values<3>().grad_scalar_shape() )
    {}

    ~PatchFETestFull() {}

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
            std::stringstream ss1, ss2;
            patch_fe_values_.print_data_tables(ss1, true, false, false);
            WarningOut() << ss1.str();
            patch_fe_values_.print_operations(ss2);
            WarningOut() << ss2.str();
        }

        for(auto dh_cell : dh_->local_range() ) {
            ElementAccessor<3> elm = dh_cell.elm();
            auto p = *( bulk_integrals_[dh_cell.dim()-1]->points(element_cache_map_.position_in_cache(dh_cell.elm_idx()), &element_cache_map_).begin() );
            double jxw = 0.0, jxw_ref = 0.0;
            double det = 0.0, det_ref = 0.0;
            double scalar_shape = 0.0, scalar_shape_ref = 0.0;
            arma::vec3 grad_scalar_dof0("0 0 0");
            arma::vec3 grad_scalar_dof0_ref("0 0 0");
            arma::vec3 grad_scalar_dof1("0 0 0");
            arma::vec3 grad_scalar_dof1_ref("0 0 0");
            switch (dh_cell.dim()) {
            case 1:
                fe_values_[0].reinit(elm);
                jxw = jxw_1d_(p);
                det = det_1d_(p);
                scalar_shape = scalar_shape_1d_.shape(0)(p);
            	grad_scalar_dof0 = grad_scalar_shape_1d_.shape(0)(p);
                jxw_ref = fe_values_[0].JxW(0);
                det_ref = fe_values_[0].determinant(0);
                scalar_shape_ref = fe_values_[0].shape_value(0, 0);
                grad_scalar_dof0_ref = fe_values_[0].shape_grad(0, 0);
                break;
            case 2:
                fe_values_[1].reinit(elm);
                jxw = jxw_2d_(p);
                det = det_2d_(p);
                scalar_shape = scalar_shape_2d_.shape(0)(p);
            	grad_scalar_dof0 = grad_scalar_shape_2d_.shape(0)(p);
            	grad_scalar_dof1 = grad_scalar_shape_2d_.shape(1)(p);
                jxw_ref = fe_values_[1].JxW(0);
                det_ref = fe_values_[1].determinant(0);
                scalar_shape_ref = fe_values_[1].shape_value(0, 0);
                grad_scalar_dof0_ref = fe_values_[1].shape_grad(0, 0);
                grad_scalar_dof1_ref = fe_values_[1].shape_grad(1, 0);
                break;
            case 3:
                fe_values_[2].reinit(elm);
                jxw = jxw_3d_(p);
                det = det_3d_(p);
                scalar_shape = scalar_shape_3d_.shape(0)(p);
            	grad_scalar_dof0 = grad_scalar_shape_3d_.shape(0)(p);
            	grad_scalar_dof1 = grad_scalar_shape_3d_.shape(1)(p);
                jxw_ref = fe_values_[2].JxW(0);
                det_ref = fe_values_[2].determinant(0);
                scalar_shape_ref = fe_values_[2].shape_value(0, 0);
                grad_scalar_dof0_ref = fe_values_[2].shape_grad(0, 0);
                grad_scalar_dof1_ref = fe_values_[2].shape_grad(1, 0);
                break;
            }
            EXPECT_DOUBLE_EQ( jxw, jxw_ref );
            EXPECT_DOUBLE_EQ( det, det_ref );
            EXPECT_DOUBLE_EQ( scalar_shape, scalar_shape_ref );
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
                    	i_dof_high++;
                    }
                    else {
                	    auto result = conc_join_shape_2d_.shape(i_dof)(p_low);
                	    auto ref = fe_values_[0].shape_value(i_dof_low, 0);
                	    EXPECT_DOUBLE_EQ( result, ref );
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
                    	i_dof_high++;
                	}
                	else {
                	    auto result = conc_join_shape_3d_.shape(i_dof)(p_low);
                	    auto ref = fe_values_[1].shape_value(i_dof_low, 0);
                	    EXPECT_DOUBLE_EQ( result, ref );
                    	i_dof_low++;
                	}
                }
                break;
            }
        }

    }

    FeQArray<Vector> grad_scalar_shape_1d_;
    FeQArray<Vector> grad_scalar_shape_2d_;
    FeQArray<Vector> grad_scalar_shape_3d_;
    FeQArray<Vector> grad_scalar_shape_side_1d_;
    FeQArray<Vector> grad_scalar_shape_side_2d_;
    FeQArray<Vector> grad_scalar_shape_side_3d_;
};


class PatchFETestCompare : public PatchFETestBase {
public:
	/// Number of repeats
	static const unsigned int n_loops = 1e05;

	PatchFETestCompare(unsigned int quad_order, std::shared_ptr<DOFHandlerMultiDim> dh)
    : PatchFETestBase(quad_order, dh)
    {}

    ~PatchFETestCompare() {}

    void reinit_patch_fe() override {
        START_TIMER("reinit_patch_fe");
        for (uint i=0; i<PatchFETestCompare::n_loops; ++i)
            patch_fe_values_.reinit_patch();
        END_TIMER("reinit_patch_fe");
    }

    void reinit_fe_values() {
        START_TIMER("reinit_fe_values");
        for (uint i=0; i<PatchFETestCompare::n_loops; ++i) {

            for(auto cell_it = this->dh_->local_range().begin(); cell_it != this->dh_->local_range().end(); ++cell_it) {
                ElementAccessor<3> elm = cell_it->elm();
                switch (elm.dim()) {
                case 1:
                    fe_values_[0].reinit(elm);
                    for (uint si=0; si<2; ++si)
                        fe_values_side_[0].reinit(*elm.side(si));
                    break;
                case 2:
                    fe_values_[1].reinit(elm);
                    for (uint si=0; si<3; ++si)
                        fe_values_side_[1].reinit(*elm.side(si));
                    break;
                case 3:
                    fe_values_[2].reinit(elm);
                    for (uint si=0; si<4; ++si)
                        fe_values_side_[2].reinit(*elm.side(si));
                    break;
                }
            }

        }
        END_TIMER("reinit_fe_values");
    }
};


void compare_evaluation_func(Mesh* mesh, unsigned int quad_order) {
    MixedPtr<FE_P_disc> fe(quad_order);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    PatchFETestFull patch_fe(quad_order, dh);
    patch_fe.initialize();
    patch_fe.test_evaluation(true);
    patch_fe.reset();
    patch_fe.test_evaluation();
}



TEST(PatchFeTest, complete_evaluation) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    std::string input_str = "{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }";
    Mesh* mesh = mesh_full_constructor(input_str);

    // two tests with different quad_order
    compare_evaluation_func(mesh, 1);
    compare_evaluation_func(mesh, 2);
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


//TEST(PatchFeTest, normalize) {
//
//    typedef Eigen::Array<double,Eigen::Dynamic,1>  ArrayDbl;
//
//    Eigen::Vector<ArrayDbl, Eigen::Dynamic> vector_table;
//    vector_table.resize(3);
//    eigen_tools::resize_table(vector_table, 4);
//
//    vector_table(0) << 3.0, 0.5, 0.0, 1.0;
//    vector_table(1) << 0.0, 0.5, 2.5, 2.0;
//    vector_table(2) << 4.0, 0.0, 0.0, 2.0;
//
//    ArrayDbl norm_vec;
//    norm_vec.resize(4);
//    Eigen::VectorXd A(3);
//
//    for (uint i=0; i<4; ++i) {
//        A(0) = vector_table(0)(i);
//        A(1) = vector_table(1)(i);
//        A(2) = vector_table(2)(i);
//        norm_vec(i) = A.norm();
//    }
//    for (uint i=0; i<3; ++i) {
//        vector_table(i) /= norm_vec;
//    }
////    std::cout << "norm_vec: " << norm_vec.transpose() << std::endl << std::endl;
//
//    for (uint i=0; i<4; ++i)
//        std::cout << "(" << vector_table(0)(i) << ", " << vector_table(1)(i) << ", " << vector_table(2)(i) << ")" << std::endl;
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
