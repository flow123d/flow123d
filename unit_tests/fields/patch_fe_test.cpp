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


class PatchFETest {
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


    PatchFETest(std::shared_ptr<DOFHandlerMultiDim> dh)
    : dh_(dh), fe_values_(CacheMapElementNumber::get(), dh_->ds()->fe()),
	  quad_({new QGauss(0, 2), new QGauss(1, 2), new QGauss(2, 2), new QGauss(3, 2)}),
	  bulk_integral_data_(20, 10),
	  edge_integral_data_(12, 6),
	  jac_det_1d_( this->fe_values_.determinant(this->quad_[1]) ),
	  jac_det_2d_( this->fe_values_.determinant(this->quad_[2]) ),
	  jac_det_3d_( this->fe_values_.determinant(this->quad_[3]) ),
	  jac_det_side_1d_( this->fe_values_.determinant_side(this->quad_[0]) ),
	  jac_det_side_2d_( this->fe_values_.determinant_side(this->quad_[1]) ),
	  jac_det_side_3d_( this->fe_values_.determinant_side(this->quad_[2]) )
    {
        eval_points_ = std::make_shared<EvalPoints>();
        // first step - create integrals, then - initialize cache and initialize PatchFEValues on all dimensions
        this->create_integrals();
        element_cache_map_.init(eval_points_);
        this->initialize();
    }

    ~PatchFETest() {}

    void create_integrals() {
        bulk_integrals_[0] = eval_points_->add_bulk<1>(*quad_[1]);
        bulk_integrals_[1] = eval_points_->add_bulk<2>(*quad_[2]);
        bulk_integrals_[2] = eval_points_->add_bulk<3>(*quad_[3]);
        edge_integrals_[0] = eval_points_->add_edge<1>(*quad_[0]);
        edge_integrals_[1] = eval_points_->add_edge<2>(*quad_[1]);
        edge_integrals_[2] = eval_points_->add_edge<3>(*quad_[2]);
    }

    void initialize() {
        UpdateFlags u = update_values | update_JxW_values | update_quadrature_points;
        this->fe_values_.initialize<1>(*this->quad_[1], u);
        this->fe_values_.initialize<2>(*this->quad_[2], u);
        this->fe_values_.initialize<3>(*this->quad_[3], u);
        this->fe_values_.initialize<1>(*this->quad_[0], u);
        this->fe_values_.initialize<2>(*this->quad_[1], u);
        this->fe_values_.initialize<3>(*this->quad_[2], u);
    }

    /// Return BulkPoint range of appropriate dimension
    inline Range< BulkPoint > bulk_points(unsigned int dim, unsigned int element_patch_idx) const {
        return bulk_integrals_[dim-1]->points(element_patch_idx, &element_cache_map_);
    }

    /// Return EdgePoint range of appropriate dimension
    inline Range< EdgePoint > edge_points(unsigned int dim, const DHCellSide &cell_side) const {
	    return edge_integrals_[dim-1]->points(cell_side, &element_cache_map_);
    }

    void add_integrals(DHCellAccessor cell) {
        // Bulk integral
        uint subset_idx = bulk_integrals_[cell.dim()-1]->get_subset_idx();
        bulk_integral_data_.emplace_back(cell, subset_idx);

        unsigned int reg_idx = cell.elm().region_idx().idx();
        // Different access than in other integrals: We can't use range method CellIntegral::points
        // because it passes element_patch_idx as argument that is not known during patch construction.
        for (uint i=uint( eval_points_->subset_begin(cell.dim(), subset_idx) );
                  i<uint( eval_points_->subset_end(cell.dim(), subset_idx) ); ++i) {
            element_cache_map_.add_eval_point(reg_idx, cell.elm_idx(), i, cell.local_idx());
        }

        for( DHCellSide cell_side : cell.side_range() ) {
            if ( (cell_side.n_edge_sides() >= 2) && (cell_side.edge_sides().begin()->element().idx() == cell.elm_idx())) {
                auto range = cell_side.edge_sides();
                edge_integral_data_.emplace_back(range, edge_integrals_[range.begin()->dim()-1]->get_subset_idx());

                for( DHCellSide edge_side : range ) {
                    unsigned int reg_idx = edge_side.element().region_idx().idx();
                    for (auto p : edge_integrals_[range.begin()->dim()-1]->points(edge_side, &element_cache_map_) ) {
                        element_cache_map_.add_eval_point(reg_idx, edge_side.elem_idx(), p.eval_point_idx(), edge_side.cell().local_idx());
                    }
                }
            }
        }
    }

    void update_patch() {
        for (unsigned int i=0; i<bulk_integral_data_.permanent_size(); ++i) {
            uint dim = bulk_integral_data_[i].cell.dim();
            uint element_patch_idx = element_cache_map_.position_in_cache(bulk_integral_data_[i].cell.elm_idx());
            uint elm_pos = fe_values_.register_element(bulk_integral_data_[i].cell, element_patch_idx);
            for (auto p : this->bulk_points(dim, element_patch_idx) ) {
                unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                fe_values_.register_point(bulk_integral_data_[i].cell, elm_pos, value_cache_idx);
            }
        }
        for (unsigned int i=0; i<edge_integral_data_.permanent_size(); ++i) {
        	auto range = edge_integral_data_[i].edge_side_range;
            uint dim = range.begin()->dim();
            for( DHCellSide edge_side : range )
            {
                uint element_patch_idx = element_cache_map_.position_in_cache(edge_side.cell().elm_idx());
            	uint elm_pos = fe_values_.register_element(edge_side.cell(), element_patch_idx, 1);
                for (auto p : this->edge_points(dim, edge_side) ) {
            	    unsigned int value_cache_idx = p.elm_cache_map()->element_eval_point(p.elem_patch_idx(), p.eval_point_idx());
                    fe_values_.register_point(edge_side.cell(), elm_pos, value_cache_idx, 1);
                }
            }
        }
        fe_values_.reinit_patch();
    }


    std::shared_ptr<DOFHandlerMultiDim> dh_;
    PatchFEValues<3> fe_values_;                                  ///< Common FEValues object over all dimensions
    std::shared_ptr<EvalPoints> eval_points_;                     ///< EvalPoints object shared by all integrals
    ElementCacheMap element_cache_map_;                           ///< ElementCacheMap according to EvalPoints
    std::array<std::shared_ptr<BulkIntegral>, 3> bulk_integrals_; ///< Bulk integrals of dim 1,2,3
    std::array<std::shared_ptr<EdgeIntegral>, 3> edge_integrals_; ///< Edge integrals of dim 1,2,3
    std::array<Quadrature*, 4> quad_;                             ///< Quadratures of dim 0,1,2,3
    RevertableList<BulkIntegralData> bulk_integral_data_;         ///< Holds data for computing bulk integrals.
    RevertableList<EdgeIntegralData> edge_integral_data_;         ///< Holds data for computing edge integrals.
    ElQ<Scalar> jac_det_1d_;
    ElQ<Scalar> jac_det_2d_;
    ElQ<Scalar> jac_det_3d_;
    ElQ<Scalar> jac_det_side_1d_;
    ElQ<Scalar> jac_det_side_2d_;
    ElQ<Scalar> jac_det_side_3d_;
};

TEST(PatchFeTest, bulk_points) {
    FilePath::set_io_dirs(".",UNIT_TESTS_SRC_DIR,"",".");
    Profiler::instance();
    PetscInitialize(0,PETSC_NULL,PETSC_NULL,PETSC_NULL);

    std::string input_str = "{ mesh_file=\"mesh/simplest_cube.msh\", optimize_mesh=false }";
    Mesh* mesh = mesh_full_constructor(input_str);

    MixedPtr<FE_P_disc> fe(0);
    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh, fe);
    std::shared_ptr<DOFHandlerMultiDim> dh = std::make_shared<DOFHandlerMultiDim>(*mesh);
    dh->distribute_dofs(ds);

    PatchFETest patch_fe(dh);
    for(auto cell_it = dh->local_range().begin(); cell_it != dh->local_range().end(); ++cell_it) {
    	patch_fe.add_integrals(*cell_it);
    }
    patch_fe.bulk_integral_data_.make_permanent();
    patch_fe.edge_integral_data_.make_permanent();
    patch_fe.element_cache_map_.make_paermanent_eval_points();
    patch_fe.element_cache_map_.create_patch(); // simplest_cube.msh contains 4 bulk regions, 9 bulk elements and 32 bulk points
    patch_fe.update_patch();

    patch_fe.fe_values_.print(true, true, true, false);
}
