/*
 * observe.cc
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#include <string>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <queue>

#include "system/global_defs.h"
#include "input/accessors.hh"
#include "input/input_type.hh"
#include "system/armadillo_tools.hh"

#include "mesh/mesh.h"
#include "mesh/bih_tree.hh"
#include "mesh/region.hh"
#include "mesh/accessors.hh"
#include "io/observe.hh"
#include "io/element_data_cache.hh"
#include "fem/mapping_p1.hh"
#include "tools/unit_si.hh"


namespace IT = Input::Type;


/**
 * Helper class allows to work with ObservePoint above elements with different dimensions
 *
 * Allows:
 *  - calculate projection of points by dimension
 *  - snap to subelement
 */
template<unsigned int dim>
class ProjectionHandler {
public:
	/// Constructor
	ProjectionHandler() {};

	ObservePointData projection(arma::vec3 input_point, unsigned int i_elm, ElementAccessor<3> elm) {
		arma::mat::fixed<3, dim+1> elm_map = MappingP1<dim,3>::element_map(elm);
		arma::vec::fixed<dim+1> projection = MappingP1<dim,3>::project_real_to_unit(input_point, elm_map);
		projection = MappingP1<dim,3>::clip_to_element(projection);

		ObservePointData data;
		data.element_idx_ = i_elm;
		data.local_coords_ = projection.rows(1, elm.dim());
		data.global_coords_ = elm_map*projection;
		data.distance_ = arma::norm(data.global_coords_ - input_point, 2);
		data.proc_ = elm.proc();

		return data;
	}

    /**
     * Snap local coords to the subelement. Called by the ObservePoint::snap method.
     */
	void snap_to_subelement(ObservePointData & observe_data, ElementAccessor<3> elm, unsigned int snap_dim)
	{
		if (snap_dim <= dim) {
			double min_dist = 2.0; // on the ref element the max distance should be about 1.0, smaler then 2.0
			arma::vec min_center;
			for(auto &center : RefElement<dim>::centers_of_subelements(snap_dim))
			{
				double dist = arma::norm(center - observe_data.local_coords_, 2);
				if ( dist < min_dist) {
					min_dist = dist;
					min_center = center;
				}
			}
			observe_data.local_coords_ = min_center;
		}

		arma::mat::fixed<3, dim+1> elm_map = MappingP1<dim,3>::element_map(elm);
        observe_data.global_coords_ =  elm_map * RefElement<dim>::local_to_bary(observe_data.local_coords_);
	}

};

template class ProjectionHandler<1>;
template class ProjectionHandler<2>;
template class ProjectionHandler<3>;


/**
 * Helper struct, used as comparator of priority queue in ObservePoint::find_observe_point.
 */
struct CompareByDist
{
  bool operator()(const ObservePointData& lhs, const ObservePointData& rhs) const
  {
    return lhs.distance_ > rhs.distance_;
  }
};


/*******************************************************************
 * implementation of ObservePoint
 */

const Input::Type::Record & ObservePoint::get_input_type() {
    return IT::Record("ObservePoint", "Specification of the observation point.\n"
            "The actual observation element and the observation point on it is determined as follows:\n\n"
            "1. Find an initial element containing the initial point. If no such element exists, we report an error.\n"
            "2. Use BFS (Breadth-first search) starting from the inital element to find the 'observe element'. The observe element is the closest element.\n"
            "3. Find the closest projection of the inital point on the observe element and snap this projection according to the ``snap_dim``.\n")
        .allow_auto_conversion("point")
        .declare_key("name", IT::String(),
                IT::Default::read_time(
                        "Default name have the form 'obs_<id>', where 'id' "
                        "is the rank of the point on the input."),
                "Optional point name, which has to be unique.\n"
                "Any string that is a valid YAML key in record without any quoting can be used, however, "
                "using just alpha-numerical characters, and underscore instead of the space, is recommended."
                )
        .declare_key("point", IT::Array( IT::Double(), 3, 3 ), IT::Default::obligatory(),
                "Initial point for the observe point search.")
        .declare_key("snap_dim", IT::Integer(0, 4), IT::Default("4"),
                "The dimension of the sub-element to which center we snap. For value 4 no snapping is done. "
                "For values 0 up to 3 the element containing the initial point is found and then the observe"
                "point is snapped to the nearest center of the sub-element of the given dimension. "
                "E.g. for dimension 2 we snap to the nearest center of the face of the initial element."
                )
        .declare_key("snap_region", IT::String(), IT::Default("\"ALL\""),
                "The region of the initial element for snapping. Without snapping we make a projection to the initial element.")
        .declare_key("search_radius", IT::Double(0.0),
                IT::Default::read_time("Maximal distance of the observe point from the mesh relative to the mesh diameter. "),
                "Global value is defined in mesh record by the key global_snap_radius.")
        .close();
}

ObservePoint::ObservePoint()
{}


ObservePoint::ObservePoint(Input::Record in_rec, Mesh &mesh, unsigned int point_idx)
{
    in_rec_ = in_rec;

    string default_label = string("obs_") + std::to_string(point_idx);
    name_ = in_rec.val<string>("name", default_label );

    vector<double> tmp_coords;
    in_rec.val<Input::Array>("point").copy_to(tmp_coords);
    input_point_= arma::vec(tmp_coords);

    snap_dim_ = in_rec.val<unsigned int>("snap_dim");

    snap_region_name_ = in_rec.val<string>("snap_region");

    const BoundingBox &main_box = mesh.get_bih_tree().tree_box();
    double max_mesh_size = arma::max(main_box.max() - main_box.min());
    max_search_radius_ = in_rec_.val<double>("search_radius", mesh.global_snap_radius()) * max_mesh_size;
}



bool ObservePoint::have_observe_element() {
    return observe_data_.distance_ < numeric_limits<double>::infinity();
}



void ObservePoint::snap(Mesh &mesh)
{
    ElementAccessor<3> elm = mesh.element_accessor(observe_data_.element_idx_);
    switch (elm.dim()) {
        case 1:
        {
        	ProjectionHandler<1> ph;
        	ph.snap_to_subelement(observe_data_, elm, snap_dim_);
            break;
        }
        case 2:
        {
        	ProjectionHandler<2> ph;
        	ph.snap_to_subelement(observe_data_, elm, snap_dim_);
            break;
        }
        case 3:
        {
        	ProjectionHandler<3> ph;
        	ph.snap_to_subelement(observe_data_, elm, snap_dim_);
            break;
        }
        default: ASSERT(false).error("Clipping supported only for dim=1,2,3.");
    }
}



void ObservePoint::find_observe_point(Mesh &mesh) {
    RegionSet region_set = mesh.region_db().get_region_set(snap_region_name_);
    if (region_set.size() == 0)
        THROW( RegionDB::ExcUnknownSet() << RegionDB::EI_Label(snap_region_name_) << in_rec_.ei_address() );


    const BIHTree &bih_tree=mesh.get_bih_tree();
    vector<unsigned int> candidate_list;
    std::unordered_set<unsigned int> closed_elements(1023);
    std::priority_queue< ObservePointData, std::vector<ObservePointData>, CompareByDist > candidate_queue;

    // search for the initial element
    auto projected_point = bih_tree.tree_box().project_point(input_point_);
    bih_tree.find_point(projected_point, candidate_list, true);

    // closest element
    ObservePointData min_observe_point_data;
    
    for (unsigned int i_candidate=0; i_candidate<candidate_list.size(); ++i_candidate) {
        unsigned int i_elm=candidate_list[i_candidate];
        ElementAccessor<3> elm = mesh.element_accessor(i_elm);

        // project point, add candidate to queue
        auto observe_data = point_projection( i_elm, elm );
        
        // save the closest element for later diagnostic
        if(observe_data.distance_ < min_observe_point_data.distance_)
            min_observe_point_data = observe_data;
        
        // queue only the elements in the maximal search radius
        if (observe_data.distance_ <= max_search_radius_)
        	candidate_queue.push(observe_data);
        closed_elements.insert(i_elm);
    }

    // no candidates found -> exception
    if (candidate_queue.empty()) {
        THROW(ExcNoObserveElementCandidates()
            << EI_PointName(name_)
            << EI_Point(input_point_)
            << EI_ClosestEle(min_observe_point_data));
    }
    
    while (!candidate_queue.empty())
    {
        auto candidate_data = candidate_queue.top();
        candidate_queue.pop();

        unsigned int i_elm=candidate_data.element_idx_;
        ElementAccessor<3> elm = mesh.element_accessor(i_elm);

        // test if candidate is in region and update projection
        if (elm.region().is_in_region_set(region_set)) {
            ASSERT_LE(candidate_data.distance_, observe_data_.distance_).error();

			observe_data_.distance_ = candidate_data.distance_;
			observe_data_.element_idx_ = candidate_data.element_idx_;
			observe_data_.local_coords_ = candidate_data.local_coords_;
			observe_data_.global_coords_ = candidate_data.global_coords_;
			observe_data_.proc_ = candidate_data.proc_;
            break;
        }

        // add candidates to queue
		for (unsigned int n=0; n < elm->n_nodes(); n++)
			for(unsigned int i_node_ele : mesh.node_elements()[elm.node(n).idx()]) {
				if (closed_elements.find(i_node_ele) == closed_elements.end()) {
					ElementAccessor<3> neighbor_elm = mesh.element_accessor(i_node_ele);
					auto observe_data = point_projection( i_node_ele, neighbor_elm );
			        if (observe_data.distance_ <= max_search_radius_)
			        	candidate_queue.push(observe_data);
			        closed_elements.insert(i_node_ele);
				}
			}
    }

    if (! have_observe_element()) {
        THROW(ExcNoObserveElement()
            << EI_RegionName(snap_region_name_)
            << EI_PointName(name_)
            << EI_Point(input_point_)
            << EI_ClosestEle(min_observe_point_data));
    }
    snap( mesh );
    ElementAccessor<3> elm = mesh.element_accessor(observe_data_.element_idx_);
    double dist = arma::norm(elm.centre() - input_point_, 2);
    double elm_norm = arma::norm(elm.bounding_box().max() - elm.bounding_box().min(), 2);
    if (dist > 2*elm_norm)
    	WarningOut().fmt("Observe point ({}) is too distant from the mesh.\n", name_);
}




void ObservePoint::output(ostream &out, unsigned int indent_spaces, unsigned int precision)
{
    out << setw(indent_spaces) << "" << "- name: " << name_ << endl;
    out << setw(indent_spaces) << "" << "  init_point: " << field_value_to_yaml(input_point_, precision) << endl;
    out << setw(indent_spaces) << "" << "  snap_dim: " << snap_dim_ << endl;
    out << setw(indent_spaces) << "" << "  snap_region: " << snap_region_name_ << endl;
    out << setw(indent_spaces) << "" << "  observe_point: " << field_value_to_yaml(observe_data_.global_coords_, precision) << endl;
}



ObservePointData ObservePoint::point_projection(unsigned int i_elm, ElementAccessor<3> elm) {
	switch (elm.dim()) {
	case 1:
	{
		ProjectionHandler<1> ph;
		return ph.projection(input_point_, i_elm, elm);
		break;
	}
	case 2:
	{
		ProjectionHandler<2> ph;
		return ph.projection(input_point_, i_elm, elm);
		break;
	}
	case 3:
	{
		ProjectionHandler<3> ph;
		return ph.projection(input_point_, i_elm, elm);
		break;
	}
	default:
		ASSERT(false).error("Invalid element dimension!");
	}

	return ObservePointData(); // Should not happen.
}




/*******************************************************************
 * implementation of Observe
 */

const unsigned int Observe::max_observe_value_time = 1000;


Observe::Observe(string observe_name, Mesh &mesh, Input::Array in_array, unsigned int precision, std::string unit_str)
: observe_name_(observe_name),
  precision_(precision),
  point_ds_(nullptr),
  observe_time_idx_(0)
{
    observe_values_time_.reserve(max_observe_value_time);
    observe_values_time_.push_back(numeric_limits<double>::signaling_NaN());

    unsigned int global_point_idx=0, local_point_idx=0;

    // in_rec is Output input record.
    for(auto it = in_array.begin<Input::Record>(); it != in_array.end(); ++it) {
        ObservePoint point(*it, mesh, points_.size());
        point.find_observe_point(mesh);
        point.observe_data_.global_idx_ = global_point_idx++;
        if (point.observe_data_.proc_ == mesh.get_el_ds()->myp()) {
        	point.observe_data_.local_idx_ = local_point_idx++;
        	point_4_loc_.push_back(point.observe_data_.global_idx_);
        }
        else
        	point.observe_data_.local_idx_ = -1;
        points_.push_back( point );
        observed_element_indices_.push_back(point.observe_data_.element_idx_);
    }
    // make local to global map, distribution
    point_ds_ = new Distribution(Observe::max_observe_value_time * point_4_loc_.size(), PETSC_COMM_WORLD);

    // make indices unique
    std::sort(observed_element_indices_.begin(), observed_element_indices_.end());
    auto last = std::unique(observed_element_indices_.begin(), observed_element_indices_.end());
    observed_element_indices_.erase(last, observed_element_indices_.end());

    time_unit_str_ = unit_str;
    time_unit_seconds_ = UnitSI().s().convert_unit_from(unit_str);

    if (points_.size() == 0) return;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    if (rank_==0) {
        FilePath observe_file_path(observe_name_ + "_observe.yaml", FilePath::output_file);
        try {
            observe_file_path.open_stream(observe_file_);
            //observe_file_.setf(std::ios::scientific);
            observe_file_.precision(this->precision_);

        } INPUT_CATCH(FilePath::ExcFileOpen, FilePath::EI_Address_String, in_array)
        output_header();
    }
}

Observe::~Observe() {
    if (points_.size()>0 && observe_field_values_.size()>0)
        flush_values();

    observe_file_.close();
    if (point_ds_!=nullptr) delete point_ds_;
}


template <typename T>
ElementDataCache<T> & Observe::prepare_compute_data(std::string field_name, double field_time, unsigned int n_rows,
		unsigned int n_cols)
{
    if ( std::isnan(observe_values_time_[observe_time_idx_]) )
        observe_values_time_[observe_time_idx_] = field_time / time_unit_seconds_;
    else
        ASSERT(fabs(field_time / time_unit_seconds_ - observe_values_time_[observe_time_idx_]) < 2*numeric_limits<double>::epsilon())
              (field_time)(observe_values_time_[observe_time_idx_]);

    OutputDataFieldMap::iterator it=observe_field_values_.find(field_name);
    if (it == observe_field_values_.end()) {
        observe_field_values_[field_name]
					= std::make_shared< ElementDataCache<T> >(field_name, n_rows * n_cols, point_ds_->lsize());
        it=observe_field_values_.find(field_name);
    }
    return dynamic_cast<ElementDataCache<T> &>(*(it->second));
}

// explicit instantiation of template method
#define OBSERVE_PREPARE_COMPUTE_DATA(TYPE) \
template ElementDataCache<TYPE> & Observe::prepare_compute_data<TYPE>(std::string field_name, double field_time, \
		unsigned int n_rows, unsigned int n_cols)

OBSERVE_PREPARE_COMPUTE_DATA(int);
OBSERVE_PREPARE_COMPUTE_DATA(unsigned int);
OBSERVE_PREPARE_COMPUTE_DATA(double);


void Observe::output_header() {
    unsigned int indent = 2;
    observe_file_ << "# Observation file: " << observe_name_ << endl;
    observe_file_ << "time_unit: " << time_unit_str_ << endl;
    observe_file_ << "time_unit_in_seconds: " << time_unit_seconds_ << endl;
    observe_file_ << "points:" << endl;
    for(auto &point : points_)
        point.output(observe_file_, indent, precision_);
    observe_file_ << "data:" << endl;   

}

void Observe::flush_values() {
    std::vector<LongIdx> local_to_global(Observe::max_observe_value_time*point_4_loc_.size());
    for (unsigned int i=0; i<Observe::max_observe_value_time; ++i)
        for (unsigned int j=0; j<point_4_loc_.size(); ++j) local_to_global[i*point_4_loc_.size()+j] = i*points_.size()+point_4_loc_[j];

    for(auto &field_data : observe_field_values_) {
        auto serial_data = field_data.second->gather(point_ds_, &(local_to_global[0]));
        if (rank_==0) field_data.second = serial_data;
    }

    if (rank_ == 0) {
        unsigned int indent = 2;
        DebugOut() << "Observe::output_time_frame WRITE\n";
        for (unsigned int i_time=0; i_time<observe_time_idx_; ++i_time) {
            observe_file_ << setw(indent) << "" << "- time: " << observe_values_time_[i_time] << endl;
            for(auto &field_data : observe_field_values_) {
                observe_file_ << setw(indent) << "" << "  " << field_data.second->field_input_name() << ": ";
                field_data.second->print_yaml_subarray(observe_file_, precision_, i_time*points_.size(), (i_time+1)*points_.size());
                observe_file_ << endl;
            }
        }
    }

    observe_values_time_.clear();
    observe_values_time_.reserve(max_observe_value_time);
    observe_values_time_.push_back(numeric_limits<double>::signaling_NaN());
    observe_time_idx_ = 0;
}

void Observe::output_time_frame(bool flush) {
    if (points_.size() == 0) return;
    
    if ( ! no_fields_warning ) {
        no_fields_warning=true;
        // check that observe fields are set
        if (std::isnan(observe_values_time_[observe_time_idx_])) {
            // first call and no fields
            ASSERT(observe_field_values_.size() == 0);
            WarningOut() << "No observe fields for the observation stream: " << observe_name_ << endl;
        }
    }
    
    if (std::isnan(observe_values_time_[observe_time_idx_])) {
        ASSERT(observe_field_values_.size() == 0);
        return;        
    }
    
    observe_time_idx_++;
    if ( (observe_time_idx_==Observe::max_observe_value_time) || flush ) {
    	flush_values();
    } else {
        observe_values_time_.push_back( numeric_limits<double>::signaling_NaN() );
    }

}

Range<ObservePointAccessor> Observe::local_range() const {
	auto bgn_it = make_iter<ObservePointAccessor>( ObservePointAccessor(this, 0) );
	auto end_it = make_iter<ObservePointAccessor>( ObservePointAccessor(this, point_4_loc_.size()) );
    return Range<ObservePointAccessor>(bgn_it, end_it);
}

