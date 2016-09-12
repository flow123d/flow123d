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

#include "system/global_defs.h"
#include "input/accessors.hh"
#include "input/input_type.hh"
#include "system/armadillo_tools.hh"

#include "mesh/mesh.h"
#include "mesh/bih_tree.hh"
#include "mesh/region.hh"
#include "io/observe.hh"
#include "io/output_data.hh"


namespace IT = Input::Type;


const Input::Type::Record & ObservePoint::get_input_type() {
    return IT::Record("ObservePoint", "Specification of the observation point. The actual observe element and the observe point on it is determined as follows:\n\n"
            "1. Find an initial element containing the initial point. If no such element exists we report the error.\n"
            "2. Use BFS starting from the inital element to find the 'observe element'. The observe element is the closest element "
            "3. Find the closest projection of the inital point on the observe element and snap this projection according to the 'snap_dim'.\n")
        .allow_auto_conversion("point")
        .declare_key("name", IT::String(),
                IT::Default::read_time(
                        "Default name have the form 'obs_<id>', where 'id' "
                        "is the rank of the point on the input."),
                "Optional point name. Has to be unique. Any string that is valid YAML key in record without any quoting can be used however"
                "using just alpha-numerical characters and underscore instead of the space is recommended. "
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
        .declare_key("n_search_levels", IT::Integer(0), IT::Default("1"),
                "Maximum number of levels of the breadth first search used to find the observe element from the initial element. Value zero means to search only the initial element itself.")
        .close();
}

ObservePoint::ObservePoint()
{}


ObservePoint::ObservePoint(Input::Record in_rec, unsigned int point_idx)
:  distance_(numeric_limits<double>::infinity())
{
    in_rec_ = in_rec;

    string default_label = string("obs_") + std::to_string(point_idx);
    name_ = in_rec.val<string>("name", default_label );

    vector<double> tmp_coords;
    in_rec.val<Input::Array>("point").copy_to(tmp_coords);
    input_point_= arma::vec(tmp_coords);

    snap_dim_ = in_rec.val<unsigned int>("snap_dim");

    snap_region_name_ = in_rec.val<string>("snap_region");

    max_levels_ =  in_rec_.val<unsigned int>("n_search_levels");
}



void ObservePoint::update_projection(unsigned int i_elm, arma::vec local_coords, arma::vec3 global_coords)
{
    double dist = arma::norm(global_coords - input_point_, 2);
    //cout << "dist: " << dist << endl;
    if (dist < distance_) {
        distance_ = dist;
        element_idx_ = i_elm;
        local_coords_ = local_coords;
        global_coords_ = global_coords;
    }
}



bool ObservePoint::have_observe_element() {
    return distance_ < numeric_limits<double>::infinity();
}



template <int ele_dim>
void ObservePoint::snap_to_subelement()
{
    if (this->snap_dim_ >  ele_dim) return;

    double min_dist = 2.0; // on the ref element the max distance should be about 1.0, smaler then 2.0
    arma::vec min_center;
    for(auto &center : RefElement<ele_dim>::centers_of_subelements(this->snap_dim_))
    {
        double dist = arma::norm(center-local_coords_, 2);
        if ( dist < min_dist) {
            min_dist = dist;
            min_center = center;
        }
    }
    this->local_coords_ = min_center;
}


void ObservePoint::snap(Mesh &mesh)
{
    Element & elm = mesh.element[element_idx_];
    switch (elm.dim()) {
             case 1: snap_to_subelement<1>(); break;
             case 2: snap_to_subelement<2>(); break;
             case 3: snap_to_subelement<3>(); break;
             default: ASSERT(false).error("Clipping supported only for dim=1,2,3.");
    }
    this->global_coords_ =  elm.element_map() * arma::join_cols(this->local_coords_, arma::ones(1));
}



void ObservePoint::find_observe_point(Mesh &mesh) {
    RegionSet region_set = mesh.region_db().get_region_set(snap_region_name_);
    if (region_set.size() == 0)
        THROW( RegionDB::ExcUnknownSet() << RegionDB::EI_Label(snap_region_name_) << in_rec_.ei_address() );


    const BIHTree &bih_tree=mesh.get_bih_tree();
    vector<unsigned int> candidate_list, process_list;
    std::unordered_set<unsigned int> closed_elements(1023);

    // search for the initial element
    bih_tree.find_point(input_point_, candidate_list);
    process_list.swap(candidate_list);
    candidate_list.clear();
    for(unsigned int i_elm : process_list) {
        Element & elm = mesh.element[i_elm];
        arma::mat map = elm.element_map();
        arma::vec projection = elm.project_point(input_point_, map);

        // check that point is on the element
        if (projection.min() >=0.0) {
            // This is initial element.
            //input_point_.print(cout, "input_point");
            //cout << "i_el: " << i_elm << endl;
            //projection.print(cout, "projection");

            // if element match region filter store it as observe element to the obs. point
            if (elm.region().is_in_region_set(region_set)) {
                projection[elm.dim()] = 1.0; // use last coordinates for translation
                arma::vec global_coord = map*projection;
                update_projection(i_elm, projection.rows(0, elm.dim()-1), global_coord);
            }

            closed_elements.insert(i_elm);
            // add all node neighbours to the next level list
            for (unsigned int n=0; n < elm.n_nodes(); n++) {
                for(unsigned int i_node_ele : mesh.node_elements[mesh.node_vector.index(elm.node[n])])
                    candidate_list.push_back(i_node_ele);
            }
        }
    }

    if (candidate_list.size() == 0) THROW( ExcNoInitialPoint() << in_rec_.ei_address() );

    // Try to snap to the observe element with required snap_region
    for(unsigned int i_level=0; i_level < max_levels_; i_level++) {
        if (have_observe_element()) break;
        process_list.swap(candidate_list);
        candidate_list.clear();
        for(unsigned int i_elm : process_list) {
            if (closed_elements.find(i_elm) != closed_elements.end()) continue;
            Element & elm = mesh.element[i_elm];

            // if element match region filter, update the obs. point
            if (elm.region().is_in_region_set(region_set)) {
                arma::mat map = elm.element_map();
                arma::vec projection = elm.project_point(input_point_, map);
                arma::vec point_on_element = elm.clip_to_element(projection);

                point_on_element[elm.dim()] = 1.0; // use last coordinates for translation
                arma::vec global_coord = map*point_on_element;
                update_projection(i_elm, point_on_element.rows(0, elm.dim()-1), global_coord);
             }
            // add all node neighbours to the next level list
            for (unsigned int n=0; n < elm.n_nodes(); n++) {
                for(unsigned int i_node_ele : mesh.node_elements[mesh.node_vector.index(elm.node[n])])
                    candidate_list.push_back(i_node_ele);
            }
        }
    }
    if (! have_observe_element()) {
        THROW(ExcNoObserveElement() << EI_RegionName(snap_region_name_) << EI_NLevels(max_levels_) );
    }
    snap( mesh );
}



void ObservePoint::output(ostream &out, unsigned int indent_spaces, unsigned int precision)
{
    out << setw(indent_spaces) << "" << "- name: " << name_ << endl;
    out << setw(indent_spaces) << "" << "  init_point: " << field_value_to_yaml(input_point_) << endl;
    out << setw(indent_spaces) << "" << "  snap_dim: " << snap_dim_ << endl;
    out << setw(indent_spaces) << "" << "  snap_region: " << snap_region_name_ << endl;
    out << setw(indent_spaces) << "" << "  observe_point: " << field_value_to_yaml(global_coords_, precision) << endl;
}




Observe::Observe(string observe_name, Mesh &mesh, Input::Array in_array, unsigned int precision)
: mesh_(&mesh),
  observe_values_time_(numeric_limits<double>::signaling_NaN()),
  precision_(precision)
{
    // in_rec is Output input record.

    for(auto it = in_array.begin<Input::Record>(); it != in_array.end(); ++it) {
        ObservePoint point(*it, points_.size());
        point.find_observe_point(*mesh_);
        points_.push_back( point );
        observed_element_indices_.push_back(point.element_idx_);
    }
    // make indices unique
    std::sort(observed_element_indices_.begin(), observed_element_indices_.end());
    auto last = std::unique(observed_element_indices_.begin(), observed_element_indices_.end());
    observed_element_indices_.erase(last, observed_element_indices_.end());

    time_unit_str_ = "s";
    time_unit_seconds_ = 1.0;

    if (points_.size() == 0) return;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    if (rank_==0) {
        FilePath observe_file_path(observe_name + "_observe.yaml", FilePath::output_file);
        observe_file_path.open_stream(observe_file_);
        output_header(observe_name);
    }
}

Observe::~Observe() {
    observe_file_.close();
}


template<int spacedim, class Value>
void Observe::compute_field_values(Field<spacedim, Value> &field)
{
    if (points_.size() == 0) return;

    double field_time = field.time();
    if ( std::isnan(observe_values_time_) )
        observe_values_time_ = field_time;
    else
        ASSERT(fabs(field_time - observe_values_time_) < 2*numeric_limits<double>::epsilon())
              (field_time)(observe_values_time_);

    OutputDataFieldMap::iterator it=observe_field_values_.find(field.name());
    if (it == observe_field_values_.end()) {
        observe_field_values_[field.name()] = std::make_shared< OutputData<Value> >(field, points_.size());
        it=observe_field_values_.find(field.name());
    }
    OutputData<Value> &output_data = dynamic_cast<OutputData<Value> &>(*(it->second));

    unsigned int i_data=0;
    for(ObservePoint &o_point : points_) {
        unsigned int ele_index = o_point.element_idx_;
        const Value &obs_value =
                Value( const_cast<typename Value::return_type &>(
                        field.value(o_point.global_coords_,
                                ElementAccessor<spacedim>(this->mesh_, ele_index,false)) ));
        output_data.store_value(i_data,  obs_value);
        i_data++;
    }

}

// Instantiation of the method template for particular dimension.
#define INSTANCE_DIM(dim) \
template void Observe::compute_field_values(Field<dim, FieldValue<0>::Enum> &); \
template void Observe::compute_field_values(Field<dim, FieldValue<0>::Integer> &); \
template void Observe::compute_field_values(Field<dim, FieldValue<0>::Scalar> &); \
template void Observe::compute_field_values(Field<dim, FieldValue<dim>::VectorFixed> &); \
template void Observe::compute_field_values(Field<dim, FieldValue<dim>::TensorFixed> &);

// Make all instances for both dimensions.
INSTANCE_DIM(2)
INSTANCE_DIM(3)


void Observe::output_header(string observe_name) {
    unsigned int indent = 2;
    observe_file_ << "# Observation file: " << observe_name << endl;
    observe_file_ << "time_unit: " << time_unit_str_ << endl;
    observe_file_ << "time_unit_in_secodns: " << time_unit_seconds_ << endl;
    observe_file_ << "points:" << endl;
    for(auto &point : points_)
        point.output(observe_file_, indent, precision_);
    observe_file_ << "data:" << endl;

}

void Observe::output_time_frame(double time) {
    if (points_.size() == 0) return;

    if (rank_ == 0) {
        unsigned int indent = 2;
        observe_file_ << setw(indent) << "" << "- time: " << observe_values_time_ << endl;
        for(auto &field_data : observe_field_values_) {
            observe_file_ << setw(indent) << "" << "  " << field_data.second->field_name << ": ";
            field_data.second->print_all_yaml(observe_file_, precision_);
            observe_file_ << endl;
        }
    }

    observe_values_time_ = numeric_limits<double>::signaling_NaN();

}
