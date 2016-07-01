/*
 * observe.cc
 *
 *  Created on: Jun 28, 2016
 *      Author: jb
 */

#include <string>
#include "system/global_defs.h"
#include "input/accessors.hh"
#include "input/input_type.hh"

#include "mesh/region.hh"

#include "io/observe.hh"
#include "io/output_data.hh"


namespace IT = Input::Type;


ObservePoint::ObservePoint(const string &name, const arma::vec3 &coords, unsigned int snap_dim, const std::string &snap_region_name)
: name_(name), input_point_(coords), snap_dim_(snap_dim), snap_region_name_(snap_region_name)
{}




const Input::Type::Record & ObservePoint::get_input_type() {
    return IT::Record("ObservePoint", "Specification of the observation point.")
        .declare_key("name", IT::String(),
                IT::Default::read_time(
                        "Default label have form 'obs_<id>', where 'id' "
                        "is the rank of the point at input."),
                "Optional point name. Has to be unique. Any string that is valid YAML key in record without any quoting can be used however"
                "using just alpha-numerical characters and underscore instead of the space is recommended. "
                )
        .declare_key("point", IT::Array( IT::Double(), 3, 3 ), IT::Default::obligatory(),
                "Position of the observation point before snapping (see snap_to_center_dim below).")
        .declare_key("snap_dim", IT::Integer(0, 4), IT::Default("4"),
                "The dimension of the sub-element to which center we snap. For value 4 no snapping is done. "
                "For values 0 up to 3 the element containing the initial point is found and then the observe"
                "point is snapped to the nearest center of the sub-element of the given dimension. "
                "E.g. for dimension 2 we snap to the nearest center of the face of the initial element."
                "Let us note that we do not search for the nearest center in larger neighborhood of the initial element,"
                "in particular for snap_to_center_dim=0 there could be a node which is closer "
                "to the initial point then the snapped observe point.")
        .declare_key("snap_region", IT::String(), IT::Default("\"ALL\""),
                "The region of the initial element for snapping. Without snapping we make a projection to the initial element.")
        .declare_key("snapping_levels", IT::Integer(1), IT::Default("1"),
                "Maximum number of levels of the breath first search used to find initial element from the element containing the given point.")
        .close();
}


Observe::Observe(string observe_name, Mesh &mesh, Input::Record in_rec)
: mesh_(&mesh)
{
    // in_rec is Output input record.


    auto observe_fields = in_rec.val<Input::Array>("observe_fields");
    for(auto it = observe_fields.begin<Input::FullEnum>(); it != observe_fields.end(); ++it ) {
        this->field_names_.insert((std::string)*it);
    }

    auto op_input_array = in_rec.val<Input::Array>("observe_points");
    for(auto it = op_input_array.begin<Input::Record>(); it != op_input_array.end(); ++it ) {

        string default_label = string("obs_") + std::to_string(points_.size());
        string label = it->val<string>("name", default_label );

        vector<double> tmp_coords;
        it->val<Input::Array>("point").copy_to(tmp_coords);
        arma::vec3 coords= arma::vec(tmp_coords);

        string region_name = it->val<string>("snap_region");
        RegionSet snap_region = mesh_->region_db().get_region_set(region_name);
        if (snap_region.size() == 0)
            THROW( RegionDB::ExcUnknownSet() << RegionDB::EI_Label(region_name) << it->ei_address() );

        points_.push_back(ObservePoint(
                label,
                coords,
                it->val<unsigned int>("snap_dim"),
                region_name));
    }
    find_observe_points();

    observe_file_.open((observe_name + "_observe.yaml").c_str());
}

Observe::~Observe() {
    observe_file_.close();
}

void Observe::find_observe_points() {

}

template<int spacedim, class Value>
void Observe::compute_field_values(Field<spacedim, Value> &field) {

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
template void Observe::compute_field_values(Field<dim, FieldValue<3>::VectorFixed> &); \
template void Observe::compute_field_values(Field<dim, FieldValue<3>::TensorFixed> &);

// Make all instances for both dimensions.
INSTANCE_DIM(2)
INSTANCE_DIM(3)




void Observe::output_time_frame() {

}
