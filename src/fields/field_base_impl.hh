/*
 * function_base_impl.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#ifndef FIELD_BASE_IMPL_HH_
#define FIELD_BASE_IMPL_HH_

#include <string>
using namespace std;

#include <boost/type_traits.hpp>
#include <boost/format.hpp>
#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>


#include "fields/field_base.hh"
#include "fields/field_interpolated_p0.hh"
#include "fields/field_python.hh"
#include "fields/field_constant.hh"
#include "fields/field_formula.hh"
#include "fields/field_elementwise.hh"

#include "fields/field_values.hh"

#include "input/input_type.hh"
#include "input/json_to_storage.hh"


namespace it = Input::Type;







// ************************************************************************************************** implementation of FieldBase<...>


template <int spacedim, class Value>
FieldBase<spacedim, Value>::FieldBase(unsigned int n_comp)
: time_(0.0), value_(r_value_)
{
    value_.set_n_comp(n_comp);
}




template <int spacedim, class Value>
string FieldBase<spacedim, Value>::template_name() {
    return boost::str(boost::format("R%i -> %s") % spacedim % Value::type_name() );
}


template <int spacedim, class Value>
it::AbstractRecord FieldBase<spacedim, Value>::input_type
    = it::AbstractRecord("Field:"+template_name(), "Abstract record for all time-space functions.")
          .allow_auto_conversion("FieldConstant");


template <int spacedim, class Value>
Input::Type::AbstractRecord FieldBase<spacedim, Value>::get_input_type(typename Value::ElementInputType *element_input_type) {
    it::AbstractRecord type= it::AbstractRecord("Field:"+template_name(), "Abstract record for all time-space functions.");
    type.allow_auto_conversion("FieldConstant");

    FieldConstant<spacedim,Value>::get_input_type(type, element_input_type);
    FieldFormula<spacedim,Value>::get_input_type(type, element_input_type);
#ifdef HAVE_PYTHON
    FieldPython<spacedim,Value>::get_input_type(type, element_input_type);
#endif
    //FieldInterpolatedP0<spacedim,Value>::get_input_type(type, element_input_type);
    FieldElementwise<spacedim,Value>::get_input_type(type, element_input_type);

    return type;
}



template <int spacedim, class Value>
boost::shared_ptr< FieldBase<spacedim, Value> >
FieldBase<spacedim, Value>::function_factory(const Input::AbstractRecord &rec, unsigned int n_comp )
{
    boost::shared_ptr< FieldBase<spacedim, Value> > func;

    if (rec.type() == FieldInterpolatedP0<spacedim,Value>::input_type ) {
//        func= new FieldInterpolatedP0<spacedim,Value>(n_comp);
#ifdef HAVE_PYTHON
    } else if (rec.type() == FieldPython<spacedim,Value>::input_type ) {
        func=boost::make_shared< FieldPython<spacedim, Value> >(n_comp);
#endif
    } else if (rec.type() == FieldConstant<spacedim, Value>::input_type ) {
        func=boost::make_shared< FieldConstant<spacedim,Value> >(n_comp);
    } else if (rec.type() == FieldFormula<spacedim,Value>::input_type ) {
        func=boost::make_shared< FieldFormula<spacedim,Value> >(n_comp);
    } else if (rec.type() == FieldElementwise<spacedim,Value>::input_type ) {
        func=boost::make_shared< FieldElementwise<spacedim,Value> >(n_comp);
    } else {
        xprintf(PrgErr,"TYPE of Field is out of set of descendants. SHOULD NOT HAPPEN.\n");
    }
    func->init_from_input(rec);
    return func;
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::init_from_input(const Input::Record &rec) {
    xprintf(PrgErr, "The field '%s' do not support initialization from input.\n",
            typeid(this).name());
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::set_time(double time) {
    time_ = time;
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::set_mesh(Mesh *mesh) {
}



template<int spacedim, class Value>
unsigned int FieldBase<spacedim, Value>::n_comp() const {
    return (Value::NRows_ ? 0 : value_.n_rows());
}


/*
template <int spacedim, class Value>
typename Value::return_type &FieldBase<spacedim, Value>::value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm) {
    this->value(p,elm, this->r_value_);
    return this->r_value_;
}
*/

/*
template <int spacedim, class Value>
void FieldBase<spacedim, Value>::value_list(const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list, std::vector<FieldResult> &result_list)
{
    ASSERT_SIZES( point_list.size(), value_list.size() );
    for(unsigned int i=0; i< point_list.size(); i++)
        value_list[i] = value(point_list[i], elm);
}
*/



// ************************************************************************************************** implementation of Field<...>


template<int spacedim, class Value>
Field<spacedim,Value>::Field()
: FieldCommonBase(false), no_check_control_field_(NULL)
{
    this->enum_valued_ = boost::is_same<typename Value::element_type, FieldEnum>::value;
}



template<int spacedim, class Value>
void Field<spacedim, Value>::disable_where(const Field<spacedim, typename FieldValue<spacedim>::Enum > *control_field, const vector<FieldEnum> &value_list) {
    no_check_control_field_=control_field;
    no_check_values_=value_list;
}



template<int spacedim, class Value>
boost::shared_ptr< typename Field<spacedim,Value>::FieldBaseType >
Field<spacedim,Value>::operator[] (Region reg)
{
    ASSERT_LESS(reg.idx(), this->region_fields_.size());
    return this->region_fields_[reg.idx()];
}



template<int spacedim, class Value>
void Field<spacedim, Value>::set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) {
    ASSERT( this->mesh_, "Null mesh pointer, set_mesh() has to be called before set_from_input().\n");
    // initialize table if it is empty, we assume that the RegionDB is closed at this moment
    if (region_fields_.size() == 0)
        region_fields_.resize( this->mesh_->region_db().size() );

    boost::shared_ptr<FieldBaseType> field = FieldBaseType::function_factory(rec, this->n_comp_);
    // following dosn't lead to a memory leak since we use shared_ptr,
    // if the previous pointer was the last one the pointed field is correctly freed
    BOOST_FOREACH(Region reg, domain) region_fields_[reg.idx()] = field;
}



template<int spacedim, class Value>
void Field<spacedim, Value>::set_field(const RegionSet &domain, boost::shared_ptr< FieldBaseType > field) {
    // initialize table if it is empty, we assume that the RegionDB is closed at this moment
    if (region_fields_.size() == 0)
        region_fields_.resize( this->mesh_->region_db().size() );

    ASSERT_SIZES( field->n_comp() , this->n_comp_);
    field->set_mesh( this->mesh_ );
    BOOST_FOREACH(Region reg, domain) region_fields_[reg.idx()] = field;
}


template<int spacedim, class Value>
void Field<spacedim, Value>::set_time(double time) {
    if (region_fields_.size() == 0)
        region_fields_.resize( this->mesh_->region_db().size() );

    // check there are no empty field pointers, collect regions to be initialized from default value
    RegionSet regions_to_init; // empty vector

    BOOST_FOREACH(const Region &reg, this->mesh_->region_db().get_region_set("ALL") )
        if (reg.is_boundary() == this->bc_) {      // for regions that match type of the field domain
            if (! region_fields_[reg.idx()] ) {    // null field ptr
                if (no_check_control_field_) {      // is the check turned off?
                    FieldEnum value;
                    if (no_check_control_field_->get_constant_enum_value(reg, value)
                        && ( std::find(no_check_values_.begin(), no_check_values_.end(), value)
                             != no_check_values_.end() )
                       ) continue;                  // the field is not needed on this region
                }

                if (this->default_.has_value_at_declaration()) {    // try to use default
                    regions_to_init.push_back( reg );
                } else {
                    xprintf(UsrErr, "Missing value of the field '%s' on region ID: %d label: %s.\n", name_.c_str(), reg.id(), reg.label().c_str() );
                }
            }
        }

    // possibly set from default value
    if ( regions_to_init.size() ) {
        Input::JSONToStorage reader;
        Input::Type::AbstractRecord a_rec_type = make_input_tree();
        reader.read_from_default(this->default_.value(), a_rec_type );
        set_from_input( regions_to_init, reader.get_root_interface<Input::AbstractRecord>() );
    }

    // set mesh and time
    BOOST_FOREACH(const Region &reg, this->mesh_->region_db().get_region_set("ALL") )
        if (reg.is_boundary() == this->bc_ && region_fields_[reg.idx()] ) {      // for regions that match type of the field domain
                                                                                 // NULL pointers are only on "no_check" regions
            region_fields_[reg.idx()]->set_mesh(this->mesh_);
            region_fields_[reg.idx()]->set_time(time);
        }

}


// helper functions
template<int spacedim, class FieldBaseType>
FieldEnum get_constant_enum_value_dispatch(boost::shared_ptr< FieldBaseType > region_field,  const boost::true_type&) {
    return region_field->value(Point<spacedim>(), ElementAccessor<spacedim>());
}

template<int spacedim,class FieldBaseType>
FieldEnum get_constant_enum_value_dispatch(boost::shared_ptr< FieldBaseType > region_field,  const boost::false_type&) {
    return 0;
}



template<int spacedim, class Value>
bool Field<spacedim,Value>::get_constant_enum_value(RegionIdx r_idx,  FieldEnum &value) const {
    if (boost::is_same<typename Value::return_type, FieldEnum>::value) {
        boost::shared_ptr< FieldBaseType > region_field = region_fields_[r_idx.idx()];
        if (region_field && typeid(*region_field) == typeid(FieldConstant<spacedim, Value>)) {
            value = get_constant_enum_value_dispatch<spacedim>(region_field, boost::is_same<typename Value::return_type, FieldEnum>() );
            return true;
        }
    }
    return false;
}


template<int spacedim, class Value>
FieldResult Field<spacedim,Value>::field_result( ElementAccessor<spacedim> &elm) const {
    boost::shared_ptr< FieldBaseType > f = region_fields_[elm.region().idx()];
    if (f) return f->field_result();
    else return result_none;
}



template<int spacedim, class Value>
inline typename Value::return_type const & Field<spacedim,Value>::value(const Point<spacedim> &p, const ElementAccessor<spacedim> &elm)  {
    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx out of range, field: %s\n", this->name_.c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] , "Null field ptr on region id: %d, field: %s\n", elm.region().id(), this->name_.c_str());
    return region_fields_[elm.region_idx().idx()]->value(p,elm);
}



template<int spacedim, class Value>
inline void Field<spacedim,Value>::value_list(const std::vector< Point<spacedim> >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    region_fields_[elm.region_idx().idx()]->value_list(point_list,elm, value_list);
}



#endif //FUNCTION_BASE_IMPL_HH_
