/*
 * function_base_impl.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */


#ifndef FIELD_BASE_IMPL_HH_
#define FIELD_BASE_IMPL_HH_

#include <string>
#include <limits>
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


/******************************************************************************************
 * Implementation of FieldBase<...>
 */

template <int spacedim, class Value>
FieldBase<spacedim, Value>::FieldBase(unsigned int n_comp)
: time_( -numeric_limits<double>::infinity() ),
  value_(r_value_)
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
    FieldInterpolatedP0<spacedim,Value>::get_input_type(type, element_input_type);
    FieldElementwise<spacedim,Value>::get_input_type(type, element_input_type);

    return type;
}



template <int spacedim, class Value>
boost::shared_ptr< FieldBase<spacedim, Value> >
FieldBase<spacedim, Value>::function_factory(const Input::AbstractRecord &rec, unsigned int n_comp )
{
    boost::shared_ptr< FieldBase<spacedim, Value> > func;

    if (rec.type() == FieldInterpolatedP0<spacedim,Value>::input_type ) {
	//xprintf(PrgErr,"TYPE of Field currently not functional.\n");
	func=boost::make_shared< FieldInterpolatedP0<spacedim,Value> >(n_comp);
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
bool FieldBase<spacedim, Value>::set_time(double time) {
    time_ = time;
    return false; // no change
}



template <int spacedim, class Value>
void FieldBase<spacedim, Value>::set_mesh(const Mesh *mesh,  bool boundary_domain) {
}



template<int spacedim, class Value>
unsigned int FieldBase<spacedim, Value>::n_comp() const {
    return (Value::NRows_ ? 0 : value_.n_rows());
}






/******************************************************************************************
 * Implementation of Field<...>
 */

template<int spacedim, class Value>
Field<spacedim,Value>::Field()
: FieldCommonBase(false),
  read_field_descriptor_hook( &read_field_descriptor ),
  no_check_control_field_(NULL),
  is_fully_initialized_(false)
{
    this->enum_valued_ = boost::is_same<typename Value::element_type, FieldEnum>::value;
}


template<int spacedim, class Value>
it::AbstractRecord &  Field<spacedim,Value>::get_input_type() {
    return FieldBaseType::input_type;
}



/// ---------- Helper function template for make_input_tree method
template <class FieldBaseType>
IT::AbstractRecord get_input_type_resolution(Input::Type::Selection *sel,  const boost::true_type&)
{
    ASSERT( sel, "NULL pointer to selection in Field::get_input_type(), while Value==FieldEnum.\n");
    return FieldBaseType::get_input_type(sel);
}


template <class FieldBaseType>
IT::AbstractRecord get_input_type_resolution(Input::Type::Selection *sel,  const boost::false_type&)
{
    return FieldBaseType::get_input_type(NULL);
}
/// ---------- end helper function template




template<int spacedim, class Value>
it::AbstractRecord Field<spacedim,Value>::make_input_tree() {
    return get_input_type_resolution<FieldBaseType>( this->element_selection_ ,boost::is_same<typename Value::element_type, FieldEnum>());
}



template<int spacedim, class Value>
void Field<spacedim, Value>::disable_where(const Field<spacedim, typename FieldValue<spacedim>::Enum > *control_field, const vector<FieldEnum> &value_list) {
    no_check_control_field_=control_field;
    no_check_values_=value_list;
}


/*
template<int spacedim, class Value>
boost::shared_ptr< typename Field<spacedim,Value>::FieldBaseType >
Field<spacedim,Value>::operator[] (Region reg)
{
    ASSERT_LESS(reg.idx(), this->region_fields_.size());
    return this->region_fields_[reg.idx()];
}
*/


template <int spacedim, class Value>
bool Field<spacedim, Value>::get_const_accessor(Region reg, ElementAccessor<spacedim> &elm) {
	ASSERT_LESS(reg.idx(), this->region_fields_.size());
    boost::shared_ptr< FieldBaseType > region_field = this->region_fields_[reg.idx()];
    if (region_field && typeid(*region_field) == typeid(FieldConstant<spacedim, Value>)) {
        elm = ElementAccessor<spacedim>(mesh_, reg );
        return true;
    } else {
        return false;
    }
}

/*
template<int spacedim, class Value>
void Field<spacedim, Value>::set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) {
    boost::shared_ptr<FieldBaseType> field = FieldBaseType::function_factory(rec, this->n_comp_);
    set_field(domain, field);
}



template<int spacedim, class Value>
void Field<spacedim, Value>::set_field(const RegionSet &domain, boost::shared_ptr< FieldBaseType > field) {
    ASSERT( this->mesh_, "Null mesh pointer, set_mesh() has to be called before set_field().\n");
    if (domain.size() == 0) return;

    ASSERT_EQUAL( field->n_comp() , this->n_comp_);
    field->set_mesh( this->mesh_ , is_bc() );
    BOOST_FOREACH(Region reg, domain) region_fields_[reg.idx()] = field;
    changed_from_last_set_time_=true;
}
*/

template<int spacedim, class Value>
auto Field<spacedim, Value>::read_field_descriptor(Input::Record rec, const FieldCommonBase &field) -> SharedField
{
	Input::AbstractRecord field_record;
	if (rec.opt_val(field.name(), field_record))
		return FieldBaseType::function_factory(field_record, field.n_comp() );
	else
		return SharedField();
}




template<int spacedim, class Value>
bool Field<spacedim, Value>::set_time(const TimeGovernor &time, LimitSide side)
{
	static double last_set_time = -numeric_limits<double>::infinity();
	static LimitSide last_side;

    // We perform set_time only once for every time.
    if (time.t() == last_set_time && side == last_side) return false;
    last_set_time=time.t();
    last_side = side;

    // read all descriptors satifying time.ge(input_time)
    update_history(time);
    check_initialized_region_fields_();

    // set time on all regions
    // for regions that match type of the field domain
    BOOST_FOREACH(const Region &reg, this->mesh_->region_db().get_region_set("ALL") )
        if (reg.is_boundary() == this->bc_ && region_fields_[reg.idx()] ) {
        	double last_time_in_history = region_history_[reg.idx()].front().first;
        	unsigned int history_size=region_history_[reg.idx()].size();
        	unsigned int i_history;
        	if ( time.gt(last_time_in_history) ) {
        		// in smooth time
        		i_history=0;
        	} else {
        		// time .eq. input_time; i.e. jump time
        		if (side == LimitSide::right) {
        			i_history=0;
        		} else {
        			i_history=1;
        		}
        	}
        	i_history=min(i_history, history_size - 1);
        	ASSERT(i_history >= 0, "Emtpty field history.");
        	region_fields_[reg.idx()]=region_history_[reg.idx()].at(i_history).second;

//            bool changed = region_fields_[reg.idx()]->set_time(time);
//            this->changed_from_last_set_time_ = this->changed_from_last_set_time_ || changed;

        }

    this->changed_during_set_time = this->changed_from_last_set_time_;
    this->changed_from_last_set_time_ = false;
    return this->changed_during_set_time;
}


// helper functions
template<int spacedim, class FieldBaseType>
FieldEnum get_constant_enum_value_dispatch(boost::shared_ptr< FieldBaseType > region_field,  const boost::true_type&) {
    return region_field->value( typename Space<spacedim>::Point(), ElementAccessor<spacedim>());
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
void Field<spacedim,Value>::set_mesh(const Mesh &mesh) {
	this->mesh_ = &mesh;

	ASSERT_EQUAL(0, region_fields_.size());
	ASSERT_EQUAL(0, region_history_.size());

    // initialize table if it is empty, we assume that the RegionDB is closed at this moment
	region_fields_.resize( this->mesh_->region_db().size() );
    region_history_.resize( this->mesh_->region_db().size() );
}


template<int spacedim, class Value>
void Field<spacedim,Value>::update_history(const TimeGovernor &time) {
    ASSERT( this->mesh_, "Null mesh pointer, set_mesh() has to be called before.\n");

    // read input up to given time
	double input_time;
    if (input_list_.size() != 0) {
        while( list_it_ != input_list_.end()
        	   && time.ge( input_time = list_it_->val<double>("time") ) ) {

        	// get domain specification
        	RegionSet domain;
        	std::string name;
        	unsigned int id;
			if (list_it_->opt_val("r_set", name)) {
				domain = mesh_->region_db().get_region_set(name);

			} else if (list_it_->opt_val("region", name)) {
				domain.push_back( mesh_->region_db().find_label(name) );    // try find region by label
				if (! domain[0].is_valid() ) xprintf(Warn, "Unknown region with label: '%s'\n", name.c_str());

			} else if (list_it_->opt_val("rid", id)) {
				domain.push_back( mesh_->region_db().find_id(id) );         // try find region by ID
				if (! domain[0].is_valid() ) xprintf(Warn, "Unknown region with id: '%d'\n", id);

			} else {
//				THROW(ExcMissingDomain() << list_it_->ei_address() );
			}
		    if (domain.size() == 0) continue;

			// get field instance
			SharedField field_instance = read_field_descriptor_hook(*list_it_, *this);
			if (! field_instance) continue; // skip descriptors without related keys

			// add to history
		    ASSERT_EQUAL( field_instance->n_comp() , this->n_comp_);
		    field_instance->set_mesh( this->mesh_ , is_bc() );
		    BOOST_FOREACH(Region reg, domain) {
		    	region_history_[reg.idx()].push_front(
		    			HistoryPoint(input_time, field_instance)
		    	);
		    }

        	++list_it_;
        }
    }
}


template<int spacedim, class Value>
void Field<spacedim,Value>::check_initialized_region_fields_() {
	ASSERT(this->mesh_, "Null mesh pointer.");
    if (is_fully_initialized_) return;

    // check there are no empty field pointers, collect regions to be initialized from default value
    RegionSet regions_to_init; // empty vector

    BOOST_FOREACH(const Region &reg, this->mesh_->region_db().get_region_set("ALL") )
        if (reg.is_boundary() == this->bc_) {      		// for regions that match type of the field domain
            if (! region_history_[reg.idx()].size() )   // empty region history
            {
                if (no_check_control_field_) {      // is the check turned off?
                    FieldEnum value;
                    if (no_check_control_field_->get_constant_enum_value(reg, value)
                        && ( std::find(no_check_values_.begin(), no_check_values_.end(), value)
                             != no_check_values_.end() )
                       ) continue;                  // the field is not needed on this region
                }

                if (this->default_ != "") {    // try to use default
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
        reader.read_from_default(this->default_, a_rec_type );
        auto a_rec = reader.get_root_interface<Input::AbstractRecord>();
        auto field_ptr = FieldBaseType::function_factory(a_rec, this->n_comp_);
        field_ptr->set_mesh( this->mesh_ , is_bc() );
        BOOST_FOREACH(Region reg, regions_to_init) {
    		region_history_[reg.idx()]
    		                .push_front(HistoryPoint( 0.0, field_ptr) );
        }
    }
    is_fully_initialized_=true;
}



template<int spacedim, class Value>
BCField<spacedim, Value>::BCField() { this->bc_=true; }





/******************************************************************************************
 * Implementation of MultiField<...>
 */

template<int spacedim, class Value>
MultiField<spacedim, Value>::MultiField()
: FieldCommonBase(false)
{}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::init( const vector<string> &names) {
    sub_fields_.resize( names.size() );
    sub_names_ = names;
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
        sub_fields_[i_comp].set_name( this->name_ + "_" + sub_names_[i_comp] );
}



template<int spacedim, class Value>
it::AbstractRecord &  MultiField<spacedim,Value>::get_input_type() {
}



template<int spacedim, class Value>
it::AbstractRecord MultiField<spacedim,Value>::make_input_tree() {
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) {
}



template<int spacedim, class Value>
bool MultiField<spacedim, Value>::set_time(
		const TimeGovernor &time,
		LimitSide side)
{
	bool any=false;
	for( SubFieldType &field : sub_fields_) {
		any=any || field.set_time(time, side);
	}
    return any;
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_mesh(const Mesh &mesh) {
    this->mesh_ = &mesh;
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
        sub_fields_[i_comp].set_mesh(mesh);
}


#endif //FUNCTION_BASE_IMPL_HH_
