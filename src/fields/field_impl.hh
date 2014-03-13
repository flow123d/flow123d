/*
 * field_impl.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#ifndef FIELD_IMPL_HH_
#define FIELD_IMPL_HH_

#include <boost/foreach.hpp>

#include "field.hh"
#include "mesh/region.hh"
#include "input/json_to_storage.hh"


/******************************************************************************************
 * Implementation of Field<...>
 */

template<int spacedim, class Value>
Field<spacedim,Value>::Field()
: data_(std::make_shared<SharedData>()),
  read_field_descriptor_hook( &read_field_descriptor )

{
	// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
	// this invariant is kept also by n_comp setter
	shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
}


template<int spacedim, class Value>
Field<spacedim,Value>::Field(const string &name, bool bc)
: data_(std::make_shared<SharedData>()),
  read_field_descriptor_hook( &read_field_descriptor )
{
		// n_comp is nonzero only for variable size vectors Vector, VectorEnum, ..
		// this invariant is kept also by n_comp setter
		shared_->n_comp_ = (Value::NRows_ ? 0 : 1);
		shared_->bc_=bc;
		this->name( name );
}



template<int spacedim, class Value>
Field<spacedim,Value>::Field(const Field &other)
: FieldCommonBase(other),
  data_(other.data_),
  read_field_descriptor_hook( other.read_field_descriptor_hook )
{
	if (other.no_check_control_field_)
		no_check_control_field_ =  make_shared<ControlField>(*other.no_check_control_field_);

	// initialize region_fields_ vector
	// shared_is already same as other.shared_
	if (shared_->mesh_) this->set_mesh( *(shared_->mesh_) );

}


template<int spacedim, class Value>
Field<spacedim,Value> &Field<spacedim,Value>::operator=(const Field<spacedim,Value> &other)
{
	ASSERT(other.shared_->mesh_, "Must call set_mesh before assign to other field.\n");
	// check for self assignement
	if (&other == this) return *this;

	shared_ = other.shared_;
        shared_->is_fully_initialized_ = false;
	set_time_result_ = TimeStatus::unknown;
	limit_side_ = other.limit_side_;

	read_field_descriptor_hook = other.read_field_descriptor_hook;
	data_ = other.data_;

	if (other.no_check_control_field_) {
		no_check_control_field_ =  make_shared<ControlField>(*other.no_check_control_field_);
	}

	// initialize region_fields_ vector
	// shared_is already same as other.shared_
	this->set_mesh( *(shared_->mesh_) );


	return *this;
}



template<int spacedim, class Value>
it::AbstractRecord &Field<spacedim,Value>::get_input_type() {
	/*
	 * List of AbstratRecord types created by make_input_tree() in get_input_type() implementation.
	 * We have to return reference, which may be reference to not yet initialized static object.
	 *
	 * TODO: Have method to get persistent copy of an Input Type (which exists nevertheless)
	 */
	static vector<it::AbstractRecord> ar_list;

	if (is_enum_valued) {
		ar_list.push_back(make_input_tree());
		return ar_list.back();
	} else {
		return FieldBaseType::input_type;
	}
}



/// ---------- Helper function template for make_input_tree method
template <class FieldBaseType>
IT::AbstractRecord get_input_type_resolution(const Input::Type::Selection *sel,  const boost::true_type&)
{
    ASSERT( sel != nullptr,
    		"NULL pointer to selection in Field::get_input_type(), while Value==FieldEnum.\n");
    return FieldBaseType::get_input_type(sel);
}


template <class FieldBaseType>
IT::AbstractRecord get_input_type_resolution(const Input::Type::Selection *sel,  const boost::false_type&)
{
    return FieldBaseType::get_input_type(nullptr);
}
/// ---------- end helper function template




template<int spacedim, class Value>
it::AbstractRecord Field<spacedim,Value>::make_input_tree() {
	ASSERT(is_enum_valued,
			"Can not use make_input_tree() for non-enum valued fields, use get_inout_type() instead.\n" );
    return get_input_type_resolution<FieldBaseType>( shared_->element_selection_ ,boost::is_same<typename Value::element_type, FieldEnum>());
}



template<int spacedim, class Value>
auto Field<spacedim, Value>::disable_where(
		const Field<spacedim, typename FieldValue<spacedim>::Enum > &control_field,
		const vector<FieldEnum> &value_list) -> Field &
{
    no_check_control_field_=std::make_shared<ControlField>(control_field);
    shared_->no_check_values_=value_list;
    return *this;
}

template<int spacedim, class Value>
void Field<spacedim,Value>::set_mesh(const Mesh &in_mesh) {
	// since we allow copy of fields before set_mesh is called
	// we have to check that all copies set the same mesh
	if (shared_->mesh_ && shared_->mesh_ != &in_mesh) {
		THROW(ExcFieldMeshDifference() << EI_Field(name()) );
	}

	shared_->mesh_ = &in_mesh;

    // initialize table if it is empty, we assume that the RegionDB is closed at this moment
	region_fields_.resize( mesh()->region_db().size() );
	RegionHistory init_history(history_length_limit_);	// capacity
    data_->region_history_.resize( mesh()->region_db().size(), init_history );

    if (no_check_control_field_) no_check_control_field_->set_mesh(in_mesh);
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
bool Field<spacedim, Value>::is_constant(Region reg) {
	ASSERT_LESS(reg.idx(), this->region_fields_.size());
    FieldBasePtr region_field = this->region_fields_[reg.idx()];
    return (region_field && typeid(*region_field) == typeid(FieldConstant<spacedim, Value>));
}

/*
template<int spacedim, class Value>
void Field<spacedim, Value>::set_from_input(const RegionSet &domain, const Input::AbstractRecord &rec) {
    boost::shared_ptr<FieldBaseType> field = FieldBaseType::function_factory(rec, this->n_comp_);
    set_field(domain, field);
}

*/

template<int spacedim, class Value>
void Field<spacedim, Value>::set_field(
		const RegionSet &domain,
		FieldBasePtr field,
		double time)
{
    ASSERT( mesh(), "Null mesh pointer, set_mesh() has to be called before set_field().\n");
    if (domain.size() == 0) return;

    ASSERT_EQUAL( field->n_comp() , n_comp());
    field->set_mesh( mesh() , is_bc() );

    HistoryPoint hp = HistoryPoint(time, field);
    for(const Region &reg: domain) {
    	RegionHistory &region_history = data_->region_history_[reg.idx()];
    	// insert hp into descending time sequence
    	ASSERT( region_history.size() == 0 || region_history[0].first < hp.first, "Can not insert smaller time %g then last time %g in field's history.\n",
    			hp.first, region_history[0].first );
    	region_history.push_front(hp);
    }
    set_history_changed();
}



template<int spacedim, class Value>
void Field<spacedim, Value>::set_field(
		const RegionSet &domain,
		const Input::AbstractRecord &a_rec,
		double time)
{
	set_field(domain, FieldBaseType::function_factory(a_rec, n_comp()), time);
}



template<int spacedim, class Value>
auto Field<spacedim, Value>::read_field_descriptor(Input::Record rec, const FieldCommonBase &field) -> FieldBasePtr
{
	Input::AbstractRecord field_record;
	if (rec.opt_val(field.name(), field_record))
		return FieldBaseType::function_factory(field_record, field.n_comp() );
	else
		return FieldBasePtr();
}




template<int spacedim, class Value>
bool Field<spacedim, Value>::set_time(const TimeGovernor &time)
{
	ASSERT( mesh() , "NULL mesh pointer of field '%s'. set_mesh must be called before.\n",name().c_str());
	ASSERT( limit_side_ != LimitSide::unknown, "Must set limit side on field '%s' before calling set_time.\n",name().c_str());

    // We perform set_time only once for every time.
    if (time.t() == last_time_)  return changed();
    last_time_=time.t();

        // possibly update our control field
        if (no_check_control_field_) {
                no_check_control_field_->set_limit_side(limit_side_);
                no_check_control_field_->set_time(time);
        }
        
    set_time_result_ = TimeStatus::constant;
    
    // read all descriptors satisfying time.ge(input_time)
    update_history(time);
    check_initialized_region_fields_();

    // set time on all regions
    // for regions that match type of the field domain
    for(const Region &reg: mesh()->region_db().get_region_set("ALL") ) {
    	auto rh = data_->region_history_[reg.idx()];

    	// skip regions with no matching BC flag
    	// skipping regions with empty history - no_check regions
        if (reg.is_boundary() == is_bc() && !rh.empty() ) {
        	double last_time_in_history = rh.front().first;
        	unsigned int history_size=rh.size();
        	unsigned int i_history;
        	// set history index
        	if ( time.gt(last_time_in_history) ) {
        		// in smooth time
        		i_history=0;
        	} else {
        		// time .eq. input_time; i.e. jump time
        		if (limit_side_ == LimitSide::right) {
        			i_history=0;
        		} else {
        			i_history=1;
        		}
        	}
        	i_history=min(i_history, history_size - 1);
        	ASSERT(i_history >= 0, "Empty field history.");
        	// possibly update field pointer
        	auto new_ptr = rh.at(i_history).second;
        	if (new_ptr != region_fields_[reg.idx()]) {
        		region_fields_[reg.idx()]=new_ptr;
        		set_time_result_ = TimeStatus::changed;
        	}
        	// let FieldBase implementation set the time
    		//DBGMSG("Call particular set time, field: %s t: %g\n",this->name().c_str(), time.t());
    		if ( new_ptr->set_time(time.t()) )  set_time_result_ = TimeStatus::changed;

        }
    }

//    this->changed_during_set_time = this->changed_from_last_set_time_;
//    this->changed_from_last_set_time_ = false;
    return changed();
}


template<int spacedim, class Value>
void Field<spacedim, Value>::copy_from(const FieldCommonBase & other) {
	if (typeid(other) == typeid(*this)) {
		auto  const &other_field = dynamic_cast<  Field<spacedim, Value> const &>(other);
		this->operator=(other_field);
	}
}


template<int spacedim, class Value>
FieldResult Field<spacedim,Value>::field_result( ElementAccessor<spacedim> &elm) const {
    auto f = region_fields_[elm.region().idx()];
    if (f) return f->field_result();
    else return result_none;
}



template<int spacedim, class Value>
void Field<spacedim,Value>::update_history(const TimeGovernor &time) {
    ASSERT( mesh(), "Null mesh pointer, set_mesh() has to be called before.\n");

    // read input up to given time
	double input_time;
    if (shared_->input_list_.size() != 0) {
        while( shared_->list_it_ != shared_->input_list_.end()
        	   && time.ge( input_time = shared_->list_it_->val<double>("time") ) ) {

        	// get domain specification
        	RegionSet domain;
        	std::string domain_name;
        	unsigned int id;
			if (shared_->list_it_->opt_val("r_set", domain_name)) {
				domain = mesh()->region_db().get_region_set(domain_name);

			} else if (shared_->list_it_->opt_val("region", domain_name)) {
				domain.push_back( mesh()->region_db().find_label(domain_name) );    // try find region by label
				if (! domain[0].is_valid() ) xprintf(Warn, "Unknown region with label: '%s'\n", domain_name.c_str());

			} else if (shared_->list_it_->opt_val("rid", id)) {
				domain.push_back( mesh()->region_db().find_id(id) );         // try find region by ID
				if (! domain[0].is_valid() ) xprintf(Warn, "Unknown region with id: '%d'\n", id);

			} else {
				THROW(ExcMissingDomain()
						<< EI_Field(this->name())
						<< shared_->list_it_->ei_address() );
			}
		    if (domain.size() == 0) continue;

			// get field instance
			FieldBasePtr field_instance = read_field_descriptor_hook(*(shared_->list_it_), *this);
			if (field_instance)  // skip descriptors without related keys
			{
				// add to history
				ASSERT_EQUAL( field_instance->n_comp() , n_comp());
				field_instance->set_mesh( mesh() , is_bc() );
				for(const Region &reg: domain) {
					data_->region_history_[reg.idx()].push_front(
							HistoryPoint(input_time, field_instance)
					);
				}
			}
        	++shared_->list_it_;
        }
    }
}

template<int spacedim, class Value>
void Field<spacedim,Value>::check_initialized_region_fields_() {
	ASSERT(mesh(), "Null mesh pointer.");
    if (shared_->is_fully_initialized_) return;

    // check there are no empty field pointers, collect regions to be initialized from default value
    RegionSet regions_to_init; // empty vector

    for(const Region &reg : mesh()->region_db().get_region_set("ALL") )
        if (reg.is_boundary() == is_bc()) {      		// for regions that match type of the field domain
            RegionHistory &rh = data_->region_history_[reg.idx()];
        	if ( rh.empty() ||	! rh[0].second)   // empty region history
            {
        		// test if check is turned on and control field is FieldConst
                if (no_check_control_field_ && no_check_control_field_->is_constant(reg) ) {
                	// get constant enum value
                	auto elm = ElementAccessor<spacedim>(mesh(), reg);
                	FieldEnum value = no_check_control_field_->value(elm.centre(),elm);
                	// check that the value is in the disable list
                    if ( std::find(shared_->no_check_values_.begin(), shared_->no_check_values_.end(), value)
                             != shared_->no_check_values_.end() )
                        continue;                  // the field is not needed on this region
                }
                if (shared_->default_ != "") {    // try to use default
                    regions_to_init.push_back( reg );
                } else {
                	xprintf(UsrErr, "Missing value of the field '%s' on region ID: %d label: %s.\n",
                			name().c_str(), reg.id(), reg.label().c_str() );
                }
            }
        }

    // possibly set from default value
    if ( regions_to_init.size() ) {
    	xprintf(Warn, "Using default value '%s' for part of field '%s'.\n", input_default().c_str(), name().c_str());

    	// has to deal with fact that reader can not deal with input consisting of simple values
    	string default_input=input_default();
    	auto input_type = get_input_type();
        Input::JSONToStorage reader( default_input, input_type );

        auto a_rec = reader.get_root_interface<Input::AbstractRecord>();
        auto field_ptr = FieldBaseType::function_factory( a_rec , n_comp() );
        field_ptr->set_mesh( mesh(), is_bc() );
        for(const Region &reg: regions_to_init) {
    		data_->region_history_[reg.idx()]
    		                .push_front(HistoryPoint( 0.0, field_ptr) );
        }
    }
    shared_->is_fully_initialized_;
}



//template<int spacedim, class Value>
//BCField<spacedim, Value>::BCField() { this->bc_=true; }





/******************************************************************************************
 * Implementation of MultiField<...>
 */

template<int spacedim, class Value>
MultiField<spacedim, Value>::MultiField()
: FieldCommonBase()
{}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::init( const vector<string> &names) {
    sub_fields_.resize( names.size() );
    sub_names_ = names;
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
        sub_fields_[i_comp].name( name() + "_" + sub_names_[i_comp] );
}



template<int spacedim, class Value>
it::AbstractRecord &  MultiField<spacedim,Value>::get_input_type() {
}


template<int spacedim, class Value>
bool MultiField<spacedim, Value>::set_time(
		const TimeGovernor &time)
{
	bool any=false;
	for( SubFieldType &field : sub_fields_) {
		any=any || field.set_time(time);
	}
    return any;
}



template<int spacedim, class Value>
void MultiField<spacedim, Value>::set_mesh(const Mesh &mesh) {
    shared_->mesh_ = &mesh;
    for(unsigned int i_comp=0; i_comp < size(); i_comp++)
        sub_fields_[i_comp].set_mesh(mesh);
}


template<int spacedim, class Value>
void MultiField<spacedim, Value>::copy_from(const FieldCommonBase & other) {
	if (typeid(other) == typeid(*this)) {
		auto  const &other_field = dynamic_cast<  MultiField<spacedim, Value> const &>(other);
		this->operator=(other_field);
	} else if (typeid(other) == typeid(SubFieldType)) {
		auto  const &other_field = dynamic_cast<  SubFieldType const &>(other);
		sub_fields_.resize(1);
		sub_fields_[0] = other_field;
	}
}



template<int spacedim, class Value>
bool MultiField<spacedim, Value>::is_constant(Region reg) {
	bool const_all=false;
	for(auto field : sub_fields_) const_all = const_all || field.is_constant(reg);
	return const_all;
}





#endif /* FIELD_IMPL_HH_ */
