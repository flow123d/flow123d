/*
 * field_base.cc
 *
 *  Created on: Dec 4, 2012
 *      Author: jb
 */



#include "fields/field_base_impl.hh"
#include "coupling/time_governor.hh"
#include "system/exceptions.hh"
#include <limits>


/****************************************************************************
 *  Implementation of FieldCommon
 */

FieldCommonBase::FieldCommonBase(bool bc)
: changed_flag_({false, false}),       // reading this variable is in fact invalid up to the first call of the set_time
  bc_(bc),
  n_comp_(0),
  element_selection_(NULL),
  default_( "" ),
  mesh_(NULL)
{

}



void FieldCommonBase::set_input_list(const Input::Array &list)
{
	input_list_ = list;

	// check that times forms ascending sequence
	double time,last_time=0.0;
	if (list.size() == 0) return;
	for( list_it_ = input_list_.begin<Input::Record>();
	     list_it_ != input_list_.end(); ++list_it_)
		if (list_it_->find<Input::AbstractRecord>(name())) {
			// field descriptor appropriate to the field
		    time = list_it_->val<double>("time");
		    if (time < last_time) {
//		    	THROW( ExcNonascendingTime()
//		    			<< EI_Time(time)
//		    			<< EI_Field(name)
//		    			<< input_list_.ei_address());
		    }
		    last_time=time;

		}

	list_it_ = input_list_.begin<Input::Record>();
}



// getters
const std::string & FieldCommonBase::name() const
{ return name_; }


const std::string  FieldCommonBase::desc() const
{ if(default_ != "")
    return "Default Field value: " + default_ + " \n " + desc_;
  else
    return desc_;
}


const string & FieldCommonBase::get_default() const
{ return default_; }


bool FieldCommonBase::is_bc() const
{ return bc_; }


bool FieldCommonBase::is_enum_valued() const
{ return enum_valued_; }


unsigned int FieldCommonBase::n_comp() const
{ return n_comp_; }


const Mesh * FieldCommonBase::mesh() const
{ return mesh_; }



void FieldCommonBase::mark_input_times(TimeMark::Type mark_type) {
	ASSERT_LESS( 0, input_list_.size());

	// pass through field descriptors containing key matching field name.
	double time,last_time=0.0;
	for( auto it = input_list_.begin<Input::Record>();
	     it != input_list_.end(); ++it)
		if (it->find<Input::AbstractRecord>(name())) {
		    time = it->val<double>("time"); // default time=0
		    TimeGovernor::marks().add( TimeMark(time, mark_type | TimeGovernor::marks().type_input() ));
		}
}



FieldCommonBase::~FieldCommonBase() {}






/****************************************************************************
 *  Implementation of MultiField
 */




/****************************************************************************
 *  Instances of field templates
 */



INSTANCE_ALL(FieldBase)
INSTANCE_ALL(Field)
INSTANCE_ALL(BCField)
INSTANCE_ALL(MultiField)

