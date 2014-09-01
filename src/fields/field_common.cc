#include "fields/field_common.hh"

/****************************************************************************
 *  Implementation of FieldCommon
 */


FieldCommon::FieldCommon()
: shared_( std::make_shared<SharedData>() ),
  limit_side_(LimitSide::unknown),
  set_time_result_(TimeStatus::unknown)
{
    shared_->bc_=false;
    shared_->input_default_="";
    shared_->n_comp_ = 0;
    shared_->mesh_ = nullptr;
    shared_->is_fully_initialized_=false;
}



FieldCommon::FieldCommon(const FieldCommon & other)
: shared_(other.shared_),
  limit_side_(LimitSide::unknown),
  set_time_result_(TimeStatus::unknown)
{
     flags_.add( FieldFlag::input_copy );
}



IT::Record FieldCommon::field_descriptor_record(const string& record_name) {
    string rec_name, description;
    description = "Record to set fields of the equation.\n"
                "The fields are set only on the domain specified by one of the keys: 'region', 'rid', 'r_set'\n"
                "and after the time given by the key 'time'. The field setting can be overridden by\n"
                " any " + record_name + " record that comes later in the boundary data array.";

    IT::Record rec = IT::Record(record_name, description)
                     .declare_key("r_set", IT::String(), "Name of region set where to set fields.")
                     .declare_key("region", IT::String(), "Label of the region where to set fields. ")
                     .declare_key("rid", IT::Integer(0), "ID of the region where to set fields." )
                     .declare_key("time", IT::Double(0.0), IT::Default("0.0"),
                             "Apply field setting in this record after this time.\n"
                             "These times have to form an increasing sequence.");

    return rec;
}


void FieldCommon::set_input_list(const Input::Array &list)
{
    if (! flags().match(FieldFlag::declare_input)) return;

    shared_->input_list_ = list;

    // check that times forms ascending sequence
    double time,last_time=0.0;
    if (list.size() == 0) return;
    for( auto it = shared_->input_list_.begin<Input::Record>();
            it != shared_->input_list_.end(); ++it)
        if (it->find<Input::AbstractRecord>(input_name())) {
            // field descriptor appropriate to the field
            time = it->val<double>("time");
            if (time < last_time) {
                cout << shared_->input_list_.address_string();
                THROW( ExcNonascendingTime()
                        << EI_Time(time)
                        << EI_Field(input_name())
                        << shared_->input_list_.ei_address());
            }
            last_time=time;

        }

    shared_->list_it_ = shared_->input_list_.begin<Input::Record>();
}



void FieldCommon::mark_input_times(TimeMark::Type mark_type) {
    ASSERT_LESS( 0, shared_->input_list_.size());

    // pass through field descriptors containing key matching field name.
    double time,last_time=0.0;
    for( auto it = shared_->input_list_.begin<Input::Record>();
         it != shared_->input_list_.end(); ++it)
        if (it->find<Input::AbstractRecord>(input_name())) {
            time = it->val<double>("time"); // default time=0
            TimeGovernor::marks().add( TimeMark(time, mark_type | TimeGovernor::marks().type_input() ));
        }
}



FieldCommon::~FieldCommon() {}
