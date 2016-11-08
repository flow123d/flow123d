/*
 * equation_output.cc
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#include "tools/time_marks.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "io/equation_output.hh"
#include "io/output_time_set.hh"
#include "input/flow_attribute_lib.hh"
#include <memory>


namespace IT = Input::Type;



IT::Record &EquationOutput::get_input_type() {

    static const IT::Record &field_output_setting =
        IT::Record("FieldOutputSetting", "Setting of the field output. The field name, output times, output interpolation (future).")
            .allow_auto_conversion("field")
            .declare_key("field", IT::Parameter("output_field_selection"), IT::Default::obligatory(),
                    "The field name (from selection).")
            .declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
                    "Output times specific to particular field.")
            //.declare_key("interpolation", ...)
            .close();

    return IT::Record("EquationOutput",
            "Output of the equation's fields."
            "The output is done through the output stream of the associated balance law equation."
            "The stream defines output format for the full space information in selected times and "
            "observe points for the full time information. The key 'fields' select the fields for the full spatial output."
            "The set of output times may be specified  per field otherwise common time set 'times' is used. If even this is not provided"
            "the time set of the output_stream is used. The initial time of the equation is automatically added "
            "to the time set of every selected field. The end time of the equation is automatically added "
            "to the common output time set.")
        .root_of_generic_subtree()
        .declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
                "Output times used for the output fields without is own time series specification.")
        .declare_key("add_input_times", IT::Bool(), IT::Default("false"),
                "Add all input time points of the equation, mentioned in the 'input_fields' list, also as the output points.")
        .declare_key("fields", IT::Array(field_output_setting), IT::Default("[]"),
                "Array of output fields and their individual output settings.")
        .declare_key("observe_fields", IT::Array( IT::Parameter("output_field_selection")), IT::Default("[]"),
                "Array of the fields evaluated in the observe points of the associated output stream.")
        .close();
}



const IT::Instance &EquationOutput::make_output_type(const string &equation_name, const string &additional_description)
{
    string selection_name = equation_name + "_output_fields";
    string description = "Selection of output fields for the " + equation_name + " model.\n" + additional_description;
    IT::Selection sel(selection_name, description );
    int i=0;
    // add value for each field excluding boundary fields
    for( FieldCommon * field : field_list)
    {
        //DebugOut().fmt("type for field: {}\n", field->name());
        if ( !field->is_bc() && field->flags().match( FieldFlag::allow_output) )
        {
            string desc = "Output of the field " + field->name() + " (($[" + field->units().format_latex()+"]$))";
            if (field->description().length() > 0)
                desc += " (" + field->description() + ").";
            else
                desc += ".";
            sel.add_value(i, field->name(), desc, { {FlowAttribute::field_value_shape(), field->get_value_attribute()} });
            i++;
        }
    }

    const IT::Selection &output_field_selection = sel.close();

    std::vector<IT::TypeBase::ParameterPair> param_vec;
    param_vec.push_back( std::make_pair("output_field_selection", std::make_shared< IT::Selection >(output_field_selection) ) );
    return IT::Instance(get_input_type(), param_vec).close();

}


void EquationOutput::initialize(std::shared_ptr<OutputTime> stream, Input::Record in_rec, const TimeGovernor & tg)
{
    stream_ = stream;
    equation_type_ = tg.equation_mark_type();
    equation_fixed_type_ = tg.equation_fixed_mark_type();
    read_from_input(in_rec, tg);
}



void EquationOutput::read_from_input(Input::Record in_rec, const TimeGovernor & tg)
{
    ASSERT(stream_).error("The 'set_stream' method must be called before the 'read_from_input'.");
    auto &marks = TimeGovernor::marks();

    Input::Array times_array;
    if (in_rec.opt_val("times", times_array) ) {
        common_output_times_.read_from_input(times_array, tg);
    } else {
        // take times from the output_stream if key times is missing
        auto times_array_it = stream_->get_time_set_array();
        if (times_array_it) {
            common_output_times_.read_from_input(*times_array_it,  tg);
        }
    }
    // always add the end time
    common_output_times_.add(tg.end_time(), equation_fixed_type_);

    if (in_rec.val<bool>("add_input_times")) {
        // copy time marks in order to prevent invalidation of the iterator
        TimeMarks marks_copy = TimeGovernor::marks();
        for(auto time_mark_it = marks_copy.begin(equation_type_ | marks.type_input());
                time_mark_it != marks_copy.end(equation_type_ | marks.type_input());
                ++time_mark_it) {
            common_output_times_.add(time_mark_it->time(), equation_fixed_type_);
        }
    }
    auto fields_array = in_rec.val<Input::Array>("fields");
    for(auto it = fields_array.begin<Input::Record>(); it != fields_array.end(); ++it) {
        string field_name = it -> val< Input::FullEnum >("field");
        Input::Array field_times_array;
        if (it->opt_val("times", field_times_array)) {
            OutputTimeSet field_times;
            field_times.read_from_input(field_times_array, tg);
            field_output_times_[field_name] = field_times;
        } else {
            field_output_times_[field_name] = common_output_times_;
        }
        // Add init time as the output time for every output field.
        field_output_times_[field_name].add(tg.init_time(), equation_fixed_type_);
    }
    auto observe_fields_array = in_rec.val<Input::Array>("observe_fields");
    for(auto it = observe_fields_array.begin<Input::FullEnum>(); it != observe_fields_array.end(); ++it) {
        observe_fields_.insert(string(*it));
    }
}

bool EquationOutput::is_field_output_time(const FieldCommon &field, TimeStep step) const
{
    auto &marks = TimeGovernor::marks();
    auto field_times_it = field_output_times_.find(field.name());
    if (field_times_it == field_output_times_.end()) return false;
    ASSERT( step.eq(field.time()) )(step.end())(field.time())(field.name()).error("Field is not set to the output time.");
    auto current_mark_it = marks.current(step, equation_type_ | marks.type_output() );
    if (current_mark_it == marks.end(equation_type_ | marks.type_output()) ) return false;
    return (field_times_it->second.contains(*current_mark_it) );
}

void EquationOutput::output(TimeStep step)
{
    // TODO: remove const_cast after resolving problems with const Mesh.
    //Mesh *field_mesh = const_cast<Mesh *>(field_list[0]->mesh());
    stream_->make_output_mesh(*this);

    for(FieldCommon * field : this->field_list) {

        if ( field->flags().match( FieldFlag::allow_output) ) {
            if (is_field_output_time(*field, step)) {
                field->output(stream_);
            }
            // observe output
            if (observe_fields_.find(field->name()) != observe_fields_.end()) {
                field->observe_output( stream_->observe() );
            }
        }
    }
}


void EquationOutput::add_output_times(double begin, double step, double end)
{
    common_output_times_.add(begin,step, end, equation_fixed_type_ );
}
