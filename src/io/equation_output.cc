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


namespace IT = Input::Type;

TimeMark OutputTimeSet::from_time_mark_storage(TimeMarkStorage tm_storage)
{
    return TimeMark(tm_storage.first, tm_storage.second);
}
auto OutputTimeSet::to_time_mark_storage(TimeMark mark) -> TimeMarkStorage
{
    return TimeMarkStorage(mark.time(), mark.mark_type());
}



const IT::Array OutputTimeSet::get_input_type()
{
    static const IT::Record &time_grid =
        IT::Record("TimeGrid", "Equally spaced grid of time points.")
            .allow_auto_conversion("begin")
            .declare_key("begin", IT::Double(0.0), IT::Default::obligatory(),
                    "The start time of the grid.")
            .declare_key("step", IT::Double(0.0), IT::Default::optional(),
                    "The step of the grid. If not specified, the grid consists only of the start time.")
            .declare_key("end", IT::Double(0.0), IT::Default::read_time("The end time of the simulation."),
                    "The time greater or equal to the last time in the grid.")
            .close();
    return IT::Array(time_grid);
}


void OutputTimeSet::read_from_input(Input::Array in_array, TimeMark::Type mark_type)
{
    auto marks = TimeGovernor::marks();

    for(auto it =in_array.begin<Input::Record>(); it != in_array.end(); ++it) {
        double t_begin = it->val<double>("begin");
        double simulation_end_time = marks.end(TimeMark::every_type)->time();
        double t_end = it->val<double>("end", simulation_end_time );
        double t_step;
        if (! it->opt_val("step", t_step) ) {
            t_end = t_begin;
            t_step = 1.0;
        }
        if ( t_begin > t_end) {
            WarningOut().fmt("Ignoring output time grid. Time begin {}  > time end {}. {}",
                    t_begin, t_end, it->address_string());
            continue;
        }
        if ( t_step < 2*numeric_limits<double>::epsilon()) {
            WarningOut().fmt("Ignoring output time grid. Time step {} < two times machine epsilon {}. {}",
                    t_step, 2*numeric_limits<double>::epsilon(),it->address_string());
            continue;
        }

        TimeMark::Type output_mark_type = mark_type | marks.type_output();

        unsigned int n_steps=((t_end - t_begin) / t_step + TimeGovernor::time_step_precision);
        for(unsigned int i = 0; i <= n_steps; i++) {
            auto mark = TimeMark(t_begin + i * t_step, output_mark_type);
            time_marks_.insert( to_time_mark_storage( marks.add(mark) ) );
        }
    }
}


bool OutputTimeSet::contains(TimeMark mark) {
    return time_marks_.find(to_time_mark_storage(mark)) != time_marks_.end();
}


const IT::Record &EquationOutput::get_input_type() {

    static const IT::Record &field_output_setting =
        IT::Record("FieldOutputSetting", "Setting of the field output. The field name, output times, output interpolation (future).")
            .allow_auto_conversion("field")
            .declare_key("field", IT::Parameter("output_field_selection"), IT::Default::obligatory(),
                    "The field name (from selection).")
            .declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
                    "Output times specific to particular field.")
            //.declare_key("interpolation", ...)
            .close();

    return IT::Record("EquationOutput", "Configuration of fields output. "
            "The output is done through the output stream of the associated balance law equation.")
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


EquationOutput::EquationOutput(std::shared_ptr<OutputTime> stream, TimeMark::Type mark_type)
{
    stream_ = stream;
    equation_type_ = mark_type;
}



void EquationOutput::read_from_input(Input::Record in_rec)
{
    auto marks = TimeGovernor::marks();

    Input::Array times_array;
    if (in_rec.opt_val("times", times_array) ) {
        common_output_times_.read_from_input(times_array, equation_type_ );
    }
    if (in_rec.val<bool>("add_input_times")) {
        marks.add_to_type_all( equation_type_ | marks.type_input(), equation_type_ | marks.type_output() );
    }
    auto fields_array = in_rec.val<Input::Array>("fields");
    for(auto it = fields_array.begin<Input::Record>(); it != fields_array.end(); ++it) {
        string field_name = it -> val< Input::FullEnum >("field");
        Input::Array field_times_array;
        if (it->opt_val("times", field_times_array)) {
            OutputTimeSet field_times;
            field_times.read_from_input(field_times_array, equation_type_);
            field_output_times_[field_name] = field_times;
        } else {
            field_output_times_[field_name] = common_output_times_;
        }
    }
    auto observe_fields_array = in_rec.val<Input::Array>("observe_fields");
    for(auto it = fields_array.begin<Input::FullEnum>(); it != fields_array.end(); ++it) {
        observe_fields_.insert(string(*it));
    }
}


void EquationOutput::output(TimeStep step)
{

    auto &marks = TimeGovernor::marks();
    for(FieldCommon * field : this->field_list) {
        // stream output
        auto field_times_it = field_output_times_.find(field->name());
        if (field_times_it != field_output_times_.end()) {
            ASSERT( step.eq(field->time()) )(step.end())(field->time())(field->name()).error("Field is not set to the output time.");
            auto current_mark_it = marks.current(step, equation_type_ | marks.type_output() );
            if (current_mark_it != marks.end(equation_type_ | marks.type_output()) ) {
                if (field_times_it->second.contains(*current_mark_it) ) {
                    field->output(stream_);
                }
            }
        }
        // observe output
        if (observe_fields_.find(field->name()) != observe_fields_.end()) {
            field->observe_output(stream_->observe());
        }
    }
}
