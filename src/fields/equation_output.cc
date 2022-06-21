/*
 * equation_output.cc
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#include "tools/time_marks.hh"
#include "input/input_type.hh"
#include "input/accessors.hh"
#include "fields/equation_output.hh"
#include "fields/field.hh"
#include "fields/assembly_output.hh"
#include "fields/assembly_observe.hh"
#include "io/output_time_set.hh"
#include "io/observe.hh"
#include "input/flow_attribute_lib.hh"
#include "fem/dofhandler.hh"
#include "fem/discrete_space.hh"
#include "fem/fe_p.hh"
#include <memory>


namespace IT = Input::Type;



IT::Record &EquationOutput::get_input_type() {

    static const IT::Selection &interpolation_sel =
        IT::Selection("Discrete_output", "Discrete type of output. Determines type of output data (element, node, native etc).")
            .add_value(OutputTime::NODE_DATA,   "P1_average", "Node data / point data.")
			.add_value(OutputTime::CORNER_DATA, "D1_value",   "Corner data.")
			.add_value(OutputTime::ELEM_DATA,   "P0_value",   "Element data / cell data.")
			.add_value(OutputTime::NATIVE_DATA, "Native",     "Native data (Flow123D data).")
			.close();

    static const IT::Record &field_output_setting =
        IT::Record("FieldOutputSetting", "Setting of the field output. The field name, output times, output interpolation (future).")
            .allow_auto_conversion("field")
            .declare_key("field", IT::Parameter("output_field_selection"), IT::Default::obligatory(),
                    "The field name (from selection).")
            .declare_key("times", OutputTimeSet::get_input_type(), IT::Default::optional(),
                    "Output times specific to particular field.")
            .declare_key("interpolation", IT::Array( interpolation_sel ), IT::Default::read_time("Interpolation type of output data."),
					"Optional value. Implicit value is given by field and can be changed.")
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



EquationOutput::EquationOutput()
: FieldSet(), output_elem_data_assembly_(nullptr), output_node_data_assembly_(nullptr), output_corner_data_assembly_(nullptr),
  observe_output_assembly_(nullptr) {
    this->add_coords_field();
}

EquationOutput::~EquationOutput() {
    if (output_elem_data_assembly_ != nullptr) delete output_elem_data_assembly_;
    if (output_node_data_assembly_ != nullptr) delete output_node_data_assembly_;
    if (output_corner_data_assembly_ != nullptr) delete output_corner_data_assembly_;
    if (observe_output_assembly_ != nullptr) delete observe_output_assembly_;
}



const IT::Selection &EquationOutput::create_output_field_selection(const string &equation_name,
                                                                   const string &additional_description)
{
    string selection_name = equation_name + ":OutputFields";
    string description = "Selection of output fields for the " + equation_name + " model.\n" + additional_description;
    IT::Selection sel(selection_name, description );
    int i=0;
    // add value for each field excluding boundary fields
    for( FieldCommon * field : field_list)
    {
        //DebugOut().fmt("type for field: {}\n", field->name());
        if ( !field->is_bc() && field->flags().match( FieldFlag::allow_output) )
        {
            string desc = "(($[" + field->units().format_latex()+"]$)) "; + "Output of: the field " + field->name() + " ";
            if (field->flags().match(FieldFlag::equation_input))
                desc += "Input field: ";
            if (field->description().length() > 0)
                desc += field->description();
            sel.add_value(i, field->name(), desc, { {FlowAttribute::field_value_shape(), field->get_value_attribute()} });
            i++;
        }
    }

    return sel.close();
}

const IT::Instance &EquationOutput::make_output_type(const string &equation_name, const string &additional_description)
{
    return make_output_type_from_record(get_input_type(), equation_name, additional_description);
}

const IT::Instance &EquationOutput::make_output_type_from_record(Input::Type::Record &in_rec,
                                                                 const string &equation_name,
                                                                 const string &additional_description)
{
    const IT::Selection &output_field_selection = create_output_field_selection(equation_name, additional_description);

    std::vector<IT::TypeBase::ParameterPair> param_vec;
    param_vec.push_back( std::make_pair("output_field_selection", std::make_shared< IT::Selection >(output_field_selection) ) );
    return IT::Instance(in_rec, param_vec).close();
}


void EquationOutput::initialize(std::shared_ptr<OutputTime> stream, Mesh *mesh, Input::Record in_rec, const TimeGovernor & tg)
{
    stream_ = stream;
    mesh_ = mesh;
    equation_type_ = tg.equation_mark_type();
    equation_fixed_type_ = tg.equation_fixed_mark_type();
    read_from_input(in_rec, tg);

    { // DOF handler of element data output
        MixedPtr<FE_P_disc> fe_p_disc(0);
        dh_ = make_shared<DOFHandlerMultiDim>(*mesh_);
	    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_p_disc);
	    dh_->distribute_dofs(ds);
    }

    { // DOF handler of node / corner data output
        MixedPtr<FE_P_disc> fe_p_disc(1);
        dh_node_ = make_shared<DOFHandlerMultiDim>(*mesh_);
	    std::shared_ptr<DiscreteSpace> ds = std::make_shared<EqualOrderDiscreteSpace>( mesh_, fe_p_disc);
	    dh_node_->distribute_dofs(ds);
    }

    output_elem_data_assembly_ = new GenericAssembly< AssemblyOutputElemData >(this, this);
    output_node_data_assembly_ = new GenericAssembly< AssemblyOutputNodeData >(this, this);
    output_corner_data_assembly_ = new GenericAssembly< AssemblyOutputNodeData >(this, this);
    observe_output_assembly_ = new GenericAssemblyObserve< AssemblyObserveOutput >(this, this->observe_fields_, stream_->observe( mesh_ ));
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
        this->init_field_item(it, tg);
    }
    auto observe_fields_array = in_rec.val<Input::Array>("observe_fields");
    for(auto it = observe_fields_array.begin<Input::FullEnum>(); it != observe_fields_array.end(); ++it) {
        observe_fields_.insert(string(*it));
    }

}

void EquationOutput::init_field_item(Input::Iterator<Input::Record> it, const TimeGovernor & tg) {
    string field_name = it -> val< Input::FullEnum >("field");
    FieldCommon *found_field = field(field_name);

    Input::Array interpolations;
    OutputTime::DiscreteSpaceFlags interpolation = OutputTime::empty_discrete_flags();
    if (it->opt_val("interpolation", interpolations)) {
        // process interpolations
        for(auto it_interp = interpolations.begin<OutputTime::DiscreteSpace>(); it_interp != interpolations.end(); ++it_interp) {
            interpolation[ *it_interp ] = true;
        }
    } else {
        OutputTime::set_discrete_flag(interpolation, found_field->get_output_type());
    }
    Input::Array field_times_array;
    FieldOutputConfig field_config;
    if (it->opt_val("times", field_times_array)) {
        OutputTimeSet field_times;
        field_times.read_from_input(field_times_array, tg);
        field_config.output_set_ = field_times;
    } else {
        field_config.output_set_ = common_output_times_;
    }
    field_config.space_flags_ = interpolation;
    // Add init time as the output time for every output field.
    field_config.output_set_.add(tg.init_time(), equation_fixed_type_);
    // register interpolation types of fields to OutputStream
    for (uint i=0; i<OutputTime::N_DISCRETE_SPACES; ++i)
	    if (interpolation[i]) used_interpolations_.insert( OutputTime::DiscreteSpace(i) );
    // Set output configuration to field_output_times_
    if (found_field->is_multifield()) {
        for (uint i_comp=0; i_comp<found_field->n_comp(); ++i_comp) {
            field_output_times_[ found_field->full_comp_name(i_comp) ] = field_config;
        }
    } else {
        field_output_times_[field_name] = field_config;
    }
}

bool EquationOutput::is_field_output_time(const FieldCommon &field, TimeStep step) const
{
    if ( !field.get_flags().match(FieldFlag::allow_output) ) return false;
    auto &marks = TimeGovernor::marks();
    auto field_times_it = field_output_times_.find(field.name());
    if (field_times_it == field_output_times_.end()) return false;
    ASSERT( step.eq(field.time()) )(step.end())(field.time())(field.name()).error("Field is not set to the output time.");
    auto current_mark_it = marks.current(step, equation_type_ | marks.type_output() );
    if (current_mark_it == marks.end(equation_type_ | marks.type_output()) ) return false;
    return (field_times_it->second.output_set_.contains(*current_mark_it) );
}


void EquationOutput::output(TimeStep step)
{
    ASSERT_PTR(mesh_).error();

    // automatically call of stream_->write_time_frame if the time in the TimeStep is higher then in output stream
    if (step.end() > stream_->registered_time()) {
        stream_->write_time_frame();
    }

    // make observe points if not already done
	auto observe_ptr = stream_->observe(mesh_);

    this->make_output_mesh( stream_->is_parallel() );

    // NODE_DATA
    {
        FieldSet used_fields;
        for(FieldListAccessor f_acc : this->fields_range()) {
            if (is_field_output_time( *(f_acc.field()), step) && field_output_times_[f_acc->name()].space_flags_[OutputTime::NODE_DATA]) {
                f_acc->set_output_data_cache(OutputTime::NODE_DATA, stream_);
                used_fields += *(f_acc.field());
            }
        }
        if (used_fields.size()>0) {
            auto mixed_assmbly = output_node_data_assembly_->multidim_assembly();
            mixed_assmbly[1_d]->set_output_data(used_fields, stream_);
            mixed_assmbly[2_d]->set_output_data(used_fields, stream_);
            mixed_assmbly[3_d]->set_output_data(used_fields, stream_);
            output_node_data_assembly_->assemble(this->dh_node_);
        }
    }

    // CORNER_DATA
    {
        FieldSet used_fields;
        for(FieldListAccessor f_acc : this->fields_range()) {
            if (is_field_output_time( *(f_acc.field()), step) && field_output_times_[f_acc->name()].space_flags_[OutputTime::CORNER_DATA]) {
                f_acc->set_output_data_cache(OutputTime::CORNER_DATA, stream_);
                used_fields += *(f_acc.field());
            }
        }
        if (used_fields.size()>0) {
            auto mixed_assmbly = output_corner_data_assembly_->multidim_assembly();
            mixed_assmbly[1_d]->set_output_data(used_fields, stream_);
            mixed_assmbly[2_d]->set_output_data(used_fields, stream_);
            mixed_assmbly[3_d]->set_output_data(used_fields, stream_);
            output_corner_data_assembly_->assemble(this->dh_node_);
        }
    }

    // ELEM_DATA
    {
        FieldSet used_fields;
        for (FieldListAccessor f_acc : this->fields_range()) {
            if (is_field_output_time( *(f_acc.field()), step) && field_output_times_[f_acc->name()].space_flags_[OutputTime::ELEM_DATA]) {
                f_acc->set_output_data_cache(OutputTime::ELEM_DATA, stream_);
                used_fields += *(f_acc.field());
            }
        }
        if (used_fields.size()>0) {
            auto mixed_assmbly = output_elem_data_assembly_->multidim_assembly();
            mixed_assmbly[1_d]->set_output_data(used_fields, stream_);
            mixed_assmbly[2_d]->set_output_data(used_fields, stream_);
            mixed_assmbly[3_d]->set_output_data(used_fields, stream_);
            output_elem_data_assembly_->assemble(this->dh_);
        }
    }

    // NATIVE_DATA
    for(FieldListAccessor f_acc : this->fields_range()) {
        if (is_field_output_time( *(f_acc.field()), step) && field_output_times_[f_acc->name()].space_flags_[OutputTime::NATIVE_DATA]) {
            f_acc->field_output(stream_, OutputTime::NATIVE_DATA);
        }
    }

    // observe output
    if (observe_fields_.size() > 0) {
        for (auto observe_field : this->observe_fields_) {
            auto *field_ptr = this->field(observe_field);
            if ( field_ptr->flags().match( FieldFlag::allow_output) ) {
                if (field_ptr->is_multifield()) {
                    for (uint i_comp=0; i_comp<field_ptr->n_comp(); ++i_comp) {
                        observe_ptr->prepare_compute_data(field_ptr->full_comp_name(i_comp), step.end(), field_ptr->n_shape());
                    }
                } else {
                    observe_ptr->prepare_compute_data(field_ptr->name(), field_ptr->time(), field_ptr->n_shape());
                }
            }
        }
        observe_output_assembly_->assemble(this->dh_);
    }
}


void EquationOutput::add_output_times(double begin, double step, double end)
{
    common_output_times_.add(begin,step, end, equation_fixed_type_ );
}


void EquationOutput::make_output_mesh(bool parallel)
{
    // already computed
    if (stream_->is_output_data_caches_init()) return;

    // Read optional error control field name
    bool need_refinment = stream_->get_output_mesh_record();

    if(need_refinment) {
        if(stream_->enable_refinement()) {
            // create output meshes from input record
        	output_mesh_ = std::make_shared<OutputMeshDiscontinuous>(*mesh_, *stream_->get_output_mesh_record());

            // possibly set error control field for refinement
            auto ecf = select_error_control_field();
            output_mesh_->set_error_control_field(ecf);

            // actually compute refined mesh
            output_mesh_->create_refined_sub_mesh();
            output_mesh_->make_serial_master_mesh();

            stream_->set_output_data_caches(output_mesh_);
            return;
        }
        else
        {
            // skip creation of output mesh (use computational one)
           	WarningOut() << "Ignoring output mesh record.\n Output in GMSH format available only on computational mesh!";
        }
    }

    // create output mesh identical with the computational one
	bool discont = need_refinment | (used_interpolations_.find(OutputTime::CORNER_DATA) != used_interpolations_.end());
	//discont |= parallel;
	if (discont) {
		output_mesh_ = std::make_shared<OutputMeshDiscontinuous>(*mesh_);
	} else {
		output_mesh_ = std::make_shared<OutputMesh>(*mesh_);
	}
	output_mesh_->create_sub_mesh();
	if (!parallel) {
		output_mesh_->make_serial_master_mesh();
	} else {
		output_mesh_->make_parallel_master_mesh();
	}
	stream_->set_output_data_caches(output_mesh_);
}


typename OutputMeshBase::ErrorControlFieldFunc EquationOutput::select_error_control_field()
{
    std::string error_control_field_name = "";
    // Read optional error control field name
    auto it = stream_->get_output_mesh_record()->find<std::string>("error_control_field");
    if(it) error_control_field_name = *it;
    
    if(error_control_field_name!="")
    {
        FieldCommon* field =  this->field(error_control_field_name);
        // throw input exception if the field is unknown
        if(field == nullptr){
            THROW(FieldSet::ExcUnknownField()
                    << FieldCommon::EI_Field(error_control_field_name));
        }

        // throw input exception if the field is not scalar
        if( typeid(*field) == typeid(Field<3,FieldValue<3>::Scalar>) ) {

        	Field<3,FieldValue<3>::Scalar>* error_control_field = static_cast<Field<3,FieldValue<3>::Scalar>*>(field);
            DebugOut() << "Error control field for output mesh set: " << error_control_field_name << ".";
            auto lambda_function =
                [error_control_field](const Armor::array &point_list, const ElementAccessor<OutputMeshBase::spacedim> &elm, std::vector<double> &value_list)->void
                { error_control_field->value_list(point_list, elm, value_list); };

            OutputMeshBase::ErrorControlFieldFunc func = lambda_function;
            return func;

        }
        else{
            THROW(ExcFieldNotScalar()
                    << FieldCommon::EI_Field(error_control_field_name));
        }
    }
    return OutputMeshBase::ErrorControlFieldFunc();
}
