/*
 * equation_output.hh
 *
 *  Created on: Jul 8, 2016
 *      Author: jb
 */

#ifndef SRC_FIELDS_EQUATION_OUTPUT_HH_
#define SRC_FIELDS_EQUATION_OUTPUT_HH_

#include <memory>                  // for shared_ptr
#include <string>                  // for string
#include <unordered_map>           // for unordered_map
#include <unordered_set>           // for unordered_set
#include "fields/field_common.hh"  // for FieldCommon, FieldCommon::EI_Field
#include "fields/field_set.hh"     // for FieldSet
#include "fields/field_values.hh"  // for FieldValue, FieldValue<>::Scalar
#include "io/output_time_set.hh"   // for OutputTimeSet
#include "io/output_mesh.hh"
#include "system/exceptions.hh"    // for ExcStream, operator<<, DECLARE_EXC...
#include "tools/time_marks.hh"     // for TimeMark, TimeMark::Type

class OutputTime;
class TimeGovernor;
class TimeStep;
class DOFHandlerMultiDim;
namespace Input {
	class Record;
	namespace Type {
		class Instance;
		class Record;
                class Selection;
	}
}
template<unsigned int dim> class AssemblyOutputElemData;
template<unsigned int dim> class AssemblyOutputNodeData;
template< template<IntDim...> class DimAssembly> class GenericAssembly;


/**
 * A class  responsible for check for output times of individual fields
 * and store their values into the connected output stream.
 */
class EquationOutput : public FieldSet {
public:

    DECLARE_EXCEPTION(ExcFieldNotScalar, << "Field '" << FieldCommon::EI_Field::qval
                                         << "' is not scalar in spacedim 3.");

    /// Configuration of output of one field. Pair of OutputTimeSet and DiscreteSpaces.
    struct FieldOutputConfig {
        OutputTimeSet output_set_;                    ///< Set of output times.
        OutputTime::DiscreteSpaceFlags space_flags_;  ///< Array of used DiscreteSpaces
    };

    /**
     * Input type of the configuration record.
     */
    static Input::Type::Record &get_input_type();
    
    /// Default constructor
    EquationOutput();

    /// Destructor
    ~EquationOutput();

    /**
     * Make Input::Type for the output record. Particular selection of output fields is created
     * from the contents of *this FieldSet using provided equation name and additional description.
     */
    const Input::Type::Instance &make_output_type(const string &equation_name, const string &aditional_description = "");
    /**
     * Make Input::Type for the output record.
     * This function enables creating output record for a field set record with additional keys provided in @p in_rec.
     */
    const Input::Type::Instance &make_output_type_from_record(Input::Type::Record &in_rec,
                                                              const string &equation_name,
                                                              const string &aditional_description = "");
    
    /**
     * Setup the object. Set output stream for field and observe output, input record for configuration of the object and
     * TimeGovernor. The time governor is used to get the equation time mark type, the initial and the end time of the equation.
     */
    void initialize(std::shared_ptr<OutputTime> stream, Mesh *mesh, Input::Record in_rec, const TimeGovernor & tg);

    /**
     * Returns true if @param field is marked for output in the given time @param step.
     */
    bool is_field_output_time(const FieldCommon &field, TimeStep step) const;

    /**
     * Performs output of the fields marked for output in the time @param step.
     */
    void output(TimeStep step);

    /// Selects the error control field out of output field set according to input record.
    typename OutputMeshBase::ErrorControlFieldFunc select_error_control_field();
    
private:
    /**
     * Creates output selection from the field set.
     */
    const Input::Type::Selection& create_output_field_selection(const string &equation_name,
                                                                const string &additional_description);
    
    /**
     * Read from the input, set output times and time marks. Must be called after set_stream.
     * TODO: add output_stream times. Optional or always?
     */
    void read_from_input(Input::Record in_rec, const TimeGovernor & tg);

    /**
     * Add a time grid to the common_output_times.
     */
    void add_output_times(double begin, double step, double end);


    /**
     * Create the output mesh of \p stream_ OutputTime object. The field set passed in is used
     * to select the field used for adaptivity of the output mesh.
     */
    void make_output_mesh(bool parallel);

    /// Initialize data of Field given by passed Input::Record
    void init_field_item(Input::Iterator<Input::Record> it, const TimeGovernor & tg);

    /// output stream (may be shared by more equation)
    std::shared_ptr<OutputTime> stream_;
    /// The time mark type of the equation.
    TimeMark::Type equation_type_;
    /// The fixed time mark type of the equation.
    TimeMark::Type equation_fixed_type_;
    /// The time set used for the fields without explicit time set.
    OutputTimeSet common_output_times_;

    /// Time sets of individual fields.
    std::unordered_map<string, FieldOutputConfig> field_output_times_;

    /// Set of observed fields. The observe points are given within the observe stream.
    std::unordered_set<string> observe_fields_;

    /**
     * Set of interpolations which are used in performed fields.
     *
     * Allow determine type of output mesh.
     */
    std::set<OutputTime::DiscreteSpace> used_interpolations_;

    /**
     * Cached pointer at computational mesh.
     */
    Mesh *mesh_;

    /// Output mesh.
    std::shared_ptr<OutputMeshBase> output_mesh_;

    /// Objects for distribution of dofs.
    std::shared_ptr<DOFHandlerMultiDim> dh_;
    std::shared_ptr<DOFHandlerMultiDim> dh_node_;

    /// general assembly objects, hold assembly objects of appropriate dimension
    GenericAssembly< AssemblyOutputElemData > * output_elem_data_assembly_;
    GenericAssembly< AssemblyOutputNodeData > * output_node_data_assembly_;
    GenericAssembly< AssemblyOutputNodeData > * output_corner_data_assembly_;

};


#endif /* SRC_FIELDS_EQUATION_OUTPUT_HH_ */
