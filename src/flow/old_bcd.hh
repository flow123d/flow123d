/*
 * old_bcd.hh
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
 */

#ifndef OLD_BCD_HH_
#define OLD_BCD_HH_


#include <map>
#include <memory>
using namespace std;

#include "fields/field.hh"
#include "fields/field_elementwise.hh"

#include "mesh/mesh.h"
#include "input/accessors.hh"


/**
 * @brief Old BC setting system for backward compatibility.
 *
 * The class provides a singleton object for support of old input files with boundary conditions. The method @p read_flow,
 * reads the flow boundary file and creates FieldElementwise objects for the BC type, the pressure, the flux and the sigma.
 * More over it fills map that assigns boundary elements to the boundary IDs. This map si used by possible call of the
 * @p read_transport method that reads the transport boundary file and creates the FieldElementwise for the boundary concentrations.
 *
 * # Flow boundary file #
 * The boundary condition file for the flow problem contains section @p $BoundaryConditions in which the first line specifies the
 * number of BC lines. One BC line consists of the boudary face ID, the BC type (1 - Dirichlet, 2 - Neumann, 3 - Robin) followed by
 * the BC data: the pressure value (Dirichlet BC), or the flux value (Neumann BC), or the pressure and the coefficient (Robin BC).
 * After the boundary data there is face specification, it consists of the specification type (the value 2 is only supported),
 * the element ID and index of its side that corresponds to the face.
 * Sides are numbered from zero to @p dim-1, @p i-th side opposing to the @p i-th
 * node of the element according to the ordering in specification of the element in the mesh file.
 *
 * # Transport boundary file #
 * Concentration data can be specified only on faces that appears in the flow boundary file. The file has to contain section
 * @p $Transport_BCD which starts with number of transport BC lines. One transport BC line consists of line ID (ignored),
 * boundary ID (the BC line ID from the flow boundary file), boundary values for all substances used by the transport module.
 *
 *
 */
class OldBcdInput {
public:

	/**
	 * Use this to declare the key with filename.
	 */
	static string flow_old_bcd_file_key() {
		return "flow_old_bcd_file";
	}
	static string transport_old_bcd_file_key() {
		return "transport_old_bcd_file";
	}

	typedef FieldElementwise<3, FieldValue<3>::Scalar> FieldScalar;
	typedef FieldElementwise<3, FieldValue<3>::Enum> FieldEnum;
	typedef FieldElementwise<3, FieldValue<3>::Vector> FieldVector;

	typedef Field<3, FieldValue<3>::Scalar> Field_Scalar;
	typedef Field<3, FieldValue<3>::Enum> Field_Enum;
	typedef Field<3, FieldValue<3>::Vector> Field_Vector;

	typedef Field_Scalar::FieldBasePtr FieldBaseScalar;
	typedef Field_Enum::FieldBasePtr FieldBaseEnum;
	typedef Field_Vector::FieldBasePtr FieldBaseVector;

    shared_ptr<FieldEnum>	flow_type;
    shared_ptr<FieldScalar>  flow_pressure;
    shared_ptr<FieldScalar>  flow_flux;
    shared_ptr<FieldScalar>  flow_sigma;
    shared_ptr<FieldVector>  trans_conc;

    static OldBcdInput * instance();

    // hooks
    static FieldBaseEnum flow_type_hook(Input::Record rec, const FieldCommonBase &field) {
    	OldBcdInput *old_bcd = OldBcdInput::instance();
    	auto field_ptr=Field_Enum::read_field_descriptor(rec, field);
    	if (field_ptr) return field_ptr;
    	else {
    		old_bcd->read_flow_record(rec, field);
    		return old_bcd->flow_type;
    	}
    }
    static FieldBaseScalar flow_pressure_hook(Input::Record rec, const FieldCommonBase &field) {
    	OldBcdInput *old_bcd = OldBcdInput::instance();
    	auto field_ptr= Field_Scalar::read_field_descriptor(rec, field);
    	if (field_ptr) return field_ptr;
    	else {
    		old_bcd->read_flow_record(rec, field);
    		return old_bcd->flow_pressure;
    	}
    }
    static FieldBaseScalar flow_flux_hook(Input::Record rec, const FieldCommonBase &field) {
    	OldBcdInput *old_bcd = OldBcdInput::instance();
    	auto field_ptr = Field_Scalar::read_field_descriptor(rec, field);
    	if (field_ptr) return field_ptr;
    	else {
    		old_bcd->read_flow_record(rec, field);
    		return old_bcd->flow_flux;
    	}
    }
    static FieldBaseScalar flow_sigma_hook(Input::Record rec, const FieldCommonBase &field)  {
    	OldBcdInput *old_bcd = OldBcdInput::instance();
    	auto field_ptr = Field_Scalar::read_field_descriptor(rec, field);
    	if (field_ptr) return field_ptr;
    	else {
    		old_bcd->read_flow_record(rec, field);
    		return old_bcd->flow_sigma;
    	}
    }
    static FieldBaseVector trans_conc_hook(Input::Record rec, const FieldCommonBase &field)  {
    	OldBcdInput *old_bcd = OldBcdInput::instance();
    	auto field_ptr = Field_Vector::read_field_descriptor(rec,field);
    	if (field_ptr) return field_ptr;
    	else {
    		old_bcd->read_transport_record( rec, field);
    		return old_bcd->trans_conc;
    	}
    }

    void read_flow_record(Input::Record rec, const FieldCommonBase &field) {
    	FilePath bcd_file;
    	if (rec.opt_val(flow_old_bcd_file_key(), bcd_file)
    			&& string(bcd_file) != flow_input_file_) {
    		ASSERT(field.mesh(),"Null mesh pointer.");
    		read_flow(*(field.mesh()), bcd_file);
    		flow_input_file_ = string(bcd_file);
    	}
    }


    inline void read_transport_record(Input::Record rec, const FieldCommonBase &field);


    /**
     * Create flow_* fields from given input file.
     */
    void read_flow(const Mesh &mesh, const FilePath &flow_bcd);

    /**
     * Create trans_conc field from given input file.
     */
    void read_transport(unsigned int n_substances, const FilePath &transport_bcd);

    /// Maps ID to index of corresponding BC element.
    map<unsigned int, unsigned int> id_2_bcd_;

private:
 /*
    template <int spacedim, class Value>
    void set_all( Field<spacedim,Value> &target, const Mesh *mesh);

    template <int spacedim, class Value>
    void set_field( Field<spacedim,Value> &target, unsigned int bcd_ele_idx, typename Value::return_type &val);
*/
    const Mesh *mesh_;
    Region  some_bc_region_;

    string flow_input_file_;
    string transport_input_file_;


};



void OldBcdInput::read_transport_record(Input::Record rec, const FieldCommonBase &field) {
	FilePath bcd_file;
	if (rec.opt_val(transport_old_bcd_file_key(), bcd_file)
			&& string(bcd_file) != transport_input_file_) {
		ASSERT(field.mesh(),"Null mesh pointer.");
		read_transport( field.n_comp(), bcd_file);
		transport_input_file_ = string(bcd_file);
	}
}

#endif /* OLD_BCD_HH_ */
