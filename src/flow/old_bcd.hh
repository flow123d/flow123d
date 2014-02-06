/*
 * old_bcd.hh
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
 */

#ifndef OLD_BCD_HH_
#define OLD_BCD_HH_

#include "fields/field_base.hh"
#include "mesh/mesh.h"
#include <map>

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
    static OldBcdInput * instance();

    void read_flow(const FilePath &flow_bcd,
        Field<3,FieldValue<3>::Enum > &flow_type,
        Field<3,FieldValue<3>::Scalar > &flow_pressure,
        Field<3,FieldValue<3>::Scalar > &flow_flux,
        Field<3,FieldValue<3>::Scalar > &flow_sigma);

    void read_transport(const FilePath &transport_bcd,
        Field<3,FieldValue<3>::Vector > &trans_conc);

    /// Maps ID to index of corresponding BC element.
    map<unsigned int, unsigned int> id_2_bcd_;

private:
    template <int spacedim, class Value>
    void set_all( Field<spacedim,Value> &target, const Mesh *mesh);

    template <int spacedim, class Value>
    void set_field( Field<spacedim,Value> &target, unsigned int bcd_ele_idx, typename Value::return_type &val);

    const Mesh *mesh_;
    Region  some_bc_region_;
};


#endif /* OLD_BCD_HH_ */
