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


class OldBcdInput {
public:
    static OldBcdInput * instance();

    void read_flow(const FilePath &flow_bcd,
        Mesh *mesh,
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
    void set_all( Field<spacedim,Value> &target, Mesh *mesh);

    template <int spacedim, class Value>
    void set_field( Field<spacedim,Value> &target, unsigned int bcd_ele_idx, typename Value::return_type &val);

    Mesh *mesh_;
};


#endif /* OLD_BCD_HH_ */
