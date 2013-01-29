/*
 * old_bcd.cc
 *
 *  Created on: Jan 28, 2013
 *      Author: jb
 */


#include "flow/old_bcd.hh"
#include "fields/field_elementwise.hh"

#include "system/tokenizer.hh"
#include "boost/lexical_cast.hpp"

OldBcdInput * OldBcdInput::instance() {
    static OldBcdInput *obcd = new OldBcdInput;
    return obcd;
}


template <int spacedim, class Value>
void OldBcdInput::set_all( Field<spacedim,Value> &target, Mesh *mesh) {
    FieldElementwise<spacedim, Value> *in_field=new FieldElementwise<spacedim, Value>(target.n_comp());
    for(unsigned int i=0; i<Region::db().size(); i++)
        if ( Region(i).is_boundary() ) target.set_field(Region(i), in_field);
    target.set_mesh(mesh);

}

template <int spacedim, class Value>
void OldBcdInput::set_field( Field<spacedim,Value> &target, unsigned int bcd_ele_idx, typename Value::return_type &val, Region bc_reg) {
    static_cast<FieldElementwise<spacedim, Value> *>(target(bc_reg))->set_data_row(bcd_ele_idx, val);
}


#define DIRICHLET   1
#define NEUMANN     2
#define NEWTON      3

void OldBcdInput::read_flow(const FilePath &flow_bcd,
        Mesh *mesh,
        Field<3,FieldValue<3>::Enum > &flow_type,
        Field<3,FieldValue<3>::Scalar > &flow_pressure,
        Field<3,FieldValue<3>::Scalar > &flow_flux,
        Field<3,FieldValue<3>::Scalar > &flow_sigma)
{
    using namespace boost;

    mesh_=mesh;
/*
 * - read one flow file, fill fields, make ID list
 * - read second file, check IDs agains ID list, fill fields
 */
    set_all(flow_type, mesh);
    set_all(flow_pressure, mesh);
    set_all(flow_flux, mesh);
    set_all(flow_sigma, mesh);

    // get any boundary region
    Region bc_reg;
    for(unsigned int i_reg=0; i_reg<Region::db().size(); i_reg++)
        if (Region(i_reg).is_boundary()) {
            bc_reg=Region(i_reg);
            break;
        }
    if (! bc_reg.is_valid()) xprintf(PrgErr, "No boundary region.\n");

    Tokenizer tok(flow_bcd);
    try {
        double scalar, flux, sigma;

        tok.skip_to("$BoundaryConditions");
        tok.next_line(false);
        unsigned int n_boundaries = lexical_cast<unsigned int>(*tok); ++tok;
        bcd_ids_.resize(n_boundaries);

        for(unsigned int i_bcd=0; i_bcd < n_boundaries; i_bcd++) {
            tok.next_line();

            bcd_ids_[i_bcd] = lexical_cast<unsigned int>(*tok); ++tok;
            unsigned int type  = lexical_cast<unsigned int>(*tok); ++tok;

            switch( type ) {
                case DIRICHLET:
                    scalar = lexical_cast<double>(*tok); ++tok;
                    flux = 0.0;
                    sigma = 0.0;
                    break;
                case NEUMANN:
                    flux   = lexical_cast<double>(*tok); ++tok;
                    sigma = 0.0;
                    scalar = 0.0;
                    break;
                case NEWTON:
                    scalar = lexical_cast<double>(*tok); ++tok;
                    sigma  = lexical_cast<double>(*tok); ++tok;
                    flux = 0.0;
                    break;
                default :
                    xprintf(UsrErr,"Unknown type of boundary condition - cond # %d, type %c\n", bcd_ids_[i_bcd], type );
                    break;
            }

            unsigned int where  = lexical_cast<unsigned int>(*tok); ++tok;

            unsigned int eid, sid, bc_ele_idx;
            ElementIter ele;
            SideIter sde;

            switch( where ) {
                case 2: // SIDE_EL - BC given by element and its local side number
                    eid = lexical_cast<unsigned int>(*tok); ++tok;
                    sid = lexical_cast<unsigned int>(*tok); ++tok;

                    // find and set the side
                    ele = mesh->element.find_id( eid );
                    if( sid < 0 || sid >= ele->n_sides() )
                         xprintf(UsrErr,"Boundary %d has incorrect reference to side %d\n", bcd_ids_[i_bcd], sid );
                    bc_ele_idx = mesh->bc_elements.index( ele->side(sid) -> cond()->element() );
                    set_field(flow_type,     bc_ele_idx, type,   bc_reg);
                    set_field(flow_pressure, bc_ele_idx, scalar, bc_reg);
                    set_field(flow_flux,     bc_ele_idx, flux,   bc_reg);
                    set_field(flow_sigma,    bc_ele_idx, sigma,  bc_reg);
                    break;
                case 3: // SIDE_E - BC given only by element, apply to all its sides

                    xprintf(UsrErr, "Element only BCD are not supported.\n");
                    /*
                    eid = atoi( xstrtok( NULL) );

                    // find and set all exterior sides, possibly add more boundaries
                    ele = mesh->element.find_id( eid );
                    n_exterior=0;
                    FOR_ELEMENT_SIDES(ele, li) {
                        sde = ele->side( li );
                        if ( bcd=sde->cond() ) {

                            if (n_exterior > 0) {
                                xprintf(UsrErr, "Implicitly setting BC %d on more then one exterior sides of the element %d.\n", bcd_id, eid);
                                //BoundaryFullIter new_bcd = mesh->boundary.add_item();
                                //*new_bcd = *bcd;
                                //bcd=new_bcd;
                            }
                            bcd->type = type;
                            bcd->flux = flux;
                            bcd->scalar = scalar;
                            bcd->sigma = sigma;
                            n_exterior++;
                        }
                    }
                    */
                    break;
                default:
                    xprintf(UsrErr,"Unknown entity for boundary condition - cond # %d, ent. %c\n", bcd_ids_[i_bcd], where );
                    break;
            }
            unsigned int n_tags  = lexical_cast<unsigned int>(*tok); ++tok;
            while (n_tags>0) ++tok, --n_tags; // skip remaining tags



            // There was possibility to set group IDs for boundary faces (like regions for boundary elements)
            // It is probably not used so we do not implement it as it is DEPRECATED
            /*
            n_tags  = atoi( xstrtok( NULL) );
            if( n_tags > 0 ) {
                int group_id = atoi( xstrtok( NULL) );
                flow::VectorId<int>::FullIter group_iter( mesh->bcd_group_id.find_id(group_id) );

                if ( group_iter == mesh->bcd_group_id.end() ) {
                    // not found -> create new group
                    group_iter = mesh->bcd_group_id.add_item(group_id);
                }
                bcd->group = group_iter.index();   // in fact we do not use integres stored in the vector, but we use index
            }
            */

        }
    } catch (bad_lexical_cast &) {
        xprintf(UsrErr, "Wrong format of number, %s.\n", tok.position_msg().c_str());
    } // flow bcd reader
}

void OldBcdInput::read_transport(const FilePath &transport_bcd,
            Field<3,FieldValue<3>::Enum > &trans_type,
            Field<3,FieldValue<3>::Vector > &trans_conc)
{
    set_all(trans_type, mesh_);
    set_all(trans_conc, mesh_);




}
