/*
 * bc_table.hh
 *
 *  Created on: Oct 1, 2012
 *      Author: jb
 */

#ifndef BC_TABLE_HH_
#define BC_TABLE_HH_

#include <vector>
#include <string>

#include "input/input_type.hh"

using namespace std;

template <class BC_Type>
class BoundarySegmentData {

    inline BC_Type *segment_data(ElementFullIter &ele) {
        return segmnet_data_[ bc_data_index_.value(ele.index()) ];
    }

private:
    DiscreteFunction bc_data_index_;
    vector<BC_Type *> segment_data_;

};

template <class BC_Type>
class BoundaryData {
public:
    static Input::Type::Record &get_input_type();

    BoundaryData( Mesh &mesh );

    void init_from_input( Input::Record &in_rec );

    /**
     * Obsolete. Only for particular BC_Type = DarcyFlowMH_BC
     */
    void init_from_flow_bc_file( FilePath &bc_file );

    inline BC_Type *bc_for_segment_and_element( ElementFullIter ele) {
        return bc_data_[ele->bc_segment()].segment_data( ele );
    }

private:
    std::vector< BoundarySegmentData<BC_Type> > bc_data_;
};

/********************************************** Implementation */

template <class BC_Type>
Input::Type::Record &BoundaryData<BC_Type>::get_input_type() {
    using namespace Input::Type;
    static Input::Type::Record rec("BoundaryData_"+BC_Type::get_input_type().name(), "Boundary data for equation: " + BC_Type::get_input_type().name());

    if (! rec.is_finished() ) {
        rec.declare_key("segment", Integer(0), Default("0"),
                "Apply boundary data to this boundary segment. Default zero segment denotes whole boundary.");
        rec.declare_key("bc_selector", DiscreteFunction::get_input_type(), Default("0"),
                "Discrete function that select boundary condition from array 'bc_data' for particular boundary element.");
        rec.declare_key("bc_data", Array( BC_Type::get_input_type() ), Default::obligatory(),
                "Array of boundary conditions, usually just one condition on the bc segment.");
    }
}



template <class BC_Type>
void  BoundaryData<BC_Type>::init_from_flow_bc_file( FilePath &bc_file ) {
    FILE    *in;          // input file
    char     line[ LINE_SIZE ]; // line of data file
//  int where;
    int bcd_id, n_tags;
    BoundaryFullIter bcd = BOUNDARY_NULL(mesh);
    ElementFullIter ele = ELEMENT_FULL_ITER_NULL(mesh);

    ASSERT(!( mesh == NULL ),"NULL as argument of function read_boundary_list()\n");
    xprintf( Msg, "Reading boundary conditions...")/*orig verb 2*/;

    in = xfopen( boundary_filename, "rt" );
    skip_to( in, "$BoundaryConditions" );
    xfgets( line, LINE_SIZE - 2, in );

    int n_boundaries = atoi( xstrtok( line) );
    INPUT_CHECK( n_boundaries >= 1 ,"Number of bounaries < 1 in function read_boundary_list()\n");
    mesh->boundary.reserve(n_boundaries);

    int group_number=0;

    for(int i_bcd=0; i_bcd < n_boundaries; i_bcd++) {
        // Read one line
        xfgets( line, LINE_SIZE - 2, in );
        // Parse the line
        bcd_id    = atoi( xstrtok( line) );
        bcd = mesh->boundary.add_item(bcd_id);
        // DBGMSG("boundary id: %d \n",bcd_id);

        bcd->type  = atoi( xstrtok( NULL) );

        // physical data - should be moved to water_linsys
        switch( bcd->type ) {
            case DIRICHLET:
                bcd->scalar = atof( xstrtok( NULL) );
                break;
            case NEUMANN:
                bcd->flux   = atof( xstrtok( NULL) );
                break;
            case NEWTON:
                bcd->scalar = atof( xstrtok( NULL) );
                bcd->sigma  = atof( xstrtok( NULL) );
                            break;
            default :
                xprintf(UsrErr,"Unknown type of boundary condition - cond # %d, type %c\n", bcd_id, bcd->type );
                break;
        }

        unsigned int where  = atoi( xstrtok( NULL) );
        int eid, sid, n_exterior;
        SideIter sde;

        switch( where ) {
            case 2: // SIDE_EL - BC given by element and its local side number
                eid = atoi( xstrtok( NULL) );
                sid = atoi( xstrtok( NULL) );

                // find and set the side
                ele = mesh->element.find_id( eid );
                if( sid < 0 || sid >= ele->n_sides() )
                     xprintf(UsrErr,"Boundary %d has incorrect reference to side %d\n", bcd_id, sid );
                sde = ele->side( sid );
                ele->boundaries_[ sid ] = bcd;
                bcd->side=sde;

                break;
            case 3: // SIDE_E - BC given only by element, apply to all its sides
                eid = atoi( xstrtok( NULL) );

                // find and set all exterior sides, possibly add more boundaries
                ele = mesh->element.find_id( eid );
                n_exterior=0;
                FOR_ELEMENT_SIDES(ele, li) {
                    sde = ele->side( li );
                    if ( sde->is_external() ) {
                        if (n_exterior > 0) {
                            xprintf(UsrErr, "Implicitly setting BC %d on more then one exterior sides of the element %d.\n", bcd_id, eid);
                            //BoundaryFullIter new_bcd = mesh->boundary.add_item();
                            //*new_bcd = *bcd;
                            //bcd=new_bcd;
                        }
                        ele->boundaries_[ sid ] = bcd;
                        bcd->side=sde;
                        n_exterior++;
                    }
                }

                break;
            default:
                xprintf(UsrErr,"Unknown entity for boundary condition - cond # %d, ent. %c\n", bcd_id, where );
                break;
        }
        //TODO: if group is necessary set it for all bcd in case where == SIDE_E
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


    }

    xfclose( in );
    xprintf( MsgVerb, " %d conditions readed. ", mesh->n_boundaries() )/*orig verb 4*/;
    xprintf( Msg, "O.K.\n")/*orig verb 2*/;
}

#endif /* BC_TABLE_HH_ */
