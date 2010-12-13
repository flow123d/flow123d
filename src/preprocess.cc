
#include "constantdb.h"
#include "mesh/ini_constants_mesh.hh"

#include "transport.h"

#include "system.hh"
#include "elements.h"
#include "problem.h"
#include "mesh.h"
#include "initials.h"
#include "preprocess.h"

static void make_ele_initials(struct Problem*);

//=============================================================================
//SET VALUES FROM INITIAL CONDITIONS ETC...
//=============================================================================
void preprocess( struct Problem *problem )
{
    F_ENTRY;

	xprintf( Msg, "Preprocessing...")/*orig verb 2*/;
    if (ConstantDB::getInstance()->getInt("Problem_type") == UNSTEADY_SATURATED) // presunout sem
    	make_ele_initials( problem );            // pp koncentraci

                /*init_concentration_list( mesh );
	FOR_CONCENTRATIONS( con ) {
		xfgets( line, LINE_SIZE - 2, in );
		parse_concentration_line( con, line );
	}
	xfclose( in );*/

	//settings for variable density


}
//=============================================================================
//SET ELEMENT'S VALUES FROM INITIAL CONDITIONS
//=============================================================================
void make_ele_initials(struct Problem *problem)
{
	xprintf( Msg, "   Setting element's values from initial conditions...")/*orig verb 2*/;

        Mesh* mesh = (Mesh*)ConstantDB::getInstance()->getObject(MESH::MAIN_INSTANCE);
        ElementIter ele;
        FOR_ELEMENTS( ele ){
                ele->pscalar = ele->initial->epress;
        }
	xprintf( Msg, "O.K.\n")/*orig verb 2*/;
        return;
}


