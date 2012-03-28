/*
 * sources.cc
 *
 *  Created on: Feb 9, 2012
 *      Author: jiri
 */
#include <petscmat.h>
#include "transport/sources.hh"

#include "system/system.hh"
//#include "system/math_fce.h"
#include "mesh/mesh.h"
#include "transport/transport.h"
//#include "materials.hh"
#include "io/read_ini.h"
#include "la/distribution.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include "xio.h"


//=============================================================================
// TRANSPORT SOURCES CONSTRUCTOR
//=============================================================================
TransportSources::TransportSources(ConvectionTransport &transport)
{
	convectiontransport = &transport;
	alloc_sources_vectors();
	read_concentration_sources();
}

//=============================================================================
// READ CONCENTRATION SOURCES
//=============================================================================
void TransportSources::read_concentration_sources() {
    F_ENTRY;

        FILE	*in;		  // input file
		char     line[ LINE_SIZE ]; // line of data file
		int sbi,index, eid, i,n_sources, global_idx;

        std::string concentration_sources_fname = OptGetFileName("Transport", "Sources", "\\");
		in = xfopen( IONameHandler::get_instance()->get_input_file_name(concentration_sources_fname), "rt" );
		skip_to( in, "$TransportSources" );
		xfgets( line, LINE_SIZE - 2, in );
		n_sources = atoi( xstrtok( line) );
		INPUT_CHECK(!( n_sources < 1 ),"Number of concentration sources < 1 in function read_concentration_sources()\n");
	    INPUT_CHECK(!( convectiontransport->mesh_->n_elements() != n_sources),"Different number of elements and concentration sources\n");


	    for (i = 0; i < n_sources; i++) {
	    	xfgets( line, LINE_SIZE - 2, in );
	    	ASSERT(!(line == NULL),"NULL as argument of function read_concentration_sources()\n");
	    	eid    = atoi( xstrtok(line) );
	    	global_idx = convectiontransport->row_4_el[convectiontransport->mesh_->element.find_id(eid).index()];
	    	if ( convectiontransport->el_ds->is_local(global_idx) ) {
	    		index = global_idx - convectiontransport->el_ds->begin();
	    		for( sbi = 0; sbi < convectiontransport->n_substances; sbi++ ){
	    			sources_density[ sbi ][index] = atof( xstrtok( NULL) );
	    			sources_sigma[ sbi ][index] = atof( xstrtok( NULL) );
	    			sources_conc[ sbi ][index] = atof( xstrtok( NULL) );
	    		}
	    	}
	    }

		xfclose( in );
}
//=============================================================================
// ALLOCATE CONCENTRATION SOURCES VECTORS
//=============================================================================
void TransportSources::alloc_sources_vectors() {

    int i, j, sbi, ierr, rank, np, n_subst;

    n_subst = convectiontransport->n_substances;

    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    sources_density = (double**) xmalloc(n_subst * sizeof(double*));
    sources_sigma = (double**) xmalloc( n_subst* sizeof(double*));
    sources_conc = (double**) xmalloc(n_subst * sizeof(double*));
    sources_corr = (double**) xmalloc(n_subst * sizeof(double*));

    for (sbi = 0; sbi < n_subst; sbi++){
        sources_density[sbi] = (double*) xmalloc(convectiontransport->el_ds->lsize() * sizeof(double));
        sources_sigma[sbi] = (double*) xmalloc(convectiontransport->el_ds->lsize() * sizeof(double));
        sources_conc[sbi] = (double*) xmalloc(convectiontransport->el_ds->lsize() * sizeof(double));
        sources_corr[sbi] = (double*) xmalloc(convectiontransport->el_ds->lsize() * sizeof(double));
    }
/*
    for (sbi = 0; sbi < n_subst; sbi++)
        for (i = 0; i < convectiontransport->el_ds->size(); i++){
            sources_density[sbi][i] = 0.0;
            sources_sigma[sbi][i] = 0.0;
            sources_conc[sbi][i] = 0.0;
            sources_corr[sbi][i] = 0.0;
        }
*/
    vsources_density = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vsources_sigma = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vsources_conc = (Vec*) xmalloc(n_subst * (sizeof(Vec)));
    vsources_corr = (Vec*) xmalloc(n_subst * (sizeof(Vec)));

    for (sbi = 0; sbi < n_subst; sbi++){
    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, convectiontransport->el_ds->lsize(), convectiontransport->mesh_->n_elements(),
    		sources_density[sbi],&vsources_density[sbi]);
    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, convectiontransport->el_ds->lsize(), convectiontransport->mesh_->n_elements(),
    		sources_sigma[sbi],&vsources_sigma[sbi]);
    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, convectiontransport->el_ds->lsize(), convectiontransport->mesh_->n_elements(),
    		sources_conc[sbi],&vsources_conc[sbi]);
    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, convectiontransport->el_ds->lsize(), convectiontransport->mesh_->n_elements(),
    		sources_corr[sbi],&vsources_corr[sbi]);
    }
}
//=============================================================================
// COMPUTE SOURCES
//=============================================================================
void TransportSources::compute_concentration_sources(int sbi) {

	int i_loc;
		for (i_loc = 0; i_loc < convectiontransport->el_ds->lsize(); i_loc++){
			if( (sources_conc[ sbi ][i_loc] - convectiontransport->conc[MOBILE][sbi][i_loc]) > 0.0)
        		sources_corr[sbi][i_loc] = sources_density[ sbi ][i_loc] +
        		(sources_conc[ sbi ][i_loc]- convectiontransport->conc[MOBILE][sbi][i_loc])*sources_sigma[ sbi ][i_loc];
        	else
        		sources_corr[sbi][i_loc] = sources_density[ sbi ][i_loc];
        }

}


