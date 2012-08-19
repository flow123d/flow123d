/*
 * sources.cc
 *
 *  Created on: Feb 9, 2012
 *      Author: jiri
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <petscvec.h>

#include "transport/sources.hh"

#include "system/system.hh"
#include "mesh/mesh.h"
#include "la/distribution.hh"
#include "system/xio.h"


//=============================================================================
// TRANSPORT SOURCES CONSTRUCTOR
//=============================================================================
TransportSources::TransportSources(unsigned int n_subst, const  Distribution &el_distr)
: n_subst_(n_subst),
  el_distr_(el_distr)
{
	alloc_sources_vectors();
}


//=============================================================================
// ALLOCATE CONCENTRATION SOURCES VECTORS
//=============================================================================
void TransportSources::alloc_sources_vectors() {

    int i, j, sbi, ierr, rank, np, n_subst;


    MPI_Barrier(PETSC_COMM_WORLD);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &np);

    unsigned int alloc_size = n_subst_ * ( sizeof(double*) + el_distr_.lsize() * sizeof(double) );

    sources_density = (double**) xmalloc(3* alloc_size);
    sources_sigma = sources_density + 1;
    sources_conc = sources_density + 2;

    double *ptr = (double *) sources_conc + 1;

    for (sbi = 0; sbi < n_subst_; sbi++){
        sources_density[sbi] = ptr; ptr += el_distr_.lsize();
        sources_sigma[sbi] = ptr; ptr += el_distr_.lsize();
        sources_conc[sbi] = ptr; ptr += el_distr_.lsize();
    }

    sources_corr = new double[el_distr_.lsize()];

    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, el_distr_.lsize(), PETSC_DECIDE,
            sources_corr, &v_sources_corr);


/*
    for (sbi = 0; sbi < n_subst; sbi++)
        for (i = 0; i < convectiontransport->el_ds->size(); i++){
            sources_density[sbi][i] = 0.0;
            sources_sigma[sbi][i] = 0.0;
            sources_conc[sbi][i] = 0.0;
            sources_corr[sbi][i] = 0.0;
        }
*/
    /*
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
    }*/
}


//=============================================================================
// READ CONCENTRATION SOURCES
//=============================================================================
void TransportSources::read_concentration_sources(const string &sources_fname,  int *row_4_el, Mesh *mesh) {
    F_ENTRY;

        FILE	*in;		  // input file
		char     line[ LINE_SIZE ]; // line of data file
		int sbi,index, eid, i,n_sources, global_idx;


		in = xfopen( sources_fname, "rt" );
		skip_to( in, "$TransportSources" );
		xfgets( line, LINE_SIZE - 2, in );
		n_sources = atoi( xstrtok( line) );
		INPUT_CHECK( n_sources >= 1 , "Number of concentration sources < 1 in function read_concentration_sources()\n");
	    INPUT_CHECK( el_distr_.size() == n_sources, "Different number of elements and concentration sources\n");


	    for (i = 0; i < n_sources; i++) {
	    	xfgets( line, LINE_SIZE - 2, in );
	    	ASSERT(!(line == NULL),"NULL as argument of function read_concentration_sources()\n");
	    	eid    = atoi( xstrtok(line) );
	    	global_idx = row_4_el[mesh->element.find_id(eid).index()];

	    	if ( el_distr_.is_local(global_idx) ) {
	    		index = global_idx - el_distr_.begin();
	    		for( sbi = 0; sbi < n_subst_; sbi++ ){
	    			sources_density[ sbi ][index] = atof( xstrtok( NULL) );
	    			sources_sigma[ sbi ][index] = atof( xstrtok( NULL) );
	    			sources_conc[ sbi ][index] = atof( xstrtok( NULL) );
	    		}
	    	}
	    }

		xfclose( in );
}



//=============================================================================
// COMPUTE SOURCES
//=============================================================================
Vec TransportSources::compute_concentration_sources(unsigned int subst_i, double *conc) {

    double conc_diff;
    for (int i_loc = 0; i_loc < el_distr_.lsize(); i_loc++) {

        conc_diff = sources_conc[subst_i][i_loc] - conc[i_loc];
        if ( conc_diff > 0.0)
            sources_corr[i_loc] = sources_density[subst_i][i_loc] + conc_diff * sources_sigma[subst_i][i_loc];
        else
            sources_corr[i_loc] = sources_density[subst_i][i_loc];

       // cout << i_loc << " c:" << conc[i_loc] << " sc:" << sources_conc[subst_i][i_loc] << " sd:"
       //      << sources_density[subst_i][i_loc] << " ss:" << sources_sigma[subst_i][i_loc] << " cr:"
       //      << sources_corr[i_loc] << endl;
    }

    return v_sources_corr;

}


