/*
 * sources.hh
 *
 *  Created on: Feb 9, 2012
 *      Author: jiri
 */

#ifndef SOURCES_HH_
#define SOURCES_HH_

#include <string>
#include "la/distribution.hh"

using namespace std;
class Mesh;

/**
 * Class TransportSources is separating class for Newman/Newtons sources.
 * Class provides input data reading and individual sources computation for every substance.
 */


class TransportSources
{
public:
	/**
	 * Constructor
	 */
	TransportSources(unsigned int n_subst, const  Distribution &el_distr);
	/**
	 * Initial allocating method
	 */
	void alloc_sources_vectors();

	/**
	 * Input data reading method
	 */
	void read_concentration_sources(const string &sources_fname,  int *row_4_el, Mesh *mesh);

	/**
	 * Main computation method
	 */
	Vec compute_concentration_sources(unsigned int subst_i, double *conc );
private:

	unsigned int n_subst_;
	Distribution el_distr_;

	/**
	 * Input data fields
	 */
    double **sources_density;
    double **sources_sigma;
    double **sources_conc;
	/**
	 * Correction vector
	 */
    double *sources_corr;
    Vec v_sources_corr;
};




#endif /* SOURCES_HH_ */
