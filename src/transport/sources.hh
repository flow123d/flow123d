/*
 * sources.hh
 *
 *  Created on: Feb 9, 2012
 *      Author: jiri
 */

#ifndef SOURCES_HH_
#define SOURCES_HH_

class ConvectionTransport;
#include "transport/transport.h"


/**
 * Class TransportSources is separating class for Newman/Newtons sources.
 * Class provides input data reading and individual sources computation for every substance.
 */


class TransportSources
{
	friend class ConvectionTransport;
public:
	/**
	 * Constructor
	 */
	TransportSources(ConvectionTransport &convection);//(ConvectionTransport &convection, int n_subst);
	/**
	 * Initial allocating method
	 */
	void alloc_sources_vectors();
	/**
	 * Input data reading method
	 */
	void read_concentration_sources();
	/**
	 * Main computation method
	 */
	void compute_concentration_sources(int sbi);
private:
	ConvectionTransport *convectiontransport;
	/**
	 * Input data fields
	 */
    double **sources_density;
    double **sources_sigma;
    double **sources_conc;
	/**
	 * Correction vector
	 */
    double **sources_corr;
	/**
	 * Encapsulation structures for data
	 */
    Vec *vsources_density;
    Vec *vsources_sigma;
    Vec *vsources_conc;
    Vec *vsources_corr;
};




#endif /* SOURCES_HH_ */
