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
//#include "equation.hh"

class TransportSources
{
	friend class ConvectionTransport;
public:
	TransportSources(ConvectionTransport &convection);//(ConvectionTransport &convection, int n_subst);
	void alloc_sources_vectors();
	void read_concentration_sources();
	void compute_concentration_sources(int sbi);
private:
	ConvectionTransport *convectiontransport;
    Vec *vsources_density;
    Vec *vsources_sigma;
    Vec *vsources_conc;
    Vec *vsources_corr;
 //   Vec *vcumulative_corr;

    double **sources_density;
    double **sources_sigma;
    double **sources_conc;

    double **sources_corr;

    double **cumulative_corr;

};




#endif /* SOURCES_HH_ */
