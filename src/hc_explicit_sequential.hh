/*
 * hc_explicit_sequential.hh
 *
 *  Created on: Jul 14, 2011
 *      Author: jb
 */

#ifndef HC_EXPLICIT_SEQUENTIAL_HH_
#define HC_EXPLICIT_SEQUENTIAL_HH_

#include "main.h"

class DarcyFlowMH;
class DarcyFlowMHOutput;
class TimeMarks;
class Mesh;
class EquationBase;
class TransportBase;
class MaterialDatabase;

/**
 * @brief Class for solution of steady or unsteady flow with sequentially coupled explicit transport.
 *
 */
class HC_ExplicitSequential {
public:
    HC_ExplicitSequential(ProblemType problem_type);
    void run_simulation();
    ~HC_ExplicitSequential();

private:
    ProblemType type_;

    /// mesh common to darcy flow and transport
    Mesh *mesh;

    /// Material database to provide various material dependent data
    MaterialDatabase *material_database;

    /// one global time marks table
    TimeMarks * main_time_marks;

    /// steady or unsteady water flow simulator based on MH scheme
    DarcyFlowMH *water;
    /// output object for water flow
    DarcyFlowMHOutput *water_output;

    /// explicit transport with chemistry through operator splitting
    TransportBase *transport_reaction;

};

#endif /* HC_EXPLICIT_SEQUENTIAL_HH_ */
