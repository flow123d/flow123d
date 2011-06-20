/*
 * Equation.hh
 *
 *  Created on: May 18, 2011
 *      Author: jb
 */

#ifndef EQUATION_HH_
#define EQUATION_HH_


#include <petscmat.h>
#include <time_governor.hh>
class Mesh;
class MaterialDatabase;
class TimeGovernor;



/**
 * Class EquationBase is abstract base class for a general time dependent model. This class should provide general interface
 * that can be used for general coupling of various particular models. By a model we mean a discrete solver of
 * an partial or ordinary differential equation. Result of the model at one discrete time level should be a discrete field class (not yet implemented).
 * Until we have field classes we only provide method @fn get_solution_vector, which returns pointer to sequential C array with linear combination of
 * base functions that represents the solution.
 *
 * This class does not implement any constructor. In particular it does not initialize mesh, mat_base, and time. This has to be done in the constructor
 * of particular child class.
 *
 * TODO: clarify initialization of data members
 *
 */
class EquationBase {
public:
    /**
     *  Child class have to implement computation of one time step.
     */
    virtual void compute_one_step() =0;

    /**
     * This method implements basic cycle for computation until a given time. But could be overwritten at child class.
     */
    virtual void compute_until( double end_time)
    {
        ASSERT(NONULL(time),"Time governor was not created.\n");
        while ( ! time->is_end() ) compute_one_step();
    }

    /**
     * Basic getter method returns constant TimeGovernor reference which provides full read access to the time information.
     */
    inline const TimeGovernor& get_time()
        {return *time;}

    /**
     * This getter method provides the computational mesh currently used by the model.
     */
    inline  Mesh *get_mesh()
        {return mesh;}

    /**
     * This getter method provides the material database of the model.
     * TODO: Maybe it is better to have a database outside and use it to produce input fields.
     */
    inline  MaterialDatabase *get_mat_base()
        {return mat_base;}

    /**
     * Child class have to implement getter for sequential solution vector.
     */
    virtual void get_solution_vector(double * &vector, unsigned int &size) =0;

    /**
     * Child class have to implement getter for parallel solution vector.
     */
    virtual void get_parallel_solution_vector(Vec &vector) =0;

protected:
    Mesh   *mesh;
    MaterialDatabase *mat_base;
    TimeMarks   *time_marks;
    TimeGovernor *time;
};


#endif /* EQUATION_HH_ */
