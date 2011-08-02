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
#include <limits>

#include "system/system.hh"

class Mesh;
class MaterialDatabase;
class TimeGovernor;



/**
 * Class EquationBase is abstract base class for a general time dependent model. This class should provide general interface
 * that can be used for general coupling of various particular models. By a model we mean a discrete solver of
 * an partial or ordinary differential equation. Result of the model at one discrete time level should be a discrete field class (not yet implemented).
 * Until we have field classes we only provide method get_solution_vector(), which returns pointer to sequential C array with linear combination of
 * base functions that represents the solution.
 *
 * Computation of one time step (method compute_one_step() )  is split into update_solution() and choose_next_time().
 *
 * This class does not implement any constructor. In particular it does not initialize mesh, mat_base, and time. This has to be done in the constructor
 * of particular child class.
 *
 * Any constructor of child class should set solved = true. We assume, that after initialization an equation object stay solve in init time. For the first time step
 * one calls method chose_next_time() which setup time frame of the first time step.
 *
 * TODO: clarify initialization of data members
 *
 */
class EquationBase {
public:
    /**
     *  Child class have to implement computation of solution in actual time.
     */
    virtual void update_solution() {
        solved=true;
    }

    /**
     *  Computation of one time step is split into update_solution() and choose_next_time() in order to allow dependency of the next time step
     *  on other coupled models.
     */
    virtual void compute_one_step() {
        update_solution();
        choose_next_time();
    }

    /**
     * Choose the next discrete time for computation. However do this only if the last planned time has been solved.
     * Can be rewritten in child class to set possible constrains
     * according to possible equation coefficients or other data which can be result of another model.
     *
     *
     */
    virtual void choose_next_time() {
        DBGMSG("solved: %d %d\n",solved, is_end() );
        if (is_end()) return;
        if (solved) {

            time_->next_time();
            solved =false;
        }
    }

    /**
     * This method implements basic cycle for computation until a given time. But could be overwritten at child class.
     */
    virtual void compute_until( double end_time)
    {
        ASSERT(NONULL(time_),"Time governor was not created.\n");
        while ( ! time_->is_end() ) compute_one_step();
    }

    /**
     * Basic getter method returns constant TimeGovernor reference which provides full read access to the time information.
     */
    inline const TimeGovernor &time()
    {
        ASSERT(NONULL(time_),"Time governor was not created.\n");
        return *time_;
    }

    /**
     * Most actual planned time for solution.
     */
    inline double planned_time()
        { return time_->t(); }

    /**
     * Time of actual solution returned by get_solution_vector().
     */
    inline double solved_time()
        { return solved ? time_->t() : time_->last_t(); }

    /**
     * Returns true if planned_time is solved_time.
     */
    inline bool is_solved()
        {return solved;}

    /**
     * Returns true if solved_time is the end point of the time interval of the time governor.
     */
    inline bool is_end()
        {
        DBGMSG("eq end: %f %d\n", time_->t(), solved);
        return time_->is_end() && solved;
        }

    /**
     * This getter method provides the computational mesh currently used by the model.
     */
    inline  Mesh &mesh()
    {
        ASSERT(NONULL(mesh_),"Mesh was not created.\n");
        return *mesh_;
    }

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
    bool solved;

    Mesh   *mesh_;
    MaterialDatabase *mat_base;
    TimeMarks   *time_marks;
    TimeGovernor *time_;
};

/**
 * Demonstration of empty equation class, which can be used if user turns off some equation in the model.
 */
class EquationNothing : public EquationBase {

public:
    EquationNothing() {}

    virtual void get_solution_vector(double * &vector, unsigned int &size) {
        vector = NULL;
        size = 0;
    }

    virtual void get_parallel_solution_vector(Vec &vector) {};

    ~EquationNothing() {};
};


#endif /* EQUATION_HH_ */
