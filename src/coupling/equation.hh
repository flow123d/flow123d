/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Abstract base class for equation clasess.
 *
 *  @author Jan Brezina
 */

#ifndef EQUATION_HH_
#define EQUATION_HH_


#include <limits>
#include "time_governor.hh"
#include "time_marks.hh"
#include "input/accessors.hh"

#include <petscvec.h>

class Mesh;
class MaterialDatabase;
namespace Input {
    class Record;
}

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
     * Default constructor. Necessary to make tests fixtures for equations.
     * TODO:
     * Replace setting all in constructor with appropriate getters and setters.
     * Make appropriate checks if key ingredients are initialized.
     */
    EquationBase() : mesh_(NULL), mat_base(NULL), time_(NULL),
    equation_mark_type_(TimeGovernor::marks().new_mark_type()), //creating mark type for new equation
    input_record_()
    {}
    /**
     * Common initialization constructor.
     */
    EquationBase(Mesh &mesh, MaterialDatabase &mat_base, const Input::Record in_rec);

    /**
     * Require virtual destructor also for child classes.
     */
    virtual ~EquationBase() {};

    /**
     *  Child class have to implement computation of solution in actual time.
     */
    virtual void update_solution() {
        // solve equation here ...
        time_->next_time();
    }

    /**
     *  Computation of one time step is split into update_solution() and choose_next_time() in order to allow dependency of the next time step
     *  on other coupled models.
     */
//    virtual void compute_one_step() {
//        update_solution();
//        choose_next_time();
//    }

    /**
     * Fix the next discrete time for computation.
     * Can be rewritten in child class to set possible constrains
     * according to possible equation coefficients or other data which can be result of another model.
     *
     */
    virtual void choose_next_time()
        {time_->fix_dt_until_mark();}

    /**
     * Set external upper time step constrain for time governor of the equation.
     */
    virtual void set_time_upper_constraint(double dt)
        {time_->set_upper_constraint(dt);}
        
    /**
     * Set external lower time step constrain for time governor of the equation.
     */
    virtual void set_time_lower_constraint(double dt)
        {time_->set_lower_constraint(dt);}

    /**
     * Basic getter method returns constant TimeGovernor reference which provides full read access to the time information.
     */
    inline TimeGovernor const &time()
    {
        ASSERT(NONULL(time_),"Time governor was not created.\n");
        return *time_;
    }

    /**
     * Most actual planned time for solution.
     */
    inline double planned_time()
        { return time_->estimate_time(); }

    /**
     * Time of actual solution returned by get_solution_vector().
     */
    inline double solved_time()
        { return time_->t(); }

    /**
     * This getter method provides the computational mesh currently used by the model.
     */
    inline  Mesh &mesh()
    {
        return *mesh_;
    }

    /**
     * This getter method provides the material database of the model.
     * TODO: Maybe it is better to have a database outside and use it to produce input fields.
     */
    inline  MaterialDatabase &material_base()
        {return *mat_base;}

    /**
     * Getter for equation time mark type.
     */
    inline TimeMark::Type mark_type()
        {return equation_mark_type_;}

    /**
     * Child class have to implement getter for sequential solution vector.
     */
    virtual void get_solution_vector(double * &vector, unsigned int &size) =0;

    /**
     * Child class have to implement getter for parallel solution vector.
     */
    virtual void get_parallel_solution_vector(Vec &vector) =0;

protected:
    Mesh * mesh_;
    MaterialDatabase * mat_base;
    TimeGovernor *time_;
    TimeMark::Type equation_mark_type_;
    Input::Record input_record_;
};

/**
 * Demonstration of empty equation class, which can be used if user turns off some equation in the model.
 */
class EquationNothing : public EquationBase {

public:
    EquationNothing(Mesh &mesh, MaterialDatabase &mat_base);

    virtual void get_solution_vector(double * &vector, unsigned int &size) {
        vector = NULL;
        size = 0;
    }

    virtual void get_parallel_solution_vector(Vec &vector) {};

    virtual ~EquationNothing() {};
};


#endif /* EQUATION_HH_ */
