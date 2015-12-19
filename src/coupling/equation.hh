/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    equation.hh
 * @brief   Abstract base class for equation clasess.
 * @author  Jan Brezina
 */

#ifndef EQUATION_HH_
#define EQUATION_HH_


#include <limits>
#include "tools/time_governor.hh"
#include "tools/time_marks.hh"
#include "input/accessors.hh"

#include <petscvec.h>

class Mesh;
class Region;
class FieldSet;
typedef std::vector<Region> RegionSet;


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
 * This class does not implement any constructor. In particular it does not initialize mesh and time. This has to be done in the constructor
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
     * Default constructor. Sets all virtual methods empty. Necessary to make tests fixtures for equations.
     * TODO:
     * Replace setting all in constructor with appropriate getters and setters.
     * Make appropriate checks if key ingredients are initialized.
     */
    EquationBase();

    /**
     * Common initialization constructor.
     */
    EquationBase(Mesh &mesh, const Input::Record in_rec);

    /**
     * Require virtual destructor also for child classes.
     */
    virtual ~EquationBase() {};

    /**
     *  Initialization of the solution in the zero time.
     *  There is lot of things that can not be done in the constructor
     *  since we have not fully initialized fields yet. Fields coming from coupling
     *  has to be set after the constructor and before zero_time_step.
     */
    virtual void zero_time_step() {
      if (equation_empty_) DBGMSG("Calling 'zero_time_step' of empty equation '%s'.\n",typeid(*this).name());
      else DBGMSG("Method 'zero_time_step' of '%s' is not implemented.\n",typeid(*this).name());
    }

    /**
     *  Calculation of the next time step and its output.
     */
    virtual void update_solution() {
      if (equation_empty_) DBGMSG("Calling 'update_solution' of empty equation '%s'.\n",typeid(*this).name());
      else DBGMSG("Method 'update_solution' of '%s' is not implemented.\n",typeid(*this).name());
    }

    ///Initialize fields.
    /** 
     * All members that are needed to set fields must be set at this moment (e.g. number of components).
     */
    virtual void initialize() {
      if (equation_empty_) DBGMSG("Calling 'initialize' of empty equation '%s'.\n",typeid(*this).name());
      else DBGMSG("Method 'initialize' of '%s' is not implemented.\n",typeid(*this).name());      
    }
    
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
     * Basic getter method returns TimeGovernor reference which provides full access to the time information.
     */
    inline TimeGovernor &time()
    {
        ASSERT( time_,"Time governor was not created.\n");
        return *time_;
    }

    /**
     * Set time governor.
     *
     * Used to set pointer to common time governor (e.g. in Transport Operator Splitting, Reaction).
     */
    virtual void set_time_governor(TimeGovernor &time);

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
     * Getter for equation time mark type.
     */
    inline TimeMark::Type mark_type()
    {
    	return time().equation_mark_type();
    }

    /**
     * Return reference to the equation data object containing all fields
     * that the equation needs or produce.
     */
    FieldSet &data()
    {
    	ASSERT(eq_data_, "The equation %s did not set eq_data_ pointer.\n", input_record_.address_string().c_str());
    	return *eq_data_;
    }

    /**
     * Child class have to implement getter for sequential solution vector.
     * DEPRECATED
     */
    virtual void get_solution_vector(double * &vector, unsigned int &size)
    { ASSERT(0, "If using, needs to be implemented in ancestors!"); };

    /**
     * Child class have to implement getter for parallel solution vector.
     * DEPRECATED
     */
    virtual void get_parallel_solution_vector(Vec &vector)
    { ASSERT(0, "If using, needs to be implemented in ancestors!"); };

    /**
     * @brief Write computed fields.
     */
    virtual void output_data() {
      if (equation_empty_) DBGMSG("Calling 'output_data' of empty equation '%s'.\n",typeid(*this).name());
      else DBGMSG("Method 'output_data' of '%s' is not implemented.\n",typeid(*this).name());
    }

protected:
    bool equation_empty_;       ///< flag is true if only default constructor was called
    Mesh * mesh_;
    TimeGovernor *time_;
    Input::Record input_record_;
    
    /**
     * Pointer to the equation data object. Every particular equation is responsible
     * to set the pointer in its constructor. This is used by the general method
     * EqData::data(). This approach is simpler than making EqData::data() a virtual method.
     */
    FieldSet *eq_data_;
    
    /// object for calculation and writing the mass balance to file.
    boost::shared_ptr<Balance> balance_;
    
};


#endif /* EQUATION_HH_ */
