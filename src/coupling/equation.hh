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
class FieldCommonBase;
class Region;
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
     * Default constructor. Necessary to make tests fixtures for equations.
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
        ASSERT( time_,"Time governor was not created.\n");
        return *time_;
    }

    /**
     * Set time governor.
     *
     * Used to set pointer to common time governor of Transport Operator Splitting.
     */
    virtual void set_time(TimeGovernor &time);

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
        {return equation_mark_type_;}

    /**
     * Child class have to implement getter for sequential solution vector.
     */
    virtual void get_solution_vector(double * &vector, unsigned int &size) =0;

    /**
     * Child class have to implement getter for parallel solution vector.
     */
    virtual void get_parallel_solution_vector(Vec &vector) =0;

    /**
     * Child class have to implement setter for local part concentration vector.
     */
    virtual void set_concentration_vector(Vec &vector) =0;

protected:
    Mesh * mesh_;
    TimeGovernor *time_;
    TimeMark::Type equation_mark_type_;
    Input::Record input_record_;
};



/**
 * Base class for a data subclasses of equations. We suggest following structure of an equation class:
 * @code
 *      class SomeEquation : public EqBase {
 *          class EqData : public EqDataBase {
 *          ...
 *          }
 *          ...
 *
 *          ///
 *          EqData data_;
 *      }
 * @endcode
 * The @p SomeEquation::EqData class should be used to introduce all necessary data of the equation in terms of bulk (Field)
 * and boundary (BCField) fields. The base class EqDataBase implements some common operations with these fields as
 * construction of the Input::Type objects, reading fields from the input and calling set_time methods.
 *
 *
 * TODO:
 * Mechanism to access fields in EqData from other equations and vice versa, set (or copy) field from somewhere else.
 * - Is it necessary to have EqDataBase as separate class from EqBase ?
 * - getter mechanism:
 *   template <int spacedim, class T>
 *   Field<spacedim, T> *get_field(string &name) {
 *      FieldCommonBase *base_field = get_field_(name); // EqDataBase should have some name database
 *      if ( typeid(base_field) == typeid(Field<spacedim, T> *)
 *      return static_cast<...>(base_field)
 *   }
 *
 * - setter mechanism:
 *   - Should it suppress field on the input? Not easy to accomplish.
 *     Need way how to modify an Input::Type tree from higher levels,
 *     e.g. Let Equation1 to make its Input::Type::Record, make its deep_copy & modify it and use it by Equation_2
 *     After all Input::Type tree is just a definition of input structure. Seems that we made it too rigid.
 */
class EqDataBase {
public:
	TYPEDEF_ERR_INFO( EI_Domain, string);
	DECLARE_INPUT_EXCEPTION(ExcUnknownDomain,
			<< "Unknown field domain: " << EI_Domain::val << "\n");

    /**
     * The only constructor. The name of the equation has to be provided by parameter @p eq_name.
     */
    EqDataBase(const std::string& eq_name);

    /**
     * Adds given field into list of fields for group operations on fields.
     * Parameters are: @p field pointer, @p name of the key in the input, @p desc - description of the key, and optional parameter
     * @p d_val with default value. This method is rather called through the macro ADD_FIELD
     */
    void add_field( FieldCommonBase *field, const string &name, const string &desc, Input::Type::Default d_val= Input::Type::Default::optional() );

    /**
     * This method returns a Record for
     * - @p eq_class_name should be name of the particular equation class, the name of bulk data record has form:
     *   'EqName_BulkData' and record for boundary data has name 'EqName_BoundaryData'. However, these names has
     *   only documentation purpose since these records are not descendants of an AbstractRecord.
     */
    Input::Type::Record generic_input_type(bool bc_regions);

    /**
     * Return Input::Type for field descriptor record. The key 'bulk_data' of an equation should be declared as
     * Array of this record.
     */
    virtual Input::Type::Record bulk_input_type();

    /**
     * Return Input::Type for field descriptor record. The key 'boundary_data' of an equation should be declared as
     * Array of this record.
     */
    virtual Input::Type::Record boundary_input_type();

    /**
     * Actualize fields for actual state of the time governor. We assume that this method is called
     * before the actual time in time governor is solved.
     *
     * - call set_time, which should read up to the first bigger time
     * - check that all fields are initialized
     *
     */
    virtual void set_time(const TimeGovernor &time);

    /**
     * Set mesh pointer in EqDataBase so that it can be set in  those fields that needs it.
     * The mesh has to be set before initialization of fields from input.
     */
    void set_mesh(Mesh *mesh);

    /**
     * Set bulk and boundary region data list accessors.
     * The method takes accessor @p bulk_list of an input array containing records with
     * definition of fields on individual regions (Region data list) and similar accessor @p
     * bc_list with Region data list for boundary fields. The method do not initialize fields
     * from the lists this is done in set_time.
     *
     * Either of accessors can be empty. You can get an empty accessor by default constructor Input::Array().
     */
    void init_from_input(Input::Array bulk_list, Input::Array bc_list);

    /**
     * Reads input from one region - one time descriptor.
     */
    virtual RegionSet read_boundary_list_item(Input::Record rec);

    virtual RegionSet read_bulk_list_item(Input::Record rec);

    virtual ~EqDataBase();

protected:
    EqDataBase();

    void check_times(Input::Array &list);

    void set_time(const TimeGovernor &time, Input::Array &list, Input::Iterator<Input::Record> &it, bool bc_regions);

    RegionSet read_list_item(Input::Record rec, bool bc_regions);

    /// Pointer to mesh where the equation data fields live.
    Mesh *mesh_;
    /// Equation name. Used to name input type records.
    std::string equation_name_;
    /// List of all fields.
    std::vector<FieldCommonBase *> field_list;

    /// Accessors to to bulk and boundary input arrays.
    Input::Array bulk_input_array_, boundary_input_array_;
    /// Iterators into these arrays pointing to the first unprocessed item.
    Input::Iterator<Input::Record> bulk_it_, boundary_it_;
};


/**
 * Macro to simplify call of EqDataBase::add_field method. Two forms are supported:
 *
 * ADD_FIELD(some_field, description);
 * ADD_FIELD(some_field, description, Default);
 *
 * The first form adds name "some_field" to the field member some_field, also adds description of the field. No default
 * value is specified, so the user must initialize the field on all regions (This is checked at the end of the method
 * EqDataBase::init_from_input.
 *
 * The second form adds also default value to the field, that is Default(".."), or Default::read_time(), other default value specifications are
 * meaningless. The automatic conversion to FieldConst is used, e.g.  Default::("0.0") is automatically converted to
 * { TYPE="FieldConst", value=[ 0.0 ] } for a vector valued field, so you get zero vector on output on regions with default value.
 */
#define ADD_FIELD(name, ...)                   add_field(&name, string(#name), __VA_ARGS__)






/**
 * Demonstration of empty equation class, which can be used if user turns off some equation in the model.
 */
class EquationNothing : public EquationBase {

public:
    EquationNothing(Mesh &mesh);

    virtual void get_solution_vector(double * &vector, unsigned int &size) {
        vector = NULL;
        size = 0;
    }

    virtual void get_parallel_solution_vector(Vec &vector) {};

    virtual ~EquationNothing() {};
};


#endif /* EQUATION_HH_ */
