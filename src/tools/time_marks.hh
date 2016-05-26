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
 * @file    time_marks.hh
 * @brief   
 * @author  Jan Brezina
 *          Created on: Jun 15, 2011
 */

#ifndef TIME_MARKS_HH_
#define TIME_MARKS_HH_

#include <ostream>
#include <vector>

#include "system/system.hh"
#include "system/global_defs.h"


/**
 * @brief Class used for marking specified times at which some events occur.
 * 
 * This class represents one record in the TimeMarks simple database.
 * Class members can not be modified after the item is created.
 */
class TimeMark {
public:

    /**
     *  MarkType is a bitmap where each bit represents one base type such as (Output, Input, BC change, Fixed...)
     *  This allow more complex queries through bitwise operations. Also one TimeMark can be shared by more events.
     *  In the context of TimeMarks the Type can be either fixed or vague. If a TimeGovernor is connected to the TimeMarks object
     *  the  TimeMarks with fixed Type are used to match exactly their times. Base Types should be created/obtained from TimeMarks class
     *  through the TimeMarks::new_mark_type method.
     *  There are three Types predefined in TimeMarks constructor:
     *  - type_fixed_time (hex 0x01)
     *  - type_output (hex 0x02)
     *  - type_bc_change (hex 0x04)
     *  @see TimeMarks
     */
    typedef unsigned long int Type;

    /// Mark Type with all bits set.
    static const Type every_type;
    /// Mark Type with all bits unset.
    static const Type none_type;
    
    /**
     * Constructor for a TimeMarks::Mark
     * @param time time of the mark
     * @param type type of the mark
     *
     * In order to create a fixed TimeMark (at time=0.1) with base TimeMark::Type my_type, use the TimeMarks class:
     * TimeMark( 0.1, timemarks.type_fixed_time() | my_type)
     */
    TimeMark(double time, Type type) :
        time_(time), mark_type_(type) {}


    /// Getter for mark type.
    inline Type mark_type() const {
        return mark_type_;
    }

    /// Getter for the time of the TimeMark.
    inline double time() const {
        return time_;
    }

    /**
     * Returns true if TimeMark's type has 1 on all positions where mask has 1.
     * @param mask {Select bits that should be 1 for matching mark types.
     */

    inline bool match_mask(const TimeMark::Type &mask) const {
        return ( mask & (~mark_type_) ) == 0;
    }

    /// Add more bits that a mark satisfies.
    /// @param type type that should be modified
    inline void add_to_type(const TimeMark::Type &type) {
        mark_type_ |= type;
    }

    /// Comparison of time marks according to their time.
    /// @param another is another Timemark which should be compared.
    bool operator<(const TimeMark& another) const
      { return time_ < another.time(); }


private:
    /// The marked time.
    double time_;
    /// The type of the TimeMark.
    Type mark_type_;
};

/**
 * Output to stream operator for TimeMark class.
 */
std::ostream& operator<<(std::ostream& stream, const TimeMark &marks);




/***************************************************************************************/





/***************************************************************************************/
class TimeGovernor;
class TimeMarksIterator;

/**
 * @brief This class is a collection of time marks to manage various events occurring during simulation time.
 *
 * <b> TimeMark and their types </b>
 * 
 * One TimeMark consists of time and type (TimeMark::Type), see the constructor TimeMark::TimeMark.
 * The type of mark is bitmap where individual bits corresponds to some base event types like changing a BC, output solution, coupling time with another
 * equation and so on. Base types can be combined by bitwise or (operator|).
 *
 * Special types are predefined in TimeMarks class. These are returned by their getters:
 * - type_fixed_time() - marks of this type are considered as fixed times by a TimeGovernor which is connected to particular TimeMarks object.
 * - type_output() - this type marks times at which solution should be computed and written to output.
 * - type_bc_change() - this type marks times at which BC is is changed and model has to be updated. 
 * 
 * <b> TimeMarks collection </b>
 * 
 * TimeMarks collect marks of various types and provides methods for iterating over stored marks. You can selectively access only marks matching given
 * type mask. See TimeMark::match_mask.
 *
 * You can add one new mark through method add or add evenly spaced marks of same type by TimeMarks::add_time_marks.
 *
 * You can allocate new TimeMark::Type in the context of one TimeMarks object by TimeMarks::new_mark_type.
 *
 * For a given TimeGovernor (not necessarily connected one) you can ask about existence of mark in current time interval (TimeMarks::is_current) and see TimeMarks
 * close to the current time (TimeMarks::next and TimeMarks::last). The current time interval is left-open and right-closed: (t,t+dt]. Repeatedly used TimeMarks::next always returns the same TimeMark if the time of the TimeGovernor is not changed.
 *
 * In most cases there will be only one TimeMarks object for the whole solved problem and used by TimeGovernors of individual equations. However
 * this is not necessary.
 * 
 * @see TimeMark
 */
class TimeMarks {

public:
    /// Iterator class for iteration over time marks of particular type. This is always const_iterator.
    typedef TimeMarksIterator iterator;

    /**
     * Default constructor.
     */
    TimeMarks();

    /**
     * Reset state after construction (through default constructor).
     * Useful for unit tests.
     */
    void reinit();

    /**
     * Add a new base mark within the context of the particular TimeMarks instance.
     * User should keep the returned value (MarkType is basically a bitmap) for further queries and
     * TimeMark insertions. ATTENTION: You can not use the TimeMark::Type with other TimeMarks instance!
     * Types are added as they are prepared in next_mark_type_. 
     * Next mark type is updated by (left) bit shifting operator.
     */
    TimeMark::Type new_mark_type();

    /// Predefined base TimeMark type that is taken into account by the TimeGovernor.
    /// Is defined by constructor as 0x01.
    inline TimeMark::Type type_fixed_time()
    { return type_fixed_time_;}

    /// Predefined base TimeMark type for output times.
    /// Is defined by constructor as 0x02.
    inline TimeMark::Type type_output()
    { return type_output_;}

    /// Predefined base TimeMark type for times when the boundary condition is changed.
    /// Is defined by constructor as 0x04.
    inline TimeMark::Type type_input()
    { return type_input_;}


    /**
     * Basic method for inserting TimeMarks.
     * @param mark    Reference to TimeMark object.
     */
    void add(const TimeMark &mark);

    /**
     * Method for creating and inserting equally spaced TimeMarks.
     * @param time    Time of the first TimeMark.
     * @param dt      Lenght of interval between equally spaced TimeMarks.
     * @param end_time  No marks after the end_time.
     * @param type    Type of inserted TimeMarks or their combinations.
     *
     * Current lazy implementation have complexity O(m*n) where m is number of inserted time marks and n number of time marks in the array.
     * TODO: O(n+m) implementation
     */
    void add_time_marks(double time, double dt, double end_time, TimeMark::Type type);

    /**
     * Find the last time mark matching given mask, and returns true if it is in the time interval of
     * current time step.
     */
    bool is_current(const TimeGovernor &tg, const TimeMark::Type &mask) const;

    /**
     * Return the first TimeMark with time strictly greater then tg.time() that match the mask.
     * The time governor tg  is used also for time comparisons.
     *
     * @param tg    the time governor
     * @param mask  mask of marks to iterate on
     *
     * TODO: have also method which accepts double (time) instead of the whole TimeGovernor.
     * and compare without safety.
     */
    TimeMarks::iterator next(const TimeGovernor &tg, const TimeMark::Type &mask) const;

    /**
     * Return the last TimeMark with time less or equal to tg.time() that match the mask.
     * The time governor tg  is used also for time comparisons.
     * @param tg    the time governor
     * @param mask  mask of marks to iterate on
     */
    TimeMarks::iterator last(const TimeGovernor &tg, const TimeMark::Type &mask) const;

    /**
     * Returns iterator to the last time mark matching given @p mask.
     */
    TimeMarks::iterator last(const TimeMark::Type &mask) const;

    /// Iterator for the begin mimics container-like  of TimeMarks
    TimeMarks::iterator begin(TimeMark::Type mask = TimeMark::every_type) const;

    /// Iterator for the end mimics container-like  of TimeMarks
    TimeMarks::iterator end(TimeMark::Type mask = TimeMark::every_type) const;

    /// Friend output operator.
    friend std::ostream& operator<<(std::ostream& stream, const TimeMarks &marks);

private:

    /// MarkType that will be used at next new_time_mark() call.
    TimeMark::Type next_mark_type_;

    /// TimeMarks list sorted according to the their time.
    std::vector<TimeMark> marks_;

    /// Predefined type for fixed time.
    TimeMark::Type type_fixed_time_;
    /// Predefined type for output.
    TimeMark::Type type_output_;
    /// Predefined type for change of boundary condition.
    TimeMark::Type type_input_;
};



/**
 * @brief Iterator over TimeMark objects in TimeMarks object (database of TimeMark objects).
 *
 * Iterator over the TimeMarks of particular mask. This is always const iterator, i.e. it points to const TimeMark.
 * While iterating over TimeMarks with different types, all non matching types are skipped.
 */
class TimeMarksIterator {
public:
    /**Constructor. It is used in TimeMarks class which has the vector of TimeMark objects.
     * @param marks is vector of TimeMark objects.
     * @param it is iterator over the vector of TimeMark objects.
     * @param mask is the type of marks over which we iterate.
     */
    TimeMarksIterator(const std::vector<TimeMark> &marks,const  std::vector<TimeMark>::const_iterator &it, const TimeMark::Type &mask)
    : marks_(marks), it_(it), mask_(mask) {}

    TimeMarksIterator &operator=(const TimeMarksIterator &it)
    {OLD_ASSERT(&marks_ == &it.marks_, "Can not assign TimeMarks::iterator of different container.\n");
     it_=it.it_;
     mask_=it.mask_;
     return *this;
    }

    /// Prefix increment. Skip non matching marks.
    TimeMarksIterator &operator++()
    {
    	while ( it_ != marks_.end() ) {
    		++it_;
    		if (it_->match_mask(mask_)) break;
    	}
    	return (*this);
    }

    /// Prefix decrement. Skip non matching marks.
    TimeMarksIterator &operator--()
    {
    	while ( it_ != marks_.begin() ) {
    		--it_;
    		if (it_->match_mask(mask_)) break;
    	}
    	return (*this);
    }

    ///  * dereference operator
    inline const TimeMark & operator *() const
    {
    	OLD_ASSERT(it_!= marks_.end(), "Out of marks vector.\n");
    	return *it_;
    }

    /// -> dereference operator
    inline const TimeMark * operator ->() const
    {
    	OLD_ASSERT(it_!= marks_.end(), "Out of marks vector.\n");
    	return &(*(it_));
    }

    inline bool operator ==(const TimeMarksIterator &other) const
        {return it_ == other.it_; }

    inline bool operator !=(const TimeMarksIterator &other) const
            {return it_ != other.it_; }

    /// Returns mask.
    TimeMark::Type mask()
    { return mask_; }

private:
    /// Reference to the vector of TimeMark objects.
    const std::vector<TimeMark> &marks_;
    /// Iterator over the vector of TimeMark objects.
    std::vector<TimeMark>::const_iterator it_;
    /// Mask type.
    TimeMark::Type mask_;
};

#endif /* TIME_MARKS_HH_ */
