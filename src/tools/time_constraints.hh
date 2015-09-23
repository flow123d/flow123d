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
 * $Id:  $
 * $Revision:  $
 * $LastChangedBy:  $
 * $LastChangedDate:  $
 *
 * @file time_constraints.hh
 * @brief List of possible time constraints.
 *  @author Pavel Exner
 */

#ifndef TIME_CONSTRAINTS_HH_
#define TIME_CONSTRAINTS_HH_


#include "system/global_defs.h"
#include <iostream>

typedef std::string TimeConstraintName;

class TimeConstraint {
public:
    TimeConstraint();
    TimeConstraint(TimeConstraintName name, std::string message, double value = 0.0);
    TimeConstraintName name() const;
    std::string message() const;
    double value() const;
    void set_value(double new_value);
    std::string to_string() const;
private:
    TimeConstraintName name_;
    std::string message_;
    double value_;
};


class TimeConstraintList{
public:
    TYPEDEF_ERR_INFO( EI_Undefined, std::string);
    DECLARE_EXCEPTION(ExcConstraintUndefined,
            << "Time step constraint '" << EI_Undefined::val << "' has not been defined.");
    
    typedef std::vector<TimeConstraint>::iterator TimeConstraintIt;
    
    TimeConstraintList();
    
    bool define_constraint(TimeConstraintName name, std::string message, double value = 0.0);
    /// Finds the constraint of the given name.
    /** @param name is the name of the constraint
     * @param time_constraint is the output time step constraint
     */
    bool find(TimeConstraintName name, TimeConstraintIt& time_constraint_it);
    TimeConstraint & get(TimeConstraintName name);
    
    /// Sets the time step constraint @p time_constraint with the type @p name and value @p value.
    void set(TimeConstraint& time_constraint, TimeConstraintName name, double value);
    
    void print_all(std::ostream &stream);
    
private:
    std::vector<TimeConstraint> time_constraints_;
};


/**************************   implemenation of TimeConstraints   *****************************/

// inline TimeConstraint& TimeConstraintList::upper()
// {   return *upper_; }
// 
// inline TimeConstraint& TimeConstraintList::lower()
// {   return *lower_; }



/**************************   implemenation of TimeConstraint   *****************************/

inline TimeConstraintName TimeConstraint::name() const
{ return name_; }

inline std::string TimeConstraint::message() const
{ return message_; }

inline double TimeConstraint::value() const
{ return value_;}

inline void TimeConstraint::set_value(double new_value)
{ value_ = new_value;}

#endif //TIME_CONSTRAINTS_HH_