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
 * @file time_constraints.cc
 * @brief List of possible time constraints.
 *  @author Pavel Exner
 */

#include "system/system.hh"

#include "tools/time_constraints.hh"

using namespace std;

/**************************   implemenation of TimeConstraints   *****************************/

TimeConstraintList::TimeConstraintList()
{}


bool TimeConstraintList::define_constraint(TimeConstraintName name, string message, double value)
{
    TimeConstraintIt tc;
    if(find(name,tc)) return false;   //already defined
    
    time_constraints_.push_back( TimeConstraint (name, message, value) );
    return true;
}


bool TimeConstraintList::find(TimeConstraintName name, TimeConstraintIt &time_constraint_it)
{
    time_constraint_it = time_constraints_.begin();
    for(; time_constraint_it != time_constraints_.end(); ++time_constraint_it)
    {   
        if(time_constraint_it->name() == name) 
        {
            return true;
        }
    }
    
    // not found
    return false;
}


TimeConstraint& TimeConstraintList::get(TimeConstraintName name)
{
    TimeConstraintIt tc_it;
    if(find(name,tc_it))
        return *tc_it;
    else
        THROW(ExcConstraintUndefined() << EI_Undefined(name));
    
    //not going to happen
    return *tc_it;
}


void TimeConstraintList::print_all(ostream& stream)
{
    for(TimeConstraint &tc: time_constraints_)
    {
        stream << tc.to_string() << std::endl;
    }
}



/**************************   implemenation of TimeConstraint   *****************************/

TimeConstraint::TimeConstraint()
: name_(""), message_(""), value_(0)
{}

TimeConstraint::TimeConstraint(TimeConstraintName name, string message, double value)
: name_(name), message_(message), value_(value)
{}

string TimeConstraint::to_string() const
{
    std::stringstream s;
    s << name_ << " [" << value_ << "]: " << message_;
    return s.str();
}
