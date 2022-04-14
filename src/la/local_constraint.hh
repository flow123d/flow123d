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
 * @file    local_constraint.hh
 * @brief
 */


#ifndef LOCAL_CONSTRAINT_HH_
#define LOCAL_CONSTRAINT_HH_


#include "system/asserts.hh"
class LocSystem;

enum class ConstraintType {
    none =0b00,
    rows =0b01,
    cols =0b10,
	both =0b11
};


class LocalConstraint
{
public:
	/// Default constructor
	LocalConstraint() : type_(ConstraintType::none) {}

	/// Constructor. Set all data of local constraint.
	LocalConstraint(uint i_elm, uint loc_dof, double solution, double diag=0.0)
	: type_(ConstraintType::both), i_element_(i_elm), loc_dof_(loc_dof), solution_(solution), diag_(diag) {}

	/// Set only row data to constraint.
	void set_row(uint i_elm, uint loc_dof, double solution, double diag=0.0)
	{
	    ASSERT(type_ == ConstraintType::none);
		type_ = ConstraintType::rows;
		i_element_ = i_elm;
		loc_dof_ = loc_dof;
		solution_ = solution;
		diag_ = diag;
	}

	/// Set only column data to constraint.
	void set_col(uint i_elm, uint loc_dof, double solution)
	{
	    ASSERT(type_ == ConstraintType::none);
		type_ = ConstraintType::cols;
		i_element_ = i_elm;
		loc_dof_ = loc_dof;
		solution_ = solution;
		diag_ = 0.0;
	}

	/// Getter to i_element
	inline uint i_element() const
	{
	    return this->i_element_;
	}

    bool operator < (const LocalConstraint &other) {
        return (i_element_ < other.i_element_);
    }

private:
	ConstraintType type_;
	uint i_element_;
	uint loc_dof_;
	double solution_;
	double diag_;

	friend class LocalSystem;
};

#endif // LOCAL_CONSTRAINT_HH_
