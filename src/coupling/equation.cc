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
 * @file    equation.cc
 * @brief   Abstract base class for equation clasess.
 * @author  Jan Brezina
 */

#include <petscmat.h>
#include "tools/time_governor.hh"


#include "equation.hh"
#include "system/system.hh"
#include "input/accessors.hh"
#include "fields/field_set.hh"

#include <boost/foreach.hpp>



/*****************************************************************************************
 * Implementation of EqBase
 */

EquationBase::EquationBase()
: equation_empty_(true),
  mesh_(NULL),
  time_(NULL),
  input_record_(),
  eq_data_(nullptr)
{}



EquationBase::EquationBase(Mesh &mesh, const  Input::Record in_rec)
: equation_empty_(false),
  mesh_(&mesh),
  time_(NULL),
  input_record_(in_rec),
  eq_data_(nullptr)
{}


void EquationBase::set_time_governor(TimeGovernor &time)
{
  time_ = &time;
}

