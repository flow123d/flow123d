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

#include <petscmat.h>
#include "time_governor.hh"


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

