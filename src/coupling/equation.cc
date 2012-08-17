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
#include "time_marks.hh"

#include "equation.hh"
#include "system/system.hh"
#include "input/accessors.hh"


EquationBase::EquationBase(TimeMarks &marks, Mesh &mesh, MaterialDatabase &mat_base,const  Input::Record in_rec)
: time_marks(&marks),
  mesh_(&mesh),
  mat_base(&mat_base),
  time_(NULL),
  equation_mark_type_(time_marks->new_mark_type()),
  input_record_(in_rec)
{}



EquationNothing::EquationNothing(TimeMarks &marks, Mesh &mesh, MaterialDatabase &mat_base)
: EquationBase(marks, mesh, mat_base, Input::Record() )
{}
