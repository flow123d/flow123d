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
 * @brief Singleton class that registers all input types that have to be lazy-evaluated.
 *  @author Jan Stebel
 */

#include "input/lazy_types.hh"
#include "input/type_base.hh"
#include "input/type_record.hh"

using namespace Input::Type;


LazyTypes::LazyTypes()
{}

LazyTypes::TypeVector::iterator LazyTypes::begin()
{
	return types.begin();
}

LazyTypes::TypeVector::iterator LazyTypes::end()
{
	return types.end();
}


void LazyTypes::addType(LazyType *type)
{
	types.push_back(type);
}

void LazyTypes::addSelection(LazyType *sel)
{
	selections.push_back(sel);
}


void LazyTypes::finish()
{
	// first finish all lazy input types
	for (TypeVector::iterator it=begin(); it!=end(); it++)
		(*it)->finish();

	// then finalize abstract records so that no type can derive from them
	for (TypeVector::iterator it=begin(); it!=end(); it++)
	{
		if (dynamic_cast<AbstractRecord *>(*it) != 0)
			dynamic_cast<AbstractRecord *>(*it)->no_more_descendants();
	}

	// at last finish all selections
	for (SelectionVector::iterator it=selections.begin(); it!=selections.end(); it++)
		(*it)->finish();

	types.clear();
	selections.clear();
}

