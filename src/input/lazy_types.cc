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


void LazyTypes::addType(const boost::shared_ptr<TypeBase> &type)
{
	types.push_back(type);
}


void LazyTypes::finish()
{
    // TODO: dynamic cast as the switch may be expensive, in such case use some notification about type

	// first finish all lazy input types save Selection (we have to leave open Selection in AbstractType key TYPE)
	for (TypeVector::iterator it=types.begin(); it!=types.end(); it++) {
	    if (dynamic_pointer_cast<Selection>(*it) == 0) {
	        (*it)->finish();
	    }
	}

	// then finalize abstract records so that no type can derive from them
	for (TypeVector::iterator it=types.begin(); it!=types.end(); it++)
	{
	    boost::shared_ptr<AbstractRecord> a_rec_ptr = dynamic_pointer_cast<AbstractRecord>(*it);
		if ( a_rec_ptr!= 0) a_rec_ptr->no_more_descendants();
	}

	// at last finish all selections (including those in AbstractRecord)
	for (TypeVector::iterator it=types.begin(); it!=types.end(); it++) (*it)->finish();

	types.clear();
}

