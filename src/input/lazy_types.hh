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

#ifndef LAZY_TYPES_HH_
#define LAZY_TYPES_HH_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>


namespace Input {

namespace Type {


class TypeBase;
class Selection;

/**
 * The class LazyType is an abstract base for input types that cannot be initialized at once but in two steps.
 * This concerns in particular Record, AbstractRecord, Array and Selection. These input types are typically defined by means
 * of static variables, whose order of initialization is not known a priori. Since e.g. a Record can link to other
 * input types through its keys, these input types cannot be accessed directly at the initialization phase.
 * The remaining part of initialization can be done later, typically from main(), by calling the method finish().
 */
/*
class LazyType
{
public:

	virtual void finish() = 0;

    /// Empty virtual destructor.
    virtual ~LazyType( void ) {}
};
*/

/**
 * The Singleton class LazyTypes serves for handling the lazy-evaluated input types, derived from the base class
 * LazyType. When all static variables are initialized, the method LazyTypes::instance().finish() can be called
 * in order to finish initialization of lazy types such as Records, AbstractRecords, Arrays and Selections.
 * Selections have to be finished after all other types since they are used by AbstractRecords to register all
 * derived types. For this reason LazyTypes contains two arrays - one for Selections, one for the rest.
 */
class LazyTypes
{
public:

	typedef std::vector< boost::shared_ptr<TypeBase> > TypeVector;
	typedef std::vector< boost::shared_ptr<TypeBase> > SelectionVector;

	/**
	 * The reference to the singleton instance.
	 */
	static LazyTypes &instance()
	{
		static LazyTypes INSTANCE;
		return INSTANCE;
	}


	/**
	 * Registers new lazy type to be finished later.
	 */
	void addType(const boost::shared_ptr<TypeBase> &type);


	/// Finishes all registered lazy types.
	void finish();

private:
	/// Private default constructor ensures that no other instances are created.
	LazyTypes();

	/// The array of registered lazy types.
	TypeVector types;

};


} // closing namespace Type
} // closing namespace Input

#endif /* LAZY_TYPES_HH_ */
