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
 * @brief Classes for storing substance data.
 *  @author Jan Stebel
 */


#ifndef SUBSTANCE_HH_
#define SUBSTANCE_HH_




#include "input/accessors.hh"
#include "system/named_vector.hh"

/**
 * Class Substance is a storage for data which are specific for a (chemical) substance.
 * The purpose is to easily share them among equations (e.g. between transport and reactions).
 */
class Substance {

public:

	/// Default constructor.
	Substance();

	/// Initialization from input tree.
	Substance(const Input::Record &in_rec);

	/// Getter for substance name.
	std::string name() const { return name_; }

	/// Getter for molar mass.
	double molar_mass() const { return molar_mass_; }


	/// Input type for a substance.
	static Input::Type::Record input_type;

protected:

	/// Name of a chemical substance.
	std::string name_;

	/// Molar mass [kg/mol] of the substance.
	double molar_mass_;

	friend class SubstanceList;
};


/**
 * SubstanceList is an envelope around a vector of substances, which provides
 * some additional functionality such as:
 * - various ways of initialization (fron JSON input, reference, or list of names)
 * - export of vector of names (required e.g. by some field classes)
 */
class SubstanceList {

public:

	/// Read from input array.
	void initialize(const Input::Array &in_array);

	/// Bind to existing list.
	void initialize(SubstanceList &list);

	/// Construct from a list of names.
	void initialize(const std::vector<std::string> &names);

	inline const Substance &operator[](unsigned int index) { return (*substances_)[index]; }

	inline const std::vector<std::string> &names() { return (*names_); }

	unsigned int size() const { return substances_->size(); }

private:

	/// The actual list of substances.
	boost::shared_ptr<std::vector<Substance> > substances_;

	/// Auxiliary list of substance names used in some classes.
	boost::shared_ptr<std::vector<std::string> > names_;
};


class SubstanceArray : public NamedVector<Substance>
{
public:
    /// Read from input array.
    void initialize(const Input::Array &in_array);
};





#endif // SUBSTANCE_HH_
