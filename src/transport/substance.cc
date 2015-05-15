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


#include <iostream>
#include <iomanip>

#include "transport/substance.hh"



using namespace Input::Type;


Record & Substance::get_input_type() {
	static Record type = Record("Substance", "Chemical substance.")
		.declare_key("name", String(), Default::obligatory(), "Name of the substance.")
		.declare_key("molar_mass", Double(0), Default("1"), "Molar mass of the substance [kg/mol].")
		.allow_auto_conversion("name")
		.close();

	return type;
}




Substance::Substance()
	: name_(""),
	  molar_mass_(1)
{
}

Substance::Substance(const Input::Record &in_rec)
{
	name_ = in_rec.val<std::string>("name");
	molar_mass_ = in_rec.val<double>("molar_mass");
}







void SubstanceList::initialize(const Input::Array &in_array)
{
	substances_ = boost::make_shared<std::vector<Substance> >();
	names_ = boost::make_shared<std::vector<std::string> >();

	for (auto it = in_array.begin<Input::Record>(); it != in_array.end(); ++it)
	{
		Substance s(*it);
		(*substances_).push_back(s);
		(*names_).push_back(s.name());
	}
}

void SubstanceList::initialize(SubstanceList &list)
{
	substances_ = list.substances_;
	names_ = list.names_;
}


void SubstanceList::initialize(const std::vector<std::string> &names)
{
	substances_ = boost::make_shared<std::vector<Substance> >();
	names_ = boost::make_shared<std::vector<std::string> >(names);

	// copy names to internal vectors
	(*substances_).resize(names.size());
	for (unsigned int i=0; i<names.size(); ++i) (*substances_)[i].name_ = names[i];
	*names_ = names;
}

