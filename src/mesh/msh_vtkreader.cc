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
 * @file    msh_vtkreader.cc
 * @brief
 * @author  dalibor
 */


#include <iostream>
#include "boost/lexical_cast.hpp"
#include "msh_vtkreader.hh"


VtkMeshReader::VtkMeshReader(const FilePath &file_name)
{
	parse_result_ = doc_.load_file( ((std::string)file_name).c_str() );
	read_nodes_elms_count();
}



VtkMeshReader::VtkMeshReader(std::istream &in)
{
	parse_result_ = doc_.load(in);
	read_nodes_elms_count();
}



void VtkMeshReader::read_nodes_elms_count()
{
	pugi::xml_node piece_node = doc_.child("VTKFile").child("UnstructuredGrid").child("Piece");
	n_nodes_ = boost::lexical_cast<unsigned int>( piece_node.attribute("NumberOfPoints").value() );
	n_elements_ = boost::lexical_cast<unsigned int>( piece_node.attribute("NumberOfCells").value() );
}

