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
#include "system/system.hh"


VtkMeshReader::VtkMeshReader(const FilePath &file_name)
{
	parse_result_ = doc_.load_file( ((std::string)file_name).c_str() );
	read_base_vtk_attributes();
}



VtkMeshReader::VtkMeshReader(std::istream &in)
{
	parse_result_ = doc_.load(in);
	read_base_vtk_attributes();
}



void VtkMeshReader::read_base_vtk_attributes()
{
	pugi::xml_node node = doc_.child("VTKFile");
	// flag of compressed data
	std::string compressor = node.attribute("compressor").as_string();
	compressed_ = (compressor == "vtkZLibDataCompressor");
	// header type of appended data
	header_type_ = node.attribute("header_type").as_string();
	// size of node and element vectors
	node = node.child("UnstructuredGrid").child("Piece");
	n_nodes_ = node.attribute("NumberOfPoints").as_uint();
	n_elements_ = node.attribute("NumberOfCells").as_uint();
}



VtkMeshReader::DataArrayAttributes VtkMeshReader::get_data_array_attr(DataSections data_section, std::string data_array_name)
{
    static std::vector<std::string> section_names = {"Points", "Cells", "CellData"};

    pugi::xml_node node = doc_.child("VTKFile").child("UnstructuredGrid").child("Piece");
    if (data_section == DataSections::points) {
    	node = node.child("Points").child("DataArray");
    } else {
    	ASSERT(data_array_name!="").error("Empty Name attribute of DataArray!\n");
    	node = node.child( section_names[data_section].c_str() ).find_child_by_attribute("DataArray", "Name", data_array_name.c_str());
    }

    DataArrayAttributes attributes;
    attributes.name_ = node.attribute("name").as_string();
    attributes.type_ = this->get_data_type( node.attribute("type").as_string() );
    attributes.n_components_ = node.attribute("NumberOfComponents").as_uint(1);
    std::string format = node.attribute("format").as_string();
    if (format=="appended") {
        attributes.format_ = DataFormat::appended;
    } else if (format=="ascii") {
        attributes.format_ = DataFormat::ascii;
    } else {
        ASSERT(false).error("Unsupported or missing VTK format.");
    }
    attributes.offset_ = node.attribute("offset").as_uint();
    return attributes;
}



VtkMeshReader::DataType VtkMeshReader::get_data_type(std::string type_str, bool only_integral) {
    if (type_str=="UInt32") {
        return DataType::uint32;
    } else if (type_str=="UInt64") {
    	return DataType::uint64;
    } else if (type_str=="Float64" && !only_integral) {
    	return DataType::float64;
    } else {
        ASSERT(false).error("Unsupported or missing VTK data type.");
        return DataType::uint32;
    }

}
