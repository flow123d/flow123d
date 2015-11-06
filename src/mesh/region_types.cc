/*
 * region_types.cc
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include "mesh/region_types.hh"


namespace IT = Input::Type;


/*******************************************************************
 * implementation of RegionBase
 */

IT::AbstractRecord & RegionBase::get_input_type() {
	return IT::AbstractRecord("Region", "Abstract record for Region.")
			.close();
}


Region RegionBase::add_region(Mesh *mesh, unsigned int id, const std::string &label) {
	return mesh->region_db_.add_region(id, label);
}


Region RegionBase::add_region(Mesh *mesh, unsigned int id, const std::string &label, unsigned int dim) {
	return mesh->region_db_.add_region(id, label, dim);
}


Region RegionBase::add_region(Mesh *mesh, unsigned int id, unsigned int dim) {
	return mesh->region_db_.add_region(id, dim);
}



/*******************************************************************
 * implementation of RegionFromId
 */

RegionFromId::RegionFromId(const Input::Record &rec, Mesh *mesh)
{
	string region_label = rec.val<string>("name");
	unsigned int region_id = rec.val<unsigned int>("id");
	this->add_region(mesh, region_id, region_label);
}



const IT::Record & RegionFromId::get_region_input_type()
{
    return IT::Record("From_Id", "Region declared by id and name.")
        .derive_from(RegionBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.close();
}



/*******************************************************************
 * implementation of RegionFromLabel
 */

RegionFromLabel::RegionFromLabel(const Input::Record &rec, Mesh *mesh)
{
	string mesh_label = rec.val<string>("mesh_label");
	string region_name = rec.val<string>("name", mesh_label);

	unsigned int region_id = mesh->region_db().find_label(mesh_label).id();
	this->add_region(mesh, region_id, region_name);
}



const IT::Record & RegionFromLabel::get_region_input_type()
{
    return IT::Record("From_Label", "Region declared by mesh_label and name.")
        .derive_from(RegionBase::get_input_type())
		.declare_key("name", IT::String(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("mesh_label", IT::String(), IT::Default::obligatory(),
				"The mesh_label is e.g. physical volume name in GMSH format.")
		.close();
}



/*******************************************************************
 * implementation of RegionFromElements
 */

RegionFromElements::RegionFromElements(const Input::Record &rec, Mesh *mesh)
{
	unsigned int region_id;
	string region_label = rec.val<string>("name");

	Input::Iterator<unsigned int> it = rec.find<unsigned int>("id");
	if (it) {
		region_id = (*it);
	} else {
		Region reg = mesh->region_db().find_label(region_label);
		if ( reg.is_valid() ) region_id = reg.id();
		else THROW( ExcNonexistingLabel() << EI_Region_Label(region_label) );
	}

	this->add_region(mesh, region_id, region_label);

	// TODO We must use MapElementIDToRegionID taken from mesh->region_db or gets as constructor argument
	RegionDB::MapElementIDToRegionID map;
    Input::Array element_list;
	if (rec.opt_val("element_list", element_list) ) {
		for (Input::Iterator<unsigned int> it_element = element_list.begin<unsigned int>();
				it_element != element_list.end();
		        ++it_element) {

			std::map<unsigned int, unsigned int>::iterator it_map = map.find((*it_element));
			if (it_map == map.end()) {
				map.insert( std::make_pair((*it_element), region_id) );
			} else {
				xprintf(Warn, "Element with id %u can't be added more than once.\n", (*it_element));
			}
		}
	}
}



const IT::Record & RegionFromElements::get_region_input_type()
{
    return IT::Record("From_Elements", "Region declared by name and enum of elements.")
        .derive_from(RegionBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::optional(),
				"The ID of the region to which you assign label.")
		.declare_key("element_list", IT::Array( IT::Integer(0) ), IT::Default::optional(),
				"Specification of the region by the list of elements. This is not recomended")
		.close();
}



/*******************************************************************
 * implementation of RegionBoundary
 * Need new implementation, will be solved later.
 */

// RegionBoundary::RegionBoundary(const Input::Record &rec, Mesh *mesh) {}

// const IT::Record & RegionBoundary::get_region_input_type() {}



/*******************************************************************
 * implementation of RegionUnion
 */

RegionUnion::RegionUnion(const Input::Record &rec, Mesh *mesh) {}



const IT::Record & RegionUnion::get_region_input_type()
{
    return IT::Record("Union", "Defines region as a union of given two or more regions.")
        .derive_from(RegionBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.declare_key("region_labels", IT::Array( IT::String() ),
				"List of labels of the regions that has to be added to the region.")
		.declare_key("region_ids", IT::Array( IT::Integer() ),
				"List of IDs of the regions that has to be added to the region.")
		.close();
}



/*******************************************************************
 * implementation of RegionDifference
 */

RegionDifference::RegionDifference(const Input::Record &rec, Mesh *mesh) {}



const IT::Record & RegionDifference::get_region_input_type()
{
    return IT::Record("Difference", "Defines region as a difference of given pair of regions.")
        .derive_from(RegionBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.declare_key("region_labels", IT::Array( IT::String(), 2, 2 ),
				"Defines given pair of regions.")
		.close();
}



/*******************************************************************
 * implementation of RegionIntersection
 */

RegionIntersection::RegionIntersection(const Input::Record &rec, Mesh *mesh) {}



const IT::Record & RegionIntersection::get_region_input_type()
{
    return IT::Record("Intersection", "Defines region as an intersection of given two or more regions.")
        .derive_from(RegionBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.declare_key("region_labels", IT::Array( IT::String() ),
				"List of labels of the regions that define intersection of the region.")
		.declare_key("region_ids", IT::Array( IT::Integer() ),
				"List of IDs of the regions that define intersection of the region.")
		.close();
}

