/*
 * region_set.cc
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include "mesh/region_set.hh"


namespace IT = Input::Type;


/*******************************************************************
 * implementation of RegionSetBase
 */

IT::Abstract & RegionSetBase::get_input_type() {
	return IT::Abstract("Region", "Abstract record for Region.")
			.close();
}


Region RegionSetBase::add_region(Mesh *mesh, unsigned int id, const std::string &label) {
	return mesh->region_db_.add_region(id, label);
}


Region RegionSetBase::add_region(Mesh *mesh, unsigned int id, const std::string &label, unsigned int dim) {
	return mesh->region_db_.add_region(id, label, dim);
}


Region RegionSetBase::add_region(Mesh *mesh, unsigned int id, unsigned int dim) {
	return mesh->region_db_.add_region(id, dim);
}


void RegionSetBase::add_set(Mesh *mesh, const string& set_name, const RegionSet & set) {
	mesh->region_db_.add_set(set_name, set);
}



/*******************************************************************
 * implementation of RegionSetFromId
 */

RegionSetFromId::RegionSetFromId(const Input::Record &rec, Mesh *mesh)
{
	string region_label = rec.val<string>("name");
	unsigned int region_id = rec.val<unsigned int>("id");
	this->add_region(mesh, region_id, region_label);
}



const IT::Record & RegionSetFromId::get_region_input_type()
{
    return IT::Record("From_Id", "Region declared by id and name.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.close();
}



const int RegionSetFromId::registrar =
		Input::register_class< RegionSetFromId, const Input::Record &, Mesh * >("From_Id") +
		RegionSetFromId::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetFromLabel
 */

RegionSetFromLabel::RegionSetFromLabel(const Input::Record &rec, Mesh *mesh)
{
	string region_name;
	string mesh_label = rec.val<string>("mesh_label");
	if ( !rec.opt_val<string>("name", region_name) ) {
		region_name = mesh_label;
	}

	Region reg = mesh->region_db().find_label(mesh_label);
	if (reg == Region()) {
		xprintf(Warn, "Unknown region in mesh with label '%s'\n", mesh_label.c_str());
	} else {
		unsigned int region_id = reg.id();
		this->add_region(mesh, region_id, region_name);
	}
}



const IT::Record & RegionSetFromLabel::get_region_input_type()
{
    return IT::Record("From_Label", "Region declared by mesh_label and name.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("mesh_label", IT::String(), IT::Default::obligatory(),
				"The mesh_label is e.g. physical volume name in GMSH format.")
		.close();
}



const int RegionSetFromLabel::registrar =
		Input::register_class< RegionSetFromLabel, const Input::Record &, Mesh * >("From_Label") +
		RegionSetFromLabel::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetFromElements
 */

RegionSetFromElements::RegionSetFromElements(const Input::Record &rec, Mesh *mesh)
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



const IT::Record & RegionSetFromElements::get_region_input_type()
{
    return IT::Record("From_Elements", "Region declared by name and enum of elements.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::optional(),
				"The ID of the region to which you assign label.")
		.declare_key("element_list", IT::Array( IT::Integer(0) ), IT::Default::optional(),
				"Specification of the region by the list of elements. This is not recomended")
		.close();
}



const int RegionSetFromElements::registrar =
		Input::register_class< RegionSetFromElements, const Input::Record &, Mesh * >("From_Elements") +
		RegionSetFromElements::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetBoundary
 * Need new implementation, will be solved later.
 */

// RegionSetBoundary::RegionSetBoundary(const Input::Record &rec, Mesh *mesh) {}

// const IT::Record & RegionSetBoundary::get_region_input_type() {}

//const int RegionSetBoundary::registrar =
//		Input::register_class< RegionSetBoundary, const Input::Record &, Mesh * >("Region_Boundary") +
//		RegionSetBoundary::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetUnion
 */

RegionSetUnion::RegionSetUnion(const Input::Record &rec, Mesh *mesh)
{
	string set_name = rec.val<string>("name");
	Input::Iterator<Input::Array> labels = rec.find<Input::Array>("region_labels");

	pair<string,string> set_names = mesh->region_db().get_and_check_operands(*labels);
	RegionSet region_set = mesh->region_db().union_sets( set_names.first, set_names.second );
	add_set(mesh, set_name, region_set);
}



const IT::Record & RegionSetUnion::get_region_input_type()
{
    return IT::Record("Union", "Defines region as a union of given two or more regions.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("region_labels", IT::Array( IT::String(), 2,2 ), IT::Default::obligatory(),
				"Defines region as a union of given pair of regions.")
		.close();
}



const int RegionSetUnion::registrar =
		Input::register_class< RegionSetUnion, const Input::Record &, Mesh * >("Union") +
		RegionSetUnion::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetDifference
 */

RegionSetDifference::RegionSetDifference(const Input::Record &rec, Mesh *mesh)
{
	string set_name = rec.val<string>("name");
	Input::Iterator<Input::Array> labels = rec.find<Input::Array>("region_labels");

	pair<string,string> set_names = mesh->region_db().get_and_check_operands(*labels);
	RegionSet region_set = mesh->region_db().difference( set_names.first, set_names.second );
	add_set(mesh, set_name, region_set);
}



const IT::Record & RegionSetDifference::get_region_input_type()
{
    return IT::Record("Difference", "Defines region as a difference of given pair of regions.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("region_labels", IT::Array( IT::String(), 2, 2 ), IT::Default::obligatory(),
				"Defines region as a difference of given pair of regions.")
		.close();
}



const int RegionSetDifference::registrar =
		Input::register_class< RegionSetDifference, const Input::Record &, Mesh * >("Difference") +
		RegionSetDifference::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetIntersection
 */

RegionSetIntersection::RegionSetIntersection(const Input::Record &rec, Mesh *mesh)
{
	string set_name = rec.val<string>("name");
	Input::Iterator<Input::Array> labels = rec.find<Input::Array>("region_labels");

	pair<string,string> set_names = mesh->region_db().get_and_check_operands(*labels);
	RegionSet region_set = mesh->region_db().intersection( set_names.first, set_names.second );
	add_set(mesh, set_name, region_set);
}



const IT::Record & RegionSetIntersection::get_region_input_type()
{
    return IT::Record("Intersection", "Defines region as an intersection of given two or more regions.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("region_labels", IT::Array( IT::String(), 2,2 ), IT::Default::obligatory(),
				"Defines region as an intersection of given pair of regions.")
		.close();
}


const int RegionSetIntersection::registrar =
		Input::register_class< RegionSetIntersection, const Input::Record &, Mesh * >("Intersection") +
		RegionSetIntersection::get_region_input_type().size();

