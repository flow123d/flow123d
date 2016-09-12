/*
 * region_set.cc
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include "mesh/region_set.hh"
#include "mesh/mesh.h"
#include <string>
#include <sstream>


namespace IT = Input::Type;


/*******************************************************************
 * implementation of RegionSetBase
 */

RegionSetBase::RegionSetBase(Mesh *mesh)
: region_db_(mesh->region_db_),
  el_to_reg_map_(mesh->region_db_.el_to_reg_map_) {}

IT::Abstract & RegionSetBase::get_input_type() {
	return IT::Abstract("Region", "Abstract record for Region.")
			.close();
}



/*******************************************************************
 * implementation of RegionSetFromId
 */

RegionSetFromId::RegionSetFromId(const Input::Record &rec, Mesh *mesh)
: RegionSetBase(mesh)
{
	Region reg;
	string region_label = rec.val<string>("name");
	unsigned int region_id = rec.val<unsigned int>("id");

	try {
		reg = mesh->region_db().find_id(region_id);
	} catch(RegionDB::ExcUniqueRegionId &e) {
		e << rec.ei_address();
	}
	if ( reg.is_valid() ) {
		try {
			region_db_.rename_region(reg, region_label);
		} catch (RegionDB::ExcNonuniqueLabel &e) {
			e << rec.ei_address();
			throw;
		}
	} else {
		unsigned int dim = rec.val<unsigned int>("dim", RegionDB::undefined_dim);
		try {
			region_db_.add_region(region_id, region_label, dim, rec.address_string() );
		} catch (RegionDB::ExcNonuniqueLabel &e) {
			e << rec.ei_address();
			throw;
		} catch (RegionDB::ExcAddingIntoClosed &e) {
			e << rec.ei_address();
			throw;
		} catch (RegionDB::ExcCantAdd &e) {
			e << rec.ei_address();
			throw;
		}
	}
}



const IT::Record & RegionSetFromId::get_region_input_type()
{
    return IT::Record("From_Id", "Region declared by id and name.\n"
    					"Allows to create new region with given id and label\n"
    					"or specify existing region by id which will be renamed.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.declare_key("dim", IT::Integer(0), IT::Default::optional(),
				"The dim of the region to which you assign label. "
				"Value is taken into account only if new region is created.")
		.close();
}



const int RegionSetFromId::registrar =
		Input::register_class< RegionSetFromId, const Input::Record &, Mesh * >("From_Id") +
		RegionSetFromId::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetFromLabel
 */

RegionSetFromLabel::RegionSetFromLabel(const Input::Record &rec, Mesh *mesh)
: RegionSetBase(mesh)
{
	string new_name = rec.val<string>("name");
	string mesh_label = rec.val<string>("mesh_label");

	Region reg = mesh->region_db().find_label(mesh_label);
	if ( reg.is_valid() ) {
		try {
			region_db_.rename_region(reg, new_name);
		} catch (RegionDB::ExcNonuniqueLabel &e) {
			e << rec.ei_address();
			throw;
		}
	} else {
		WarningOut().fmt("Unknown region in mesh with label '%s'\n", mesh_label);
	}
}



const IT::Record & RegionSetFromLabel::get_region_input_type()
{
    return IT::Record("From_Label", "Allows to rename existing region specified by mesh_label.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"New label (name) of the region. Has to be unique in one mesh.")
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
: RegionSetBase(mesh)
{
	unsigned int region_id;
	string region_label = rec.val<string>("name");
	Input::Iterator<unsigned int> it = rec.find<unsigned int>("id");

	if (it) {
		region_id = (*it);
	} else {
		region_id = this->get_max_region_id();
	}

	try {
		region_db_.add_region(region_id, region_label, RegionDB::undefined_dim, rec.address_string() );
	} catch (RegionDB::ExcNonuniqueLabel &e) {
		e << rec.ei_address();
		throw;
	} catch (RegionDB::ExcAddingIntoClosed &e) {
		e << rec.ei_address();
		throw;
	} catch (RegionDB::ExcCantAdd &e) {
		e << rec.ei_address();
		throw;
	}

	Input::Array element_list = rec.val<Input::Array>("element_list");
	std::vector<unsigned int> element_ids;
	for (Input::Iterator<unsigned int> it_element = element_list.begin<unsigned int>();
			it_element != element_list.end();
	        ++it_element) {
		std::map<unsigned int, unsigned int>::iterator it_map = el_to_reg_map_.find((*it_element));
		if (it_map != el_to_reg_map_.end()) {
			WarningOut().fmt("Region assigned to element with id {} will be rewritten.\n", (*it_element));
		}
		el_to_reg_map_.insert( std::make_pair((*it_element), region_id) );

	}
}



const IT::Record & RegionSetFromElements::get_region_input_type()
{
    return IT::Record("From_Elements", "Region declared by name, ID and enum of elements.\n"
    								   "Allows to create new region and assign elements to its.\n"
    								   "Elements are specified by ids.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::optional(),
				"The ID of the region to which you assign label.\n"
				"If new region is created and ID is not set, unique ID will be generated automatically.")
		.declare_key("element_list", IT::Array( IT::Integer(0), 1 ), IT::Default::obligatory(),
				"Specification of the region by the list of elements.")
		.close();
}



const int RegionSetFromElements::registrar =
		Input::register_class< RegionSetFromElements, const Input::Record &, Mesh * >("From_Elements") +
		RegionSetFromElements::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetBoundary
 * Need new implementation, will be solved later.
 */

// RegionSetBoundary::RegionSetBoundary(const Input::Record &rec, Mesh *mesh) : RegionSetBase(mesh) {}

// const IT::Record & RegionSetBoundary::get_region_input_type() {}

//const int RegionSetBoundary::registrar =
//		Input::register_class< RegionSetBoundary, const Input::Record &, Mesh * >("Region_Boundary") +
//		RegionSetBoundary::get_region_input_type().size();



/*******************************************************************
 * implementation of RegionSetUnion
 */

RegionSetUnion::RegionSetUnion(const Input::Record &rec, Mesh *mesh)
: RegionSetBase(mesh)
{
	string name_of_set = rec.val<string>("name");
	Input::Iterator<Input::Array> region_ids = rec.find<Input::Array>("region_ids");
	Input::Iterator<Input::Array> regions = rec.find<Input::Array>("regions");
	std::vector<string> set_names;

	if (regions) {
		set_names = region_db_.get_and_check_operands(*regions);
	}
	if (region_ids) {
		for (Input::Iterator<unsigned int> it_ids = region_ids->begin<unsigned int>();
				it_ids != region_ids->end();
		        ++it_ids) {
			Region reg;
			try {
				reg = region_db_.find_id(*it_ids);
			} catch(RegionDB::ExcUniqueRegionId &e) {
				e << region_ids->ei_address();
				throw;
			}
			if (!reg.is_valid()) {
				std::string label = region_db_.create_label_from_id(*it_ids);
				reg = region_db_.add_region(*it_ids, label, RegionDB::undefined_dim, rec.address_string() );
			}
			set_names.push_back( reg.label() );
		}

	}

	RegionSet region_set = region_db_.union_set(set_names);
	if (region_set.size() == 0) {
		THROW( ExcEmptyRegionSetResult() << EI_Operation_Type("Union") << rec.ei_address() );
	}
	region_db_.add_set(name_of_set, region_set);
}



const IT::Record & RegionSetUnion::get_region_input_type()
{
    return IT::Record("Union", "Defines region as a union of given two or more regions.\n"
    						   "Regions can be given by names or IDs or both ways together.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("region_ids", IT::Array( IT::Integer(0)),
				"List of region ID numbers that has to be added to the region set.")
		.declare_key("regions", IT::Array( IT::String() ),
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
: RegionSetBase(mesh)
{
	string name_of_set = rec.val<string>("name");
	Input::Iterator<Input::Array> labels = rec.find<Input::Array>("regions");

	std::vector<string> set_names = mesh->region_db().get_and_check_operands(*labels);
	OLD_ASSERT( set_names.size() == 2, "Wrong number of operands. Expect 2.\n" );

	RegionSet set_1 = region_db_.get_region_set( set_names[0] );
	RegionSet set_2 = region_db_.get_region_set( set_names[1] );
	RegionSet set_diff;

	std::stable_sort(set_1.begin(), set_1.end(), Region::comp);
	std::stable_sort(set_2.begin(), set_2.end(), Region::comp);
	set_diff.resize(set_1.size() + set_2.size());

	RegionSet::iterator it = std::set_difference(set_1.begin(), set_1.end(), set_2.begin(), set_2.end(), set_diff.begin(), Region::comp);
	set_diff.resize(it - set_diff.begin());

	if (set_diff.size() == 0) {
		THROW( ExcEmptyRegionSetResult() << EI_Operation_Type("Difference") << rec.ei_address() );
	}
	region_db_.add_set(name_of_set, set_diff);
}



const IT::Record & RegionSetDifference::get_region_input_type()
{
    return IT::Record("Difference", "Defines region as a difference of given pair of regions.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("regions", IT::Array( IT::String(), 2, 2 ), IT::Default::obligatory(),
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
: RegionSetBase(mesh)
{
	string name_of_set = rec.val<string>("name");
	Input::Iterator<Input::Array> regions = rec.find<Input::Array>("regions");
	std::vector<string> set_names = mesh->region_db().get_and_check_operands(*regions);

	RegionSet region_set = region_db_.get_region_set(set_names[0]);
	for (unsigned int i=1; i<set_names.size(); i++) {
		region_set = this->intersection( region_set, set_names[i] );
	}

	if (region_set.size() == 0) {
		THROW( ExcEmptyRegionSetResult() << EI_Operation_Type("Intersection") << rec.ei_address() );
	}
	region_db_.add_set(name_of_set, region_set);
}



const IT::Record & RegionSetIntersection::get_region_input_type()
{
    return IT::Record("Intersection", "Defines region as an intersection of given two or more regions.")
        .derive_from(RegionSetBase::get_input_type())
		.declare_key("name", IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("regions", IT::Array( IT::String(), 2 ), IT::Default::obligatory(),
				"Defines region as an intersection of given pair of regions.")
		.close();
}


const int RegionSetIntersection::registrar =
		Input::register_class< RegionSetIntersection, const Input::Record &, Mesh * >("Intersection") +
		RegionSetIntersection::get_region_input_type().size();


RegionSet RegionSetIntersection::intersection( RegionSet target_set, const string & source_set_name) const {
	RegionSet set_insec;
	RegionSet source_set = region_db_.get_region_set( source_set_name );
	RegionSet::iterator it;

	std::stable_sort(target_set.begin(), target_set.end(), Region::comp);
	std::stable_sort(source_set.begin(), source_set.end(), Region::comp);

	set_insec.resize(target_set.size() + source_set.size());
	it = std::set_intersection(target_set.begin(), target_set.end(), source_set.begin(), source_set.end(), set_insec.begin(), Region::comp);
	set_insec.resize(it - set_insec.begin());

	return set_insec;
}
