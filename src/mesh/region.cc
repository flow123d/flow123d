/*
 * material_dispatch.cc
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include <string>
#include <sstream>

#include "mesh/region.hh"
#include "system/exceptions.hh"

#include "input/input_type.hh"
#include "input/type_base.hh"
#include "input/accessors.hh"
#include <boost/foreach.hpp>

using namespace std;



std::string Region::label() const
    { return db_->get_label(idx_); }



unsigned int Region::id() const
    { return db_->get_id(idx_); }



unsigned int Region::dim() const
    { return db_->get_dim(idx_); }

/**************************************************************************************************
 * Implementation of     RegionDB
 */


namespace IT=Input::Type;



IT::Record RegionDB::region_input_type =
        IT::Record("Region", "Definition of region of elements.")
        .declare_key("name",IT::String(), IT::Default::obligatory(),
                "Label (name) of the region. Has to be unique in one mesh.\n")
        .declare_key("id", IT::Integer(0), IT::Default::obligatory(),
                "The ID of the region to which you assign label.")
        .declare_key("element_list", IT::Array( IT::Integer(0) ), IT::Default::optional(),
                "Specification of the region by the list of elements. This is not recomended")
        .close();

IT::Record RegionDB::region_set_input_type =
        IT::Record("RegionSet", "Definition of one region set.")
        .declare_key("name", IT::String(), IT::Default::obligatory(),
                "Unique name of the region set.")
        .declare_key("region_ids", IT::Array( IT::Integer(0)),
                "List of region ID numbers that has to be added to the region set.")
        .declare_key("region_labels", IT::Array( IT::String()),
                "List of labels of the regions that has to be added to the region set.")
        .declare_key("union", IT::Array( IT::String(), 2,2),
                "Defines region set as a union of given pair of sets. Overrides previous keys.")
        .declare_key("intersection", IT::Array( IT::String(), 2,2),
                "Defines region set as an intersection of given pair of sets. Overrides previous keys.")
        .declare_key("difference", IT::Array( IT::String(), 2,2),
                "Defines region set as a difference of given pair of sets. Overrides previous keys.")
        .close();


/// Default constructor
RegionDB::RegionDB()
: closed_(false), n_boundary_(0), n_bulk_(0)  {

    // adding implicit boundary and bulk regions
    // How to deal with dimension, clean solution is to have implicit region for every
    // dimension, or we can allow regions of mixed dimension
    implicit_bulk_ = add_region(Region::undefined-1, "IMPLICIT BULK", 0, Region::bulk);
    implicit_boundary_ = add_region(Region::undefined-2, "IMPLICIT BOUNDARY", 0, Region::boundary);

}


Region RegionDB::add_region( unsigned int id, const std::string &label, unsigned int dim, bool boundary) {
    if (closed_) xprintf(PrgErr, "Can not add to closed region DB.\n");

    Region r_id = find_id(id);
    Region r_label=find_label(label);

    if (r_id.is_valid() && r_label.is_valid()) {
        // both iterators are valid; check they are same, and match new values
        if (r_id.idx() != r_label.idx() ) THROW(ExcInconsistentAdd()
                << EI_Label(label) << EI_ID(id) << EI_LabelOfOtherID(r_id.label()) << EI_IDOfOtherLabel(r_label.id())
                );

        if ( r_id.dim() != dim || r_id.is_boundary() != boundary ) {
            //DBGMSG("dims: %d %d, bc: %d %d\n", r_id.dim(),dim, r_id.is_boundary(), boundary);
            THROW(ExcInconsistentAdd() << EI_Label(label) << EI_ID(id) );
        }
        return r_id;
    } else {
        // check that both finds failed
        if (! r_id.is_valid() && ! r_label.is_valid() ) {
            // if DB is open add new entry
            if (closed_)  THROW( ExcAddingIntoClosed() << EI_Label(label) <<EI_ID(id) );
            else {
                unsigned int index;
                if (boundary) (index = (n_boundary_ <<1)), n_boundary_++;
                else (index = (n_bulk_ << 1)+1),  n_bulk_++;
                if (index >= max_n_regions) xprintf(UsrErr, "Too many regions, more then %d\n", max_n_regions);
                if ( ! region_set_.insert( RegionItem(index, id, label, dim) ).second )
                   THROW( ExcCantAdd()  << EI_Label(label) <<EI_ID(id) );
                return Region(index, *this);
            }
        } else {
            if ( r_id.is_valid() ) THROW(ExcInconsistentAdd()
                    << EI_Label(label) << EI_ID(id) << EI_LabelOfOtherID(r_id.label())
                    );
            else
                if (r_label.is_valid() ) THROW(ExcInconsistentAdd()
                        << EI_Label(label) << EI_ID(id) << EI_IDOfOtherLabel(r_label.id())
                        );
                else xprintf(PrgErr, "Rotten inconsistency.\n");
        }
    }
    return Region(); // should not happen
}



Region RegionDB::add_region(unsigned int id, unsigned int dim) {
    Region r_id = find_id(id);
    if (r_id.is_valid()) return add_region(id, r_id.label(), dim, r_id.is_boundary() );
    // else
    stringstream ss;
    ss << "region_" << id;
    return add_region(id, ss.str(), dim, false);
}



Region RegionDB::find_label(const std::string &label) const
{
    RegionTable::index<Label>::type::iterator it_label = region_set_.get<Label>().find(label);
    if (it_label==region_set_.get<Label>().end()  ) return Region();
    return Region(it_label->index, *this);
}



Region RegionDB::find_id(unsigned int id) const
{
    RegionTable::index<ID>::type::iterator it_id = region_set_.get<ID>().find(id);
    if ( it_id==region_set_.get<ID>().end() ) return Region();
    return Region(it_id->index, *this);
}




const std::string & RegionDB::get_label(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_set_.get<Index>().find(idx);
    ASSERT( it!= region_set_.get<Index>().end(), "No region with index: %u\n", idx);
    return  it->label;
}



unsigned int RegionDB::get_id(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_set_.get<Index>().find(idx);
    ASSERT( it!= region_set_.get<Index>().end(), "No region with index: %u\n", idx);
    return  it->id;
}



unsigned int RegionDB::get_dim(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_set_.get<Index>().find(idx);
    ASSERT( it!= region_set_.get<Index>().end(), "No region with index: %u\n", idx);
    return  it->dim_;
}



void RegionDB::close() {
    closed_=true;
    // set default sets
   for(unsigned int i=0; i< size(); i++) {
       Region reg(i, *this);

       if (reg.is_boundary() && (reg.boundary_idx() < boundary_size()) ) {
           add_to_set("BOUNDARY", reg );
           add_to_set("ALL", reg);
       } else
       if ( (! reg.is_boundary()) && (reg.bulk_idx() < bulk_size()) ) {
           add_to_set("BULK", reg );
           add_to_set("ALL", reg);
       }
    }
}


unsigned int RegionDB::size() const {
    ASSERT(closed_, "RegionDB not closed yet.\n");
    return 2* max(n_boundary_, n_bulk_);
}



unsigned int RegionDB::boundary_size() const {
    ASSERT(closed_, "RegionDB not closed yet.\n");
    return n_boundary_;
}



unsigned int RegionDB::bulk_size() const {
    ASSERT(closed_, "RegionDB not closed yet.\n");
    return n_bulk_;
}




void RegionDB::add_to_set( const string& set_name, Region region) {
	std::map<std::string, RegionSet>::iterator it = sets_.find(set_name);

	if (it == sets_.end()) {
		RegionSet set;
		set.push_back(region);

		sets_.insert( std::make_pair(set_name, set) );
	} else {
		RegionSet & set = (*it).second;
		if ( std::find(set.begin(), set.end(), region)==set.end() ) {
			set.push_back(region); // add region if doesn't exist
		}
	}
}


void RegionDB::add_set( const string& set_name, const RegionSet & set) {
	if (sets_.find(set_name) != sets_.end()) {
		sets_.erase(sets_.find(set_name));
	}

	sets_.insert( std::make_pair(set_name, set) );
}


RegionSet RegionDB::union_sets( const string & set_name_1, const string & set_name_2) {
	RegionSet set_union;
	RegionSet set_1, set_2;
	RegionSet::iterator it;

	prepare_sets(set_name_1, set_name_2, set_1, set_2);
	set_union.resize(set_1.size() + set_2.size());
	it = std::set_union(set_1.begin(), set_1.end(), set_2.begin(), set_2.end(), set_union.begin(), Region::comp);
	set_union.resize(it - set_union.begin());

	return set_union;
}


RegionSet RegionDB::intersection( const string & set_name_1, const string & set_name_2) {
	RegionSet set_insec;
	RegionSet set_1, set_2;
	RegionSet::iterator it;

	prepare_sets(set_name_1, set_name_2, set_1, set_2);
	set_insec.resize(set_1.size() + set_2.size());
	it = std::set_intersection(set_1.begin(), set_1.end(), set_2.begin(), set_2.end(), set_insec.begin(), Region::comp);
	set_insec.resize(it - set_insec.begin());

	return set_insec;
}


RegionSet RegionDB::difference( const string & set_name_1, const string & set_name_2) {
	RegionSet set_diff;
	RegionSet set_1, set_2;
	RegionSet::iterator it;

	prepare_sets(set_name_1, set_name_2, set_1, set_2);
	set_diff.resize(set_1.size() + set_2.size());
	it = std::set_difference(set_1.begin(), set_1.end(), set_2.begin(), set_2.end(), set_diff.begin(), Region::comp);
	set_diff.resize(it - set_diff.begin());

	return set_diff;
}


void RegionDB::prepare_sets( const string & set_name_1, const string & set_name_2,
		RegionSet & set_1, RegionSet & set_2) {
	std::map<std::string, RegionSet>::iterator it_1 = sets_.find(set_name_1);
	std::map<std::string, RegionSet>::iterator it_2 = sets_.find(set_name_2);
	ASSERT( it_1 == sets_.end(), "No region set with name: %s\n", set_name_1.c_str());
	ASSERT( it_2 == sets_.end(), "No region set with name: %s\n", set_name_2.c_str());

	set_1 = (*it_1).second;
	set_2 = (*it_2).second;

	std::stable_sort(set_1.begin(), set_1.end(), Region::comp);
	std::stable_sort(set_2.begin(), set_2.end(), Region::comp);
}


const RegionSet & RegionDB::get_region_set(const string & set_name) const {
	std::map<std::string, RegionSet>::const_iterator it = sets_.find(set_name);
	ASSERT( it != sets_.end(), "No region set with name: %s\n", set_name.c_str());
	return (*it).second;
}


void RegionDB::read_sets_from_input(Input::Record rec) {
	string set_name = rec.val<string>("name");
	Input::Array region_ids = rec.val<Input::Array>("region_ids");
	Input::Array region_labels = rec.val<Input::Array>("region_labels");

	ASSERT( region_ids.size() != region_labels.size(), "Size of arrays region_ids and region_labels must be same\n");

	vector<int> id_values; // ? unsigned int
	vector<string> label_values;
	RegionSet region_set;

	region_ids.copy_to(id_values);
	region_labels.copy_to(label_values);
	for (unsigned int i=0; i<id_values.size(); ++i) {
		Region reg(id_values[i], *this);
		region_set.push_back(reg);
	}

	add_set(set_name, region_set);
}
