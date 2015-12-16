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
 * @file    region.cc
 * @brief   
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



const IT::Record & RegionDB::get_region_input_type() {
	return IT::Record("Region", "Definition of region of elements.")
		.declare_key("name",IT::String(), IT::Default::obligatory(),
				"Label (name) of the region. Has to be unique in one mesh.\n")
		.declare_key("id", IT::Integer(0), IT::Default::obligatory(),
				"The ID of the region to which you assign label.")
		.declare_key("element_list", IT::Array( IT::Integer(0) ), IT::Default::optional(),
				"Specification of the region by the list of elements. This is not recomended")
		.close();
}

const IT::Record & RegionDB::get_region_set_input_type() {
	return IT::Record("RegionSet", "Definition of one region set.")
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
}

const unsigned int RegionDB::undefined_dim = 10;


/// Default constructor
RegionDB::RegionDB()
: closed_(false), n_boundary_(0), n_bulk_(0)  {

    // adding implicit boundary and bulk regions
    // How to deal with dimension, clean solution is to have implicit region for every
    // dimension, or we can allow regions of mixed dimension
    //implicit_bulk_ = add_region(Region::undefined-1, "IMPLICIT BULK", 0, Region::bulk);
    //implicit_boundary_ = ;

}


Region RegionDB::implicit_boundary_region() {
	DimIDIter it_id = region_table_.get<DimId>().find(DimID(undefined_dim, Region::undefined-2));
    if ( it_id!=region_table_.get<DimId>().end() ) {
        return Region(it_id->index, *this);
    }

    return insert_region(Region::undefined-2, "IMPLICIT BOUNDARY", undefined_dim, Region::boundary);
}


Region RegionDB::add_region( unsigned int id, const std::string &label, unsigned int dim) {
	bool boundary = is_boundary(label);
    DimIDIter it_id = region_table_.get<DimId>().find(DimID(dim,id));
    if (it_id != region_table_.get<DimId>().end() ) {
    	// Find existing region
    	return find_by_dimid(it_id, id, label, boundary);
    }

    DimIDIter it_undef_dim = region_table_.get<DimId>().find(DimID(undefined_dim,id));
	if (it_undef_dim != region_table_.get<DimId>().end() ) {
		// Region with same ID and undefined_dim exists, replace undefined_dim
		return replace_region_dim(it_undef_dim, dim, boundary);
    }

    LabelIter it_label = region_table_.get<Label>().find(label);
    if (it_label != region_table_.get<Label>().end() ) {
        // ID is free, not label
        THROW(ExcNonuniqueLabel() << EI_Label(label) << EI_ID(id) << EI_IDOfOtherLabel(it_label->get_id()) );
    }

    return insert_region(id, label, dim, boundary);
}


Region RegionDB::add_region(unsigned int id, const std::string &label) {
	bool boundary = is_boundary(label);
	if (label.size() == 0) create_label_from_id(label, id);

    OnlyIDIter it_only_id = region_table_.get<OnlyID>().find(id);
    if (it_only_id != region_table_.get<OnlyID>().end() ) {
    	// replace region label
    	unsigned int index = it_only_id->index;

    	RegionItem item(index, it_only_id->get_id(), label, it_only_id->dim());
    	region_table_.replace(
    			region_table_.get<Index>().find(index),
                item);

    	return Region(index, *this);
    }

    LabelIter it_label = region_table_.get<Label>().find(label);
    if (it_label != region_table_.get<Label>().end() ) {
        // ID is free, not label
        THROW(ExcNonuniqueLabel() << EI_Label(label) << EI_ID(id) << EI_IDOfOtherLabel(it_label->get_id()) );
    }

    return insert_region(id, label, undefined_dim, boundary);
}


Region RegionDB::add_region(unsigned int id, unsigned int dim) {
	DimIDIter it_id = region_table_.get<DimId>().find(DimID(dim,id));
    if ( it_id!=region_table_.get<DimId>().end() ) {
        return Region(it_id->index, *this);
    }

    DimIDIter it_undef_dim = region_table_.get<DimId>().find(DimID(undefined_dim,id));
	if (it_undef_dim != region_table_.get<DimId>().end() ) {
		// Region with same ID and undefined_dim exists, replace undefined_dim
		bool boundary = is_boundary(it_undef_dim->label);
		return replace_region_dim(it_undef_dim, dim, boundary);
    }

    // else
    stringstream ss;
    ss << "region_" << id;
    return insert_region(id, ss.str(), dim, false);
}



Region RegionDB::find_label(const std::string &label) const
{
    LabelIter it_label = region_table_.get<Label>().find(label);
    if (it_label==region_table_.get<Label>().end()  ) return Region();
    return Region(it_label->index, *this);
}





Region RegionDB::find_id(unsigned int id, unsigned int dim) const
{
	DimIDIter it_id = region_table_.get<DimId>().find(DimID(dim, id));
    if ( it_id==region_table_.get<DimId>().end() ) return Region();
    return Region(it_id->index, *this);
}




Region RegionDB::find_id(unsigned int id) const
{
	if (region_table_.get<OnlyID>().count(id) > 1) {
		THROW( ExcUniqueRegionId() << EI_ID( id ) );
	}
	OnlyIDIter it_id = region_table_.get<OnlyID>().find(id);
	if ( it_id==region_table_.get<OnlyID>().end() ) return Region();
	return Region(it_id->index, *this);
}




const std::string & RegionDB::get_label(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_table_.get<Index>().find(idx);
    ASSERT( it!= region_table_.get<Index>().end(), "No region with index: %u\n", idx);
    return  it->label;
}



unsigned int RegionDB::get_id(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_table_.get<Index>().find(idx);
    ASSERT( it!= region_table_.get<Index>().end(), "No region with index: %u\n", idx);
    return  it->get_id();
}



unsigned int RegionDB::get_dim(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_table_.get<Index>().find(idx);
    ASSERT( it!= region_table_.get<Index>().end(), "No region with index: %u\n", idx);
    return  it->dim();
}



void RegionDB::close() {
    closed_=true;
    // set default sets
   for(unsigned int i=0; i< size(); i++) {
       Region reg(i, *this);
   	   RegionSet region_set;
   	   region_set.push_back( reg );


       if (reg.is_boundary() && (reg.boundary_idx() < boundary_size()) ) {
           add_to_set("BOUNDARY", reg );
           add_to_set("ALL", reg);
       	   add_set(reg.label(), region_set);
       } else
       if ( (! reg.is_boundary()) && (reg.bulk_idx() < bulk_size()) ) {
           add_to_set("BULK", reg );
           add_to_set("ALL", reg);
       	   add_set(reg.label(), region_set);
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
    // add region only if it is not in the set
    if (sets_.find(set_name) == sets_.end()) {
	    sets_.insert( std::make_pair(set_name, set) );
	}
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

	if ( it_1 == sets_.end() ) { THROW(ExcUnknownSet() << EI_Label(set_name_1)); }
	if ( it_2 == sets_.end() ) { THROW(ExcUnknownSet() << EI_Label(set_name_2)); }

	set_1 = (*it_1).second;
	set_2 = (*it_2).second;

	std::stable_sort(set_1.begin(), set_1.end(), Region::comp);
	std::stable_sort(set_2.begin(), set_2.end(), Region::comp);
}



pair<string,string> RegionDB::get_and_check_operands(const Input::Array & operands)
{
	vector<string> names;
	operands.copy_to(names);
	if ( names.size() != 2 ) THROW(ExcWrongOpNumber() << EI_NumOp(names.size()) << operands.ei_address() );
	auto ret_names = pair<string,string>(names[0], names[1]);
	if ( sets_.find( ret_names.first ) == sets_.end() )
		THROW( ExcUnknownSet()  << EI_Label( ret_names.first )
								<< operands.ei_address() );
	if ( sets_.find( ret_names.second ) == sets_.end() )
		THROW( ExcUnknownSet()  << EI_Label( ret_names.second )
								<< operands.ei_address() );
	return ret_names;
}



RegionSet RegionDB::get_region_set(const string & set_name) const {
	std::map<std::string, RegionSet>::const_iterator it = sets_.find(set_name);
	if ( it == sets_.end() ) {
		return RegionSet();
	}
	return (*it).second;
}


void RegionDB::read_sets_from_input(Input::Array arr) {
	for (Input::Iterator<Input::Record> it = arr.begin<Input::Record>();
			it != arr.end();
			++it) {

		Input::Record rec = (*it);
		string set_name;
		Input::Array region_ids, region_labels;
		Input::Array union_names, intersection_names, difference_names;
		RegionSet region_set;

		rec.opt_val("name", set_name);

		if (rec.opt_val("region_ids", region_ids) ) {
			for (Input::Iterator<unsigned int> it_ids = region_ids.begin<unsigned int>();
					it_ids != region_ids.end();
			        ++it_ids) {
				try {
					Region reg = find_id(*it_ids);
					if (reg.is_valid()) {
						if ( std::find(region_set.begin(), region_set.end(), reg)==region_set.end() ) {
							region_set.push_back(reg); // add region if doesn't exist
						}
					} else {
						xprintf(Warn, "Region with id %d doesn't exist. Skipping\n", (*it_ids));
					}
				} catch(ExcUniqueRegionId &e) {
					e << region_ids.ei_address();
					throw;
				}
			}
		}

		if (rec.opt_val("region_labels", region_labels) ) {
			for (Input::Iterator<string> it_labels = region_labels.begin<string>();
					it_labels != region_labels.end();
			        ++it_labels) {
				Region reg = find_label(*it_labels);
				if (reg.is_valid()) {
					if ( std::find(region_set.begin(), region_set.end(), reg)==region_set.end() ) {
						region_set.push_back(reg); // add region if doesn't exist
					}
				} else {
					xprintf(Warn, "Region with label %s doesn't exist. Skipping\n", (*it_labels).c_str());
				}
			}
		}

		Input::Iterator<Input::Array> operands = rec.find<Input::Array>("union");
		if ( operands ) {

			if (region_set.size() != 0) {
				xprintf(Warn, "Overwriting previous initialization of region set '%s' by union operation.\n", set_name.c_str());
			}

			pair<string,string> set_names = get_and_check_operands(*operands);
			region_set = union_sets( set_names.first, set_names.second );
		}

		operands = rec.find<Input::Array>("intersection");
		if (operands) {

			if (region_set.size() != 0) {
				xprintf(Warn, "Overwriting previous initialization of region set '%s' by intersection operation.\n", set_name.c_str());
			}

			pair<string,string> set_names = get_and_check_operands(*operands);
			region_set = intersection( set_names.first, set_names.second );
		}

		operands = rec.find<Input::Array>("difference");
		if (operands) {

			if (region_set.size() != 0) {
				xprintf(Warn, "Overwriting previous initialization of region set '%s' by difference operation.\n", set_name.c_str());
			}

			pair<string,string> set_names = get_and_check_operands(*operands);
			region_set = difference( set_names.first, set_names.second );
		}

		add_set(set_name, region_set);
	}
}

void RegionDB::read_regions_from_input(Input::Array region_list, MapElementIDToRegionID &map) {
	map.clear();

	for (Input::Iterator<Input::Record> it = region_list.begin<Input::Record>();
				it != region_list.end();
				++it) {

		Input::Record rec = (*it);
		string region_name = rec.val<string>("name");
		unsigned int region_id = rec.val<unsigned int>("id");
		add_region(region_id, region_name);

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
}

void RegionDB::create_label_from_id(const string & label, unsigned int id) {
	stringstream ss;
	ss << "region_" << id;
	ss.str(label);
}

Region RegionDB::insert_region(unsigned int id, const std::string &label, unsigned int dim, bool boundary) {
	if (closed_) THROW( ExcAddingIntoClosed() << EI_Label(label) << EI_ID(id) );

	unsigned int index;
	if (boundary) {
		index = (n_boundary_ <<1);
		n_boundary_++;
	} else  {
		index = (n_bulk_ << 1)+1;
		n_bulk_++;
	}
	if (index >= max_n_regions) xprintf(UsrErr, "Too many regions, more then %d\n", max_n_regions);
	if ( ! region_table_.insert( RegionItem(index, id, label, dim) ).second )
	   THROW( ExcCantAdd() << EI_Label(label) << EI_ID(id) );
	return Region(index, *this);
}

Region RegionDB::replace_region_dim(DimIDIter it_undef_dim, unsigned int dim, bool boundary) {
	ASSERT( it_undef_dim->dim() == undefined_dim,
			"Dimension of replaced region with id=%u must be undefined_dim, actually is: %u\n", it_undef_dim->get_id(), it_undef_dim->dim());

	unsigned int index = it_undef_dim->index;

	RegionItem item(index, it_undef_dim->get_id(), it_undef_dim->label, dim);
	region_table_.replace(
			region_table_.get<Index>().find(index),
            item);

	Region r_id=Region(index, *this);
    // check boundary
    if ( r_id.is_boundary() != boundary )
        THROW(ExcInconsistentBoundary() << EI_Label(it_undef_dim->label) << EI_ID(it_undef_dim->get_id()) );

    return r_id;
}

Region RegionDB::find_by_dimid(DimIDIter it_id, unsigned int id, const std::string &label, bool boundary) {
    unsigned int index = it_id->index;
    LabelIter it_label = region_table_.get<Label>().find(label);
    if ( it_label == region_table_.get<Label>().end() || index != it_label->index )
    	THROW(ExcNonuniqueID() << EI_Label(label) << EI_ID(id) << EI_LabelOfOtherID(it_id->label) );

    Region r_id=Region(index, *this);
    // check boundary
    if ( r_id.is_boundary() != boundary )
        THROW(ExcInconsistentBoundary() << EI_Label(label) << EI_ID(id) );

    return r_id;
}
