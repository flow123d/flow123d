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



const unsigned int RegionDB::undefined_dim = 10;


/// Default constructor
RegionDB::RegionDB()
: closed_(false), n_boundary_(0), n_bulk_(0), max_id_(0)  {

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

    return insert_region(Region::undefined-2, "IMPLICIT BOUNDARY", undefined_dim, Region::boundary, "");
}


Region RegionDB::add_region( unsigned int id, const std::string &label, unsigned int dim, const std::string &address) {
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

    return insert_region(id, label, dim, boundary, address);
}


Region RegionDB::rename_region( Region reg, const std::string &new_label ) {
	ASSERT(reg.is_valid(), "Non-existing region can't be renamed to '%s'.\n", new_label.c_str());

	// test if region with new_label exists
	LabelIter it_label = region_table_.get<Label>().find(new_label);
    if (it_label != region_table_.get<Label>().end() ) {
        if ( reg.id() == it_label->index ) {
            return reg;
        } else {
            // region with same label and other ID exists
            THROW(ExcNonuniqueLabel() << EI_Label(new_label) << EI_ID(reg.id()) << EI_IDOfOtherLabel(it_label->get_id()) );
        }
    }

	// replace region label
	unsigned int index = reg.idx();
	RegionSetTable::iterator it = sets_.find(reg.label()); // rename set
	std::swap(sets_[new_label], it->second);
	sets_.erase(it);
	bool old_boundary_flag = reg.is_boundary(); // check old x new boundary flag - take account in adding to sets

	RegionItem item(index, reg.id(), new_label, reg.dim(), this->get_region_address(index));
	region_table_.replace(
			region_table_.get<Index>().find(index),
            item);

	if (old_boundary_flag != reg.is_boundary()) { // move region between BULK and BOUNDARY sets
		xprintf(Warn, "Change boundary flag of region with id %d and label %s.\n", reg.id(), new_label.c_str());
		if (old_boundary_flag) {
			erase_from_set("BOUNDARY", reg );
			add_to_set("BULK", reg );
		} else {
			erase_from_set("BULK", reg );
			add_to_set("BOUNDARY", reg );
		}
	}

	return Region(index, *this);
}


Region RegionDB::get_region(unsigned int id, unsigned int dim) {
	DimIDIter it_id = region_table_.get<DimId>().find(DimID(dim,id));
    if ( it_id!=region_table_.get<DimId>().end() ) {
        // Region found.
    	return Region(it_id->index, *this);
    }

    DimIDIter it_undef_dim = region_table_.get<DimId>().find(DimID(undefined_dim,id));
	if (it_undef_dim == region_table_.get<DimId>().end() ) {
		// Region doesn't exist.
		return Region();
    } else {
    	// Region with same ID and undefined_dim exists, replace undefined_dim
    	bool boundary = is_boundary(it_undef_dim->label);
    	return replace_region_dim(it_undef_dim, dim, boundary);
    }
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



const std::string & RegionDB::get_region_address(unsigned int idx) const {
    RegionTable::index<Index>::type::iterator it = region_table_.get<Index>().find(idx);
    ASSERT( it!= region_table_.get<Index>().end(), "No region with index: %u\n", idx);
    return it->address;
}



void RegionDB::mark_used_region(unsigned int idx) {
    RegionTable::index<Index>::type::iterator it = region_table_.get<Index>().find(idx);
    ASSERT( it!= region_table_.get<Index>().end(), "No region with index: %u\n", idx);
    if ( !it->used ) {
    	unsigned int index = it->index;
    	RegionItem item(index, it->get_id(), it->label, it->dim(), it->address, true);
    	region_table_.replace(
    			region_table_.get<Index>().find(index),
                item);
    }
}



void RegionDB::close() {
    closed_=true;
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
	RegionSetTable::iterator it = sets_.find(set_name);

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


void RegionDB::erase_from_set( const string& set_name, Region region) {
	RegionSetTable::iterator it = sets_.find(set_name);
	ASSERT(it != sets_.end(), "Region set '%s' doesn't exist.", set_name.c_str());
	RegionSet & set = (*it).second;

	auto set_it = std::find(set.begin(), set.end(), region);
	ASSERT(set_it != set.end(), "Erased region was not found in set '%s'", set_name.c_str());
	set.erase(set_it);
}


std::vector<string> RegionDB::get_and_check_operands(const Input::Array & operands) const
{
	vector<string> names;
	operands.copy_to(names);

	for (string name : names) {
		if ( sets_.find( name ) == sets_.end() )
			THROW( ExcUnknownSet()  << EI_Label( name )
									<< operands.ei_address() );
	}

	return names;
}



RegionSet RegionDB::get_region_set(const string & set_name) const {
	RegionSetTable::const_iterator it = sets_.find(set_name);
	if ( it == sets_.end() ) {
		return RegionSet();
	}
	return (*it).second;
}


string RegionDB::create_label_from_id(unsigned int id) const {
	stringstream ss;
	ss << "region_" << id;
	return ss.str();
}

Region RegionDB::insert_region(unsigned int id, const std::string &label, unsigned int dim, bool boundary, const std::string &address) {
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
	if ( ! region_table_.insert( RegionItem(index, id, label, dim, address) ).second )
	   THROW( ExcCantAdd() << EI_Label(label) << EI_ID(id) );
	if (max_id_ < id) {
		max_id_ = id;
	}

	Region reg = Region(index, *this);
    // add region to sets
	RegionSet region_set;
	region_set.push_back( reg );
	this->add_set(reg.label(), region_set);
	add_to_set("ALL", reg);
    if (reg.is_boundary()) {
        add_to_set("BOUNDARY", reg );
    } else {
        add_to_set("BULK", reg );
    }

    return reg;
}

Region RegionDB::replace_region_dim(DimIDIter it_undef_dim, unsigned int dim, bool boundary) {
	ASSERT( it_undef_dim->dim() == undefined_dim,
			"Dimension of replaced region with id=%u must be undefined_dim, actually is: %u\n", it_undef_dim->get_id(), it_undef_dim->dim());

	unsigned int index = it_undef_dim->index;

	RegionItem item(index, it_undef_dim->get_id(), it_undef_dim->label, dim, this->get_region_address(index));
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

void RegionDB::print_region_table(ostream& stream) const {
	ASSERT(closed_, "RegionDB not closed yet.\n");

	// print header
	stream << endl << "----------- Table of all regions: -----------" << endl;
	stream << std::setfill(' ') << setw(6) << "id" << " dim label" << setw(12) << "" << "contains regions" << endl;
	// print data
	for (RegionSetTable::const_iterator it = sets_.begin(); it != sets_.end(); ++it) { // iterates through RegionSets
		unsigned int reg_id = RegionIdx::undefined;
		string rset_label = it->first;
		LabelIter label_it = region_table_.get<Label>().find(rset_label);
		if ( label_it != region_table_.get<Label>().end() ) { // checks if Region with same label as RegionSet exists
			reg_id = label_it->get_id(); // stores ID of Region - used if RegionSet contains only one Region
			stream << setw(6) << reg_id << setw(4) << label_it->dim();
		} else {
			stream << "     -   -";
		}
		stream << " " << std::left << setw(17) << rset_label << std::right;
		if (it->second.size() > 1) {
			// prints IDs of Regions (if their count is 2 or more)
			for (RegionSet::const_iterator set_it = it->second.begin(); set_it!=it->second.end(); ++set_it) {
				if (set_it != it->second.begin()) stream << ", ";
				stream << set_it->label();
			}
		} else if ( it->second.size() == 1 && (it->second)[0].id() != reg_id ) {
			// for one Region in RegionSet prints ID only if RegionSet is not wrapper of Region
			stream << (it->second)[0].label();
		} else {
			stream << "-";
		}
		stream << endl;
	}
	stream << std::setfill('-') << setw(45) << "" << std::setfill(' ') << endl << endl;
}


void RegionDB::check_regions() {
	ASSERT(closed_, "RegionDB not closed yet.\n");

	for (RegionTable::index<Index>::type::iterator it = region_table_.get<Index>().begin();
			it!= region_table_.get<Index>().end();
			++it) {
		if (!it->used)
			THROW(ExcUnusedRegion() << EI_Label(it->label) << EI_ID(it->get_id()) << Input::EI_Address(it->address) );
	}
}


RegionSet RegionDB::union_set(std::vector<string> set_names) const {
	std::set<Region, bool (*)(const Region&, const Region&)> set(Region::comp);

	for (string set_name : set_names) {
		RegionSet r_set = get_region_set(set_name);
		set.insert(r_set.begin(), r_set.end());
	}

	return RegionSet(set.begin(), set.end());
}
