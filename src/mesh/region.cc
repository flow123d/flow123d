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
    //implicit_bulk_ = add_region(Region::undefined-1, "IMPLICIT BULK", 0, Region::bulk);
    //implicit_boundary_ = ;

}


Region RegionDB::implicit_boundary_region() {
    return add_region(Region::undefined-2, "IMPLICIT BOUNDARY", undefined_dim, Region::boundary);
}


void RegionDB::check_dim_consistency(IDIter it_id, unsigned int dim) {
    // check dimension
    if (it_id->dim_ != dim) {
        // User can introduce regions through the mesh input record, however without dim
        // specification. Here we allow overwriting dimension in this case
        if (it_id->dim_ == undefined_dim) {
            RegionItem item(it_id->index, it_id->id, it_id->label, dim);
            region_set_.replace(
                    region_set_.get<Index>().find( it_id->index ),
                    item);
        }
        else THROW(ExcInconsistentDimension() << EI_Label(it_id->label) << EI_ID(it_id->id) );
    }

}


Region RegionDB::add_region( unsigned int id, const std::string &label, unsigned int dim, bool boundary) {
    if (closed_) xprintf(PrgErr, "Can not add to closed region DB.\n");

    IDIter it_id = region_set_.get<ID>().find(id);
    LabelIter it_label = region_set_.get<Label>().find(label);

    if (it_id != region_set_.get<ID>().end() ) {
        unsigned int index = it_id->index;
        if (it_id->dim_ != undefined_dim  && index != it_label->index) THROW(ExcNonuniqueID() << EI_Label(label) << EI_ID(id) << EI_LabelOfOtherID(it_id->label) );

        check_dim_consistency(it_id, dim); // possibly update DB

        Region r_id=Region(index, *this);
        // check boundary
        if ( r_id.is_boundary() != boundary )
            THROW(ExcInconsistentBoundary() << EI_Label(label) << EI_ID(id) );

        return r_id;
    } else
    if (it_label != region_set_.get<Label>().end() ) {
        // ID is free, not label
        THROW(ExcNonuniqueLabel() << EI_Label(label) << EI_ID(id) << EI_IDOfOtherLabel(it_label->id) );
    } else {
        // if DB is open add new entry
        if (closed_)
            THROW( ExcAddingIntoClosed() << EI_Label(label) <<EI_ID(id) );
        else {
            unsigned int index;
            if (boundary) {
                index = (n_boundary_ <<1);
                n_boundary_++;
            } else  {
                index = (n_bulk_ << 1)+1;
                n_bulk_++;
            }
            if (index >= max_n_regions) xprintf(UsrErr, "Too many regions, more then %d\n", max_n_regions);
            if ( ! region_set_.insert( RegionItem(index, id, label, dim) ).second )
               THROW( ExcCantAdd()  << EI_Label(label) <<EI_ID(id) );
            return Region(index, *this);
        }
    }
    return Region(); // should not happen
}


Region RegionDB::add_region(unsigned int id, const std::string &label, unsigned int dim) {
	if (label.size() != 0) {
		bool boundary = label[0] == '.';
		return add_region(id, label, dim, boundary);
	}
    // else
    stringstream ss;
    ss << "region_" << id;
	return add_region(id, ss.str(), dim, false);
}


Region RegionDB::add_region(unsigned int id, unsigned int dim) {
    RegionTable::index<ID>::type::iterator it_id = region_set_.get<ID>().find(id);
    if ( it_id!=region_set_.get<ID>().end() ) {
        // just check dimension
        check_dim_consistency(it_id, dim);
        return Region(it_id->index, *this);
    }
    // else
    stringstream ss;
    ss << "region_" << id;
    return add_region(id, ss.str(), dim, false);
}



Region RegionDB::find_label(const std::string &label) const
{
    LabelIter it_label = region_set_.get<Label>().find(label);
    if (it_label==region_set_.get<Label>().end()  ) return Region();
    return Region(it_label->index, *this);
}





Region RegionDB::find_id(unsigned int id) const
{
    IDIter it_id = region_set_.get<ID>().find(id);
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
	ASSERT( it_1 != sets_.end(), "No region set with name: %s\n", set_name_1.c_str());
	ASSERT( it_2 != sets_.end(), "No region set with name: %s\n", set_name_2.c_str());

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
				Region reg = find_id(*it_ids);
				if (reg.is_valid()) {
					if ( std::find(region_set.begin(), region_set.end(), reg)==region_set.end() ) {
						region_set.push_back(reg); // add region if doesn't exist
					}
				} else {
					xprintf(Err, "Region with id %d doesn't exist.\n", (*it_ids));
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
					xprintf(Err, "Region with label %s doesn't exist.\n", (*it_labels).c_str());
				}
			}
		}

		if (rec.opt_val("union", union_names) ) {

			if (region_set.size() != 0) {
				xprintf(Warn, "Overwriting previous initialization of region set '%s' by union operation.\n", set_name.c_str());
			}

			Input::Iterator<string> union_set_1 = union_names.begin<string>();
			Input::Iterator<string> union_set_2 = ++(union_names.begin<string>());
			region_set = union_sets( (*union_set_1), (*union_set_2) );
		}

		if (rec.opt_val("intersection", intersection_names) ) {

			if (region_set.size() != 0) {
				xprintf(Warn, "Overwriting previous initialization of region set '%s' by intersection operation.\n", set_name.c_str());
			}

			Input::Iterator<string> intersection_set_1 = intersection_names.begin<string>();
			Input::Iterator<string> intersection_set_2 = ++(intersection_names.begin<string>());
			region_set = intersection( (*intersection_set_1), (*intersection_set_2) );
		}

		if (rec.opt_val("difference", difference_names) ) {

			if (region_set.size() != 0) {
				xprintf(Warn, "Overwriting previous initialization of region set '%s' by difference operation.\n", set_name.c_str());
			}

			Input::Iterator<string> difference_set_1 = difference_names.begin<string>();
			Input::Iterator<string> difference_set_2 = ++(difference_names.begin<string>());
			region_set = difference( (*difference_set_1), (*difference_set_2) );
		}

		add_set(set_name, region_set);
	}

	/*
	  Input::Array region_ids;
	  if (rec.opt_val("region_ids", region_ids) ) {
	    ... add regions
	        // for int_item in region_ids -> RegionDB.find_id(int_item)
	        // ... bud prida nalezany region pomoci add_to_set, nebo chyba

	  }

	  if (rec.opt_val("region_labels", region_labels)) {
    // ... podobne pro region_labels, s pouzitim RegionDB::find_label

	  }

	  if (rec.opt_val("union", ...) ) {
	   ...
	   if (region_set.size() != 0) xprintf(Warn, "Overwriting previous initialization of region set 'NAME' by union operation.");
	  }
	 */
}

void RegionDB::read_regions_from_input(Input::Array region_list, MapElementIDToRegionID &map) {
	map.clear();

	for (Input::Iterator<Input::Record> it = region_list.begin<Input::Record>();
				it != region_list.end();
				++it) {

		Input::Record rec = (*it);
		string region_name = rec.val<string>("name");
		unsigned int region_id = rec.val<unsigned int>("id");
		add_region(region_id, region_name, undefined_dim);

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
