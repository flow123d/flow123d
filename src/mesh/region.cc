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

using namespace std;

RegionDB Region::db_;


std::string Region::label() const
    { return db_.get_label(idx_); }



unsigned int Region::id() const
    { return db_.get_id(idx_); }



unsigned int Region::dim() const
    { return db_.get_dim(idx_); }

/**************************************************************************************************
 * Implementation of     RegionDB
 */


namespace IT=Input::Type;

// adding implicit boundary and bulk regions
// How to deal with dimension, clean solution is to have implicit region for every
// dimension, or we can allow regions of mixed dimension
Region RegionDB::implicit_bulk=Region::db().add_region(Region::undefined-1, "IMPLICIT BULK", 0, Region::bulk);
Region RegionDB::implicit_boundary=Region::db().add_region(Region::undefined-2, "IMPLICIT BOUNDARY", 0, Region::boundary);


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
                return Region(index);
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
}



Region RegionDB::add_region(unsigned int id, unsigned int dim) {
    Region r_id = find_id(id);
    if (r_id.is_valid()) return add_region(id, r_id.label(), dim, r_id.is_boundary() );
    // else
    stringstream ss;
    ss << "region_" << id;
    return add_region(id, ss.str(), dim, false);
}



Region RegionDB::find_label(const std::string &label)
{
    RegionTable::index<Label>::type::iterator it_label = region_set_.get<Label>().find(label);
    if (it_label==region_set_.get<Label>().end()  ) return Region();
    return Region(it_label->index);
}



Region RegionDB::find_id(unsigned int id)
{
    RegionTable::index<ID>::type::iterator it_id = region_set_.get<ID>().find(id);
    if ( it_id==region_set_.get<ID>().end() ) return Region();
    return Region(it_id->index);
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
}


unsigned int RegionDB::size() {
    if (! closed_) close();
    return 2* max(n_boundary_, n_bulk_);
}



unsigned int RegionDB::boundary_size() {
    if (! closed_) close();
    return n_boundary_;
}



unsigned int RegionDB::bulk_size() {
    if (! closed_) close();
    return n_bulk_;
}



const RegionDB::RegionSet &RegionDB::boundary_regions() {
    if (boundary.size() == 0) {
        for(unsigned int i=0; i< size(); i++) {
            Region reg(i);
            if (reg.is_boundary() && reg.boundary_idx() < n_boundary_) boundary.push_back(reg);
        }
    }
    return boundary;
}


void RegionDB::add_to_set( const string& set_name, RegionIdx region) {
	std::map<std::string, RegionDB::RegionSet>::iterator it = sets_.find(set_name);

	if (it == sets_.end()) {
		RegionDB::RegionSet set;
		set.push_back(region);

		sets_.insert( std::make_pair(set_name, set) );
	} else {
		RegionDB::RegionSet & set = (*it).second;
		if ( std::find(set.begin(), set.end(), region)!=set.end() ) {
			set.push_back(region); // add region if doesn't exist
		}
	}
}


void RegionDB::add_set( const string& set_name, const RegionDB::RegionSet & set) {
	if (sets_.find(set_name) != sets_.end()) {
		sets_.erase(sets_.find(set_name));
	}

	sets_.insert( std::make_pair(set_name, set) );
}


RegionDB::RegionSet RegionDB::union_sets( const string & set_name_1, const string & set_name_2) {
	std::map<std::string, RegionDB::RegionSet>::iterator it_1 = sets_.find(set_name_1);
	std::map<std::string, RegionDB::RegionSet>::iterator it_2 = sets_.find(set_name_2);
	ASSERT( it_1 == sets_.end(), "No region set with name: %s\n", set_name_1);
	ASSERT( it_2 == sets_.end(), "No region set with name: %s\n", set_name_2);

	RegionDB::RegionSet set_union;
	RegionDB::RegionSet & set_1 = (*it_1).second;
	RegionDB::RegionSet & set_2 = (*it_2).second;
	RegionDB::RegionSet::iterator it;

	std::sort(set_1.begin(), set_1.end());
	std::sort(set_2.begin(), set_2.end());
	set_union.resize(set_1.size() + set_2.size());
	it = std::set_union(set_1.begin(), set_1.end(), set_2.begin(), set_2.end(), set_union.begin());
	set_union.resize(it - set_union.begin());

	return set_union;
}


RegionDB::RegionSet RegionDB::intersection( const string & set_name_1, const string & set_name_2) {
	// not implemented yet
	return RegionDB::RegionSet();
}


RegionDB::RegionSet RegionDB::difference( const string & set_name_1, const string & set_name_2) {
	// not implemented yet
	return RegionDB::RegionSet();
}


const RegionDB::RegionSet & RegionDB::get_region_set(const string & set_name) {
	std::map<std::string, RegionSet>::iterator it = sets_.find(set_name);
	ASSERT( it == sets_.end(), "No region set with name: %s\n", set_name);
	return (*it).second;
}


void RegionDB::read_sets_from_input(Input::Record rec) {
	// not implemented yet
}
