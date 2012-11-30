/*
 * material_dispatch.cc
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#include <string>

#include "mesh/region.hh"
#include "system/exceptions.hh"

using namespace std;

RegionDB Region::db_;


std::string Region::label() const
    { return db_.get_label(idx_); }



unsigned int Region::id() const
    { return db_.get_id(idx_); }



/**************************************************************************************************
 * Implementation of     RegionDB
 */

Region RegionDB::add_region( unsigned int id, const std::string &label, bool boundary) {
    RegionSet::index<ID>::type::iterator it_id = region_set_.get<ID>().find(id);
    RegionSet::index<Label>::type::iterator it_label = region_set_.get<Label>().find(label);

    if (it_id==region_set_.get<ID>().end() || it_label==region_set_.get<Label>().end()  ) {
        // check that both finds failed
        if (! ( it_id==region_set_.get<ID>().end() && it_label==region_set_.get<Label>().end() ))
            THROW(ExcInconsistentAdd() << EI_Label(label) << EI_ID(id) );
        // if DB is open add new entry
        if (closed_)  THROW( ExcAddingIntoClosed() << EI_Label(label) <<EI_ID(id) );
        else {
            unsigned int index;
            if (boundary) (index = (n_boundary_ <<1)), n_boundary_++;
            else (index = (n_bulk_ << 1)+1),  n_bulk_++;
            if ( ! region_set_.insert( RegionItem(index, id, label) ).second )
               THROW( ExcCantAdd()  << EI_Label(label) <<EI_ID(id) );
            return Region(index);
        }
    } else {
        // both iterators are valid
        if (it_id->index != it_label->index) THROW(ExcInconsistentAdd() << EI_Label(label) << EI_ID(id) );
        return Region(it_id->index);
    }
}



const std::string & RegionDB::get_label(unsigned int idx) const {
    RegionSet::index<Index>::type::iterator it = region_set_.get<Index>().find(idx);
    ASSERT( it!= region_set_.get<Index>().end(), "No region with index: %d\n", idx);
    return  it->label;
}



unsigned int RegionDB::get_id(unsigned int idx) const {
    RegionSet::index<Index>::type::iterator it = region_set_.get<Index>().find(idx);
    ASSERT( it!= region_set_.get<Index>().end(), "No region with index: %d\n", idx);
    return  it->id;
}



void RegionDB::close() {
    closed_=true;
}


unsigned int RegionDB::size() const {
    if (! closed_) THROW( ExcSizeWhileOpen() );
    return 2* max(n_boundary_, n_bulk_);
}



unsigned int RegionDB::boundary_size() const {
    if (! closed_) THROW( ExcSizeWhileOpen() );
    return n_boundary_;
}



unsigned int RegionDB::bulk_size() const {
    if (! closed_) THROW( ExcSizeWhileOpen() );
    return n_bulk_;
}


