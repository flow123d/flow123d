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

Region RegionDB::add_region( unsigned int id, const std::string &label, unsigned int dim, bool boundary) {
    if (closed_) xprintf(PrgErr, "Can not add to closed region DB.\n");

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
            if ( ! region_set_.insert( RegionItem(index, id, label, dim) ).second )
               THROW( ExcCantAdd()  << EI_Label(label) <<EI_ID(id) );
            return Region(index);
        }
    } else {
        // both iterators are valid
        if (it_id->index != it_label->index) THROW(ExcInconsistentAdd() << EI_Label(label) << EI_ID(id) );
        Region reg(it_id->index);
        if ( reg.dim() != dim || reg.is_boundary() != boundary ) THROW(ExcInconsistentAdd() << EI_Label(label) << EI_ID(id) );
        return reg;
    }
}



Region RegionDB::add_region(unsigned int id, unsigned int dim) {
    stringstream ss;
    ss << "region_%d" << id;
    return add_region(id, ss.str(), dim, false);
}



Region RegionDB::find_label(const std::string &label)
{
    RegionSet::index<Label>::type::iterator it_label = region_set_.get<Label>().find(label);
    if (it_label==region_set_.get<Label>().end()  ) return Region();
    return Region(it_label->index);
}



Region RegionDB::find_id(unsigned int id)
{
    RegionSet::index<ID>::type::iterator it_id = region_set_.get<ID>().find(id);
    if ( it_id==region_set_.get<ID>().end() ) return Region();
    return Region(it_id->index);
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



unsigned int RegionDB::get_dim(unsigned int idx) const {
    RegionSet::index<Index>::type::iterator it = region_set_.get<Index>().find(idx);
    ASSERT( it!= region_set_.get<Index>().end(), "No region with index: %d\n", idx);
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

