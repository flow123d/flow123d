/*
 * region.hh
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#ifndef REGION_HH_
#define REGION_HH_

#include <string>
#include <vector>
#include <map>

#include "system/system.hh"
#include "system/global_defs.h"
#include "system/exceptions.hh"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/member.hpp>

namespace BMI=::boost::multi_index;

class RegionDB;

/**
 * Class that represents disjoint part of computational domain (or domains). It consists of one integer value
 * but provides access to other data stored in RegionDB. In particular provides string label and integer ID (unordered)
 * further it provides fast (inlined) methods to:
 * 1) detect if the region is the bulk region or boundary region
 * 2) return index (this is used to select correct Field, possibly we can distinguish boundary_index and bulk_index)
 *
 * Implementation: currently we number bulk regions by odd indices and boundary regions by even indices.
 */
class Region {
public:

    /// Returns true if it is a Boundary region and false if it is a Bulk region.
    inline bool is_boundary() const
        { return !(idx_ & 1); }

    /// Returns a global index of the region.
    inline unsigned int idx() const
        { return idx_; }

    /// Returns index of the region in the boundary set.
    inline unsigned int boundary_idx() const {
        ASSERT( is_boundary(), "Try to get boundary index of bulk region: '%s' id: %d\n", label().c_str(), id() );
        return idx_ >> 1; }

    /// Returns index of the region in the bulk set.
    inline unsigned int bulk_idx() const {
        ASSERT( ! is_boundary(), "Try to get bulk index of boundary region: '%s' id: %d\n", label().c_str(), id() );
        return idx_ >> 1; }

    /// Returns label of the region (using RegionDB)
    std::string label() const;

    /// Returns id of the region (using RegionDB)
    unsigned int id() const;

    /**
     * Returns region database. Meant to be used for getting range of
     * global, boundary, and bulk region indices.
     */
    static RegionDB &db()
        { return db_;}

private:
    /// Global variable with information about all regions.
    static RegionDB db_;

    /**
     * Create accessor from the index. Private since implementation specific.
     */
    Region(unsigned int index)
    : idx_(index) {}

    unsigned int idx_;

    friend class RegionDB;
};


/**
 * Class representing a set of regions.
 * CAn be used  to set function(field) on more regions at once, possibly across meshes
 *
 * Desired properties:
 * - can construct itself from input, from a list
 *   of regions (given by label or id)
 * - support set operations
 * - say if an region is in it
 * - iterate through its regions
 */
class RegionSet {
};


/**
 * Class for conversion between an index and string label of an material.
 * Class contains only static members in order to provide globally consistent
 * indexing of materials across various meshes and functions.
 *
 * The conversion should be performed only through the input and output so
 * that the lookup overhead could be shadowed by IO operations.
 *
 * Taking sizes and creating region sets should be possible only after the database is closed.
 * We assume that all regions are known at beginning of the program (typically after reading all meshes)
 * however they need not be used through the whole computation.
 *
 * TODO: Use boost multi_index_set to avoid two complementary data structures.
 *
 */

class RegionDB {
public:
    TYPEDEF_ERR_INFO( EI_Label, const std::string);
    TYPEDEF_ERR_INFO( EI_ID, unsigned int);
    DECLARE_EXCEPTION( ExcAddingIntoClosed, << "Can not add label=" << EI_Label::qval << " into closed MaterialDispatch.\n");
    DECLARE_EXCEPTION( ExcSizeWhileOpen, << "Can not get size of MaterialDispatch yet open.");
    DECLARE_EXCEPTION( ExcInconsistentAdd, << "Can get region with id: " << EI_ID::val <<", label: " << EI_Label::qval
                                             << " each used by different region.");
    DECLARE_EXCEPTION( ExcCantAdd, << "Can not add new region into DB, id: " << EI_ID::val <<", label: " << EI_Label::qval);


    /// Default constructor
    RegionDB()
    : closed_(false), n_boundary_(0), n_bulk_(0)  {}

    /**
     * Introduce an artificial limit to keep all material indexed arrays
     * of reasonable size.
     */
    static const unsigned int max_n_regions = 64000;
    /**
     * Add new region into database and return its index. If the region is already in the DB,
     * check consistency of label and id and return its index.
     */
    Region add_region(unsigned int id, const std::string &label, bool boundary);
    /**
     * Return original label for given index @p idx.
     */
    const std::string &get_label(unsigned int idx) const;
    /**
     * Return original ID for given index @p idx.
     */
    unsigned int get_id(unsigned int idx) const;

    /**
     * Close this class for adding labels. This is necessary to return correct size
     * for material indexed arrays and vectors. After calling this method you can
     * call method @p size and method @p idx_of_label rise an exception for any unknown label.
     */
    void close();
    /**
     * Returns maximal index + 1
     */
    unsigned int size() const;
    /**
     * Returns total number boundary regions.
     */
    unsigned int boundary_size() const;
    /**
     * Returns total number bulk regions.
     */
    unsigned int bulk_size() const;


private:
    /// One item in region database
    struct RegionItem {
        RegionItem(unsigned int index, unsigned int id, const std::string &label)
            : index(index), id(id), label(label) {}

        unsigned int index;
        unsigned int id;
        std::string label;
    };

    // tags
    struct ID {};
    struct Label {};
    struct Index {};

    /// Region database
    typedef BMI::multi_index_container<
            RegionItem,
            BMI::indexed_by<
                // access by index (can not use random access since we may have empty (and unmodifiable) holes)
                BMI::ordered_unique< BMI::tag<Index>, BMI::member<RegionItem, unsigned int, &RegionItem::index > >,
                // ordered access (like stl::map) by id and label
                BMI::ordered_unique< BMI::tag<ID>,    BMI::member<RegionItem, unsigned int, &RegionItem::id> >,
                BMI::ordered_unique< BMI::tag<Label>, BMI::member<RegionItem, std::string, &RegionItem::label> >
            >
    > RegionSet;

    RegionSet region_set_;

    bool closed_;
    /// Number of boundary regions
    unsigned int n_boundary_;
    /// Number of bulk regions
    unsigned int n_bulk_;
};



#endif /* REGION_HH_ */
