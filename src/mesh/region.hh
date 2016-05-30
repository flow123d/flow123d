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
 * @file    region.hh
 * @brief   
 * @todo
 *  - komentar k RegionIdx
 *  - presun enum RegionType do public Region - komentar + pouzit v kodu
 *  - zkontrolovat chybove hlasky a ASSERTY, co z toho by melo byt pres exception?
 *
 *  - Seems that GMSH allows repeating ID and Label on regions of different dimension, therefore
 *    label and ID are not unique without dimension.
 *
 */

#ifndef REGION_HH_
#define REGION_HH_

#include <string>
#include <vector>
#include <map>

#include "system/system.hh"
#include "system/global_defs.h"
#include "system/exceptions.hh"

#include "input/accessors_forward.hh"
#include "input/input_exception.hh"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/functional/hash.hpp>

namespace BMI=::boost::multi_index;

// forward declarations
class Region;
class RegionDB;
namespace Input {
    namespace Type { class Record; }
    class Record;
    class Array;
}

/**
 * Base class that contains information about region:
 * 1) contains integer value that specifies region
 * 2) detects if region is valid
 * 3) detects if the region is the bulk or boundary
 * 4) detects boundary index of boundary region or bulk index of bulk region
 */
class RegionIdx {
public:

    /// Default region is undefined/invalid
	RegionIdx():idx_(undefined) {}

    /// Returns true if it is a Boundary region and false if it is a Bulk region.
    inline bool is_boundary() const
        { return !(idx_ & 1); }

    /// Returns false if the region has undefined/invalid value
    inline bool is_valid() const
        { return idx_!=undefined;}

    /// Returns a global index of the region.
    inline unsigned int idx() const
        { return idx_; }

    /// Returns index of the region in the boundary set.
    inline unsigned int boundary_idx() const {
    	OLD_ASSERT( is_boundary(), "Try to get boundary index of a bulk region with internal idx: %d\n", idx_ );
        return idx_ >> 1; }

    /// Returns index of the region in the bulk set.
    inline unsigned int bulk_idx() const {
    	OLD_ASSERT( ! is_boundary(), "Try to get bulk index of boundary region with internal idx: %d\n", idx_ );
        return idx_ >> 1; }

    /// Equality comparison operators for regions.
    inline bool operator==(const RegionIdx &other) const
        { return idx_ == other.idx_; }

    /// Equality comparison operators for regions.
    inline bool operator!=(const RegionIdx &other) const
        { return idx_ != other.idx_; }


protected:
    /**
     * Create accessor from the index. Should be private since implementation specific.
     * We need some way how to iterate over: all regions, boundary regions, bulk regions -
     * solution: have specific RegionSets for these three cases.
     */
    RegionIdx(unsigned int index)
    : idx_(index) { }

    /**
     * Internal region index. Regions of one RegionDB (corresponding to one mesh) forms more or less continuous sequence.
     */
	unsigned int idx_;

	/// index for undefined region
    static const unsigned int undefined=0xffffffff;

    friend class Region;
};



/**
 * Class that represents disjoint part of computational domain (or domains). It consists of one integer value
 * but provides access to other data stored in RegionDB. In particular provides string label and integer ID (unordered)
 * further it provides fast (inlined) methods to:
 * 1) detect if the region is the bulk region or boundary region
 * 2) return index (this is used to select correct Field, possibly we can distinguish boundary_index and bulk_index)
 *
 * Implementation: currently we number bulk regions by odd indices and boundary regions by even indices.
 *
 */
class Region : public RegionIdx {
public:

	/**
	 * Types of region in mesh (bulk or boundary)
	 */
    enum RegionType {
        bulk=false,
        boundary=true
    };


    /// Default region is undefined/invalid
    Region()
    : RegionIdx(), db_(NULL)
    {}

    /// This should be used for construction from known RegionIdx. (e.g. in Mesh)
    /// Do not use unless you can not get Region in other way.
    Region(RegionIdx r_idx, const RegionDB & db)
    : RegionIdx(r_idx), db_(&db)
    {}

    RegionIdx operator() (const Region &)
        {return RegionIdx(idx_); }

    /// Comparative method of two regions
    static bool comp(const Region &a, const Region &b)
    { return a.idx_ < b.idx_; }

    /// Returns label of the region (using RegionDB)
    std::string label() const;

    /// Returns id of the region (using RegionDB)
    unsigned int id() const;

    /// Returns dimension of the region.
    unsigned int dim() const;

    /**
     * Returns static region database. Meant to be used for getting range of
     * global, boundary, and bulk region indices.
     */
    const RegionDB &db() {
        return *db_;
    }

protected:
    /**
     * Create accessor from the index. Should be private since implementation specific.
     * We need some way how to iterate over: all regions, boundary regions, bulk regions -
     * solution: have specific RegionSets for these three cases.
     */
    Region(unsigned int index, const RegionDB &db)
    : RegionIdx(index), db_(&db)
    {}

    /// Global variable with information about all regions.
    const RegionDB *db_;

    friend class RegionDB;
    friend class Mesh;
};




/**
 * Type representing a set of regions.
 * CAn be used  to set function(field) on more regions at once, possibly across meshes
 *
 * Regions stored in region set are always unique
 */
typedef std::vector<Region> RegionSet;
/// Type representing a map of RegionSets.
typedef std::map<std::string, RegionSet > RegionSetTable;


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
 * TODO:
 * In order to support more meshes , possibly changing during the time we need better policy for RegionDB closing.
 * Currently we close RegionDB at first call to any of @p size methods. We need size information for initialization of
 * RegionFields. However, every RegionField should be used for assembly over just one mesh (or set of meshes - supermesh?),
 * surly this mesh has to be initialized before assembly so it could be initialized before initialization of RegionField which
 * lives on this mesh. So the solution can be: RegionDB keeps list of meshes that has their regions registered in RegionDB.
 * RegionField has signature of its mesh and check if the mesh is registered in the RegionDB before initialization of REgionField.
 *
 * Update TODO:
 *
 * - make one instance of Region DB per mesh; put region indices into mesh elements make them private and provide
 *   method that returns Region accessor. ?? Speed concerns? This introduce needless copies of mesh pointer since
 *   in most cases we do not need methods label, id, dim where whole RegionnDB is needed. Thus possibly make
 *   class RegionIdx (without Mesh) and make some asking methods for label, id in RegionDB.
 * - in RegionDB make support for sets:

/// map set name to lists of indices of its regions
typedef std::vector<RegionIdx> RegionSet;
std::map<std::string, RegionSet > sets_;

/// Add region to given set. Creat the set if it does not exist.
add_to_set( const string& set_name, RegionIdx region);
/// Add a set into map, delete possible previous value, do not worry about slow copies of
/// the set array.
add_set( const string& set_name, const RegionSet & set);
RegionSet union( const string & set_name_1, const string & set_name_2); // sort + std::set_union
RegionSet intersection( const string & set_name_1, const string & set_name_2);
RegionSet difference( const string & set_name_1, const string & set_name_2);
const RegionSet & get_set(const string & set_name);
void read_sets_from_input(Input::Record rec); // see structure of RegionDB::get_region_set_input_type()

region_sets = [
   { name="set name",
     region_ids=[ int ...],
     region_labels= ["..."], // these are merger together

     union=["set_1", "set_2"], // later overwrites previous
     intersection=
     difference=
   }
]
 *
 * - MODIFICATION IN MESH AND READER:
 * Mesh reading proccess:
 * 1) Read PhysicalNames form GMSH file, populate RegionDB (DONE in GMSH reader, may need small modifications)
 * 2) Read region definitions from the input, see
 *    - Mesh::read_regions_from_input(Input::Array region_list);
 * 3) Read nodes (DONE in GMSH reader)
 * 4) Read elements, per element:
 *    - possibly modify region according map
 *    - find element ID:
 *       if found: add_region(ID, get_label, dim, get_boundary_flag) // possibly set dimension of the region if it is undefined
 *       not found: add_region(ID, default_label, dim, false)
 *    - if region is boundary put element into Mesh::bc_elements
 *      else put it into Mesh::element
 *  ---
 *  5) Setup topology - we has to connect Boundary with existing bc_elements, and add the remaining elements,
 *     after we remove support for old bCD files we may skip creation of remaining boundary elements since there will be no way how to set
 *     BC on them.
 *
 */

class RegionDB {
public:
    /**
     * Map representing the relevance of elements to regions
     */
    typedef std::map<unsigned int, unsigned int> MapElementIDToRegionID;

    TYPEDEF_ERR_INFO( EI_Label, const std::string);
    TYPEDEF_ERR_INFO( EI_ID, unsigned int);
    TYPEDEF_ERR_INFO( EI_IDOfOtherLabel, unsigned int);
    TYPEDEF_ERR_INFO( EI_LabelOfOtherID, const std::string);
    DECLARE_INPUT_EXCEPTION( ExcAddingIntoClosed, << "Can not add label=" << EI_Label::qval << " into closed MaterialDispatch.\n");
    DECLARE_EXCEPTION( ExcNonuniqueID, << "Non-unique ID during add of elementary region id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "other elementary region with same ID but different label: " << EI_LabelOfOtherID::qval << " already exists\n");
    DECLARE_INPUT_EXCEPTION( ExcNonuniqueLabel, << "Non-unique label during add of elementary region id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "other elementary region with same label but different ID: " << EI_IDOfOtherLabel::val << " already exists\n");
    DECLARE_EXCEPTION( ExcInconsistentBoundary, << "Inconsistent add of elementary region with id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "both ID and label match an existing elementary region with different boundary flag.");

    DECLARE_INPUT_EXCEPTION( ExcCantAdd, << "Can not add new elementary region into DB, id: " << EI_ID::val <<", label: " << EI_Label::qval);

    DECLARE_INPUT_EXCEPTION( ExcUnusedRegion, << "Region with id: " << EI_ID::qval << " and label: " << EI_Label::qval
    									<< " is not used in any element." );

    DECLARE_INPUT_EXCEPTION( ExcUnknownRegion, << "Unknown region with id: " << EI_ID::val );

    DECLARE_INPUT_EXCEPTION( ExcUnknownSet, << "Operation with unknown region set: " << EI_Label::qval );

    DECLARE_INPUT_EXCEPTION( ExcUnknownSetOperand, << "Operation with unknown region set: " << EI_Label::qval);

    DECLARE_INPUT_EXCEPTION(ExcUniqueRegionId, << "Id of elementary region must be unique, id: " << EI_ID::val );

    /// Default constructor
    RegionDB();

    /**
     * Introduce an artificial limit to keep all material indexed arrays
     * of reasonable size.
     */
    static const unsigned int max_n_regions = 64000;

    /// Undefined dimension for regions introduced from mesh input record.
    /// Dimensions 0,1,2,3 are valid.
    static const unsigned int undefined_dim;


    /**
     * This method adds new region into the database and returns its index. This requires full
     * specification of the region that is given in PhysicalNames section of the GMSH MSH format
     * or in Mesh input record.
     *
     * If ID or label are found in the DB, we distinguish following cases:
     * 1)  ID is found, label is not found : warning ID has already assigned label
     * 2)  ID is not found, label is found : report error - assigning same label to different IDs
     * 3)  both ID and label are found, in same region : check remaining data, return existing region
     * 4)                             , in different   : warning ID has already assigned label
     *
     * Parameter @p id is any unique non-negative integer, parameter @p label is unique string identifier of the region,
     * @p dim is dimension of reference elements in the region, @p boundary is true if the region consist of boundary elements
     * (where one can apply boundary condition) and @p address contains source of region (address in input file or section in
     * mesh file).
     */
    Region add_region(unsigned int id, const std::string &label, unsigned int dim, const std::string &address ="implicit");

    /**
     * Change label of given Region.
     */
    Region rename_region( Region reg, const std::string &new_label );

    /**
     * Returns region given the pair of id - dim.
     * If region doesn't exist, checks if exists region with given id and undefined_dim, replaces
     * its dimension and returns its. In other cases throws exception.
     */
    Region get_region(unsigned int id, unsigned int dim);

    /**
     * Returns a @p Region with given @p label. If it is not found it returns @p undefined Region.
     */
    Region find_label(const std::string &label) const;

    /**
     * Returns a @p Region with given @p id. If it is not found it returns @p undefined Region.
     * Gmsh ID numbers are unique only over one dimension, so dimension @p dim must be provided as well.
     */
    Region find_id(unsigned int id, unsigned int dim) const;

    /**
     * Slower version that tries to find region for given ID. If it is not unique it throws.
     */
    Region find_id(unsigned int id) const;

    /**
     * Return original label for given index @p idx.
     */
    const std::string &get_label(unsigned int idx) const;

    /**
     * Return original ID for given index @p idx.
     */
    unsigned int get_id(unsigned int idx) const;

    /**
     * Return dimension of region with given index @p idx.
     */
    unsigned int get_dim(unsigned int idx) const;

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

    /**
     * Returns implicit boundary region. Is used for boundary elements created by Flow123d itself.
     * This region has label "IMPLICIT_BOUNDARY" and it is obsolete, the name is not consistent
     * with boundary label notation.
     */
    Region implicit_boundary_region();

    /*
     * Add region to given set. Create the set if it does not exist.
     *
     * @param set_name Set to which it is added region
     * @param region Added region
     */
    void add_to_set( const string& set_name, Region region);

    /**
     * Add a set into map, delete possible previous value.
     *
     * @param set_name Name of added set
     * @param set Added RegionSet
     */
    void add_set( const string& set_name, const RegionSet & set);

    /**
     * Get region set of specified name. Three sets are defined by default:
     * "ALL" - set of all regions both bulk and boundary.
     * "BULK" - set of all bulk regions
     * "BOUNDARY" - set of all boundary regions
     *
     * @param set_name Name of set
     * @return RegionSet of specified name. Returns Empty vector if the set of given name doesn't exist.
     */
    RegionSet get_region_set(const string & set_name) const;

    /**
     * Read two operands from input array of strings and check if given names
     * are existing sets. Return pair of checked set names.
     */
    std::vector<string> get_and_check_operands(const Input::Array & operands) const;

    /**
     * Print table with base information of all regions stored in RegionDB.
     */
    void print_region_table(ostream& stream) const;

    /**
     * Create label of region in format: "region_"+id
     *
     * Use if label is not set.
     */
    string create_label_from_id(unsigned int id) const;

    /**
     * Return address for given index @p idx.
     */
    const std::string & get_region_address(unsigned int idx) const;

    /**
     * Mark region with given index @p idx as used.
     *
     * Use if region is assigned to element.
     */
    void mark_used_region(unsigned int idx);

    /**
     * Create union of RegionSets of given names defined in @p set_names.
     */
    RegionSet union_set(std::vector<string> set_names) const;


private:


    typedef std::pair<unsigned int, unsigned int> DimID;

    /// One item in region database
    struct RegionItem {
        RegionItem(unsigned int index, unsigned int id, const std::string &label, unsigned int dim, const std::string &address, bool used=false)
            : index(index), id(dim, id), label(label), used(used), address(address) {}

        unsigned int get_id() const {return id.second;}
        unsigned int dim() const {return id.first;}

        // unique identifiers
        unsigned int index;
        DimID id;
        std::string label;
        // Flag signed if region is assigned to element(s)
        bool used;
        // Address where region was created (address in input file or section in mesh file)
        std::string address;
    };

    // tags
    struct DimId {};
    struct OnlyID {};
    struct Label {};
    struct Index {};

    /// Region database
    typedef BMI::multi_index_container<
            RegionItem,
            BMI::indexed_by<
                // access by index
                // Can not use random access without introducing "empty" RegionItems to fill holes in the case
                // n_boundary != n_bulk. Empty items must be unsince we may have empty (and unmodifiable) holes ?? why)
                //
                // we need O(1) access
                //BMI::random_access< BMI::tag<RandomIndex > >,
                BMI::ordered_unique< BMI::tag<Index>, BMI::member<RegionItem, unsigned int, &RegionItem::index > >,
                // use hashing for IDs, to get O(1) find complexity .. necessary for large meshes
                BMI::hashed_unique< BMI::tag<DimId>,    BMI::member<RegionItem, DimID, &RegionItem::id> >,
                // non unique index for sole ID
                BMI::hashed_non_unique< BMI::tag<OnlyID>,    BMI::const_mem_fun<RegionItem, unsigned int, &RegionItem::get_id> >,
                // ordered access (like stl::map) by label
                BMI::ordered_unique< BMI::tag<Label>, BMI::member<RegionItem, std::string, &RegionItem::label> >
            >
    > RegionTable;

    // DimID and Label index iterators
    typedef RegionTable::index<Label>::type::iterator  LabelIter;
    typedef RegionTable::index<DimId>::type::iterator  DimIDIter;
    typedef RegionTable::index<OnlyID>::type::iterator OnlyIDIter;


    /// Database of all regions (both boundary and bulk).
    RegionTable region_table_;

    /// flag for closed database, no regions can be added, but you can add region sets
    bool closed_;
    /// Number of boundary regions
    unsigned int n_boundary_;
    /// Number of bulk regions
    unsigned int n_bulk_;
    /// Maximal value of Region::id()
    unsigned int max_id_;

    /// Map of region sets
    RegionSetTable sets_;

    /// Make part of general RegionSet table.
    RegionSet all, bulk, boundary;

    /**
     * Implicit bulk and boundary regions. For GMSH mesh format only implicit_boundary is used for boundary elements
     * that are not explicitly in the mesh file.
     */
    Region implicit_bulk_, implicit_boundary_;

    /**
     * Represents the relevance of elements to regions. Defined by user in input file.
     */
    MapElementIDToRegionID el_to_reg_map_;

    /**
     * Insert new region into database.
     */
    Region insert_region(unsigned int id, const std::string &label, unsigned int dim, bool boundary, const std::string &address);

    /**
     * Replace dimension of existing region with undefined_dim.
     */
    Region replace_region_dim(DimIDIter it_undef_dim, unsigned int dim, bool boundary);

    /**
     * Find existing region given by pair (dim, id).
     */
    Region find_by_dimid(DimIDIter it_id, unsigned int id, const std::string &label, bool boundary);

    /*
     * Add region to given set. Create the set if it does not exist.
     *
     * @param set_name Set from which it is erased region
     * @param region Erased region
     */
    void erase_from_set( const string& set_name, Region region);

    /**
     * Iterate all stored regions and check if regions are assigned to element(s).
     *
     * Unused region throws exception.
     */
    void check_regions();

    /**
     * Return boundary flag for given label. Label of boundary region must start by '.' symbol.
     */
    inline bool is_boundary(const std::string &label) {
    	return (label.size() != 0) && (label[0] == '.');
    }

    friend class Mesh;
    friend class RegionSetBase;
};



#endif /* REGION_HH_ */
