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

#include "input/input_type.hh"
#include "input/accessors.hh"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/member.hpp>

namespace BMI=::boost::multi_index;

/**
 *
 */
class RegionIdx {
public:
    /**
     * Create accessor from the index. Should be private since implementation specific.
     * We need some way how to iterate over: all regions, boundary regions, bulk regions -
     * solution: have specific RegionSets for these three cases.
     */
	RegionIdx(unsigned int index)
    : idx_(index) {}

    /// Default region is undefined/invalid
	RegionIdx():idx_(undefined) {}

    /// Returns false if the region has undefined/invalid value
    inline bool is_valid() const
        { return idx_!=undefined;}

    /// Returns a global index of the region.
    inline unsigned int idx() const
        { return idx_; }

    /// Equality comparison operators for regions.
    inline bool operator==(const RegionIdx &other) const
        { return idx_ == other.idx_; }

    /// Equality comparison operators for regions.
    inline bool operator!=(const RegionIdx &other) const
        { return idx_ != other.idx_; }

    /// Compare operator for regions.
    inline bool operator<(const RegionIdx &other) const
        { return idx_ < other.idx_; }

    /// Compare operator for regions.
    inline bool operator<=(const RegionIdx &other) const
        { return idx_ <= other.idx_; }

    /// Compare operator for regions.
    inline bool operator>(const RegionIdx &other) const
        { return idx_ > other.idx_; }

    /// Compare operator for regions.
    inline bool operator>=(const RegionIdx &other) const
        { return idx_ >= other.idx_; }

protected:
	unsigned int idx_;

	/// index for undefined region
    static const unsigned int undefined=0xffffffff;
};

class RegionDB;
class OldBcdInput;

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
    enum RegionType {
        bulk=false,
        boundary=true
    };

    /**
     * Create accessor from the index. Should be private since implementation specific.
     * We need some way how to iterate over: all regions, boundary regions, bulk regions -
     * solution: have specific RegionSets for these three cases.
     */
    Region(unsigned int index)
    : RegionIdx::RegionIdx(index) {}

    /// Default region is undefined/invalid
    Region() : RegionIdx::RegionIdx() {}

    /// Returns true if it is a Boundary region and false if it is a Bulk region.
    inline bool is_boundary() const
        { return !(idx_ & 1); }

/*    /// Returns false if the region has undefined/invalid value
    inline bool is_valid() const
        { return idx_!=undefined;}

    /// Returns a global index of the region.
    inline unsigned int idx() const
        { return idx_; }*/

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

    /// Returns dimension of the region.
    unsigned int dim() const;

    /// Equality comparison operators for regions.
    inline bool operator==(const Region &other) const
        { return idx_ == other.idx_; }

    /// Equality comparison operators for regions.
    inline bool operator!=(const Region &other) const
        { return idx_ != other.idx_; }

    /**
     * Returns static region database. Meant to be used for getting range of
     * global, boundary, and bulk region indices.
     */
    static RegionDB &db()
        { return db_;}

private:
    /// Global variable with information about all regions.
    static RegionDB db_;



    //unsigned int idx_;

    friend class RegionDB;
    friend class OldBcdInput;
};


/**
 * Class representing a set of regions.
 * CAn be used  to set function(field) on more regions at once, possibly across meshes
 *
 * TODO:
 * Desired properties:
 * - can construct itself from input, from a list
 *   of regions (given by label or id)
 * - support set operations
 * - say if an region is in it
 * - iterate through its regions
 */
//class RegionSet {
//};


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
void read_sets_from_input(Input::Record rec); // see structure of RegionDB::region_set_input_type

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
 * 2) Read region definitions from the input, see RegionDB::region_input_type
 *    (TODO in RegionDB, also creates (and return to Mesh) element regions modification map: std::map< unsigned int, RegionIdx>
 *     that maps element IDs to the new region names, GMSH reader should have setter method to accept this map
 *     and modify the elements during reading)
 * 3) Read region sets - TODO in RegionDB
 * 4) Read boundary key of the Mesh record and mark appropriate regions as boundary (TODO in RegionDB)
 * 5) Read nodes (DONE in GMSH reader)
 * 6) Read elements, per element:
 *    - possibly modify region according map
 *    - find element ID:
 *       if found: add_region(ID, get_label, dim, get_boundary_flag) // possibly set dimension of the region if it is undefined
 *       not found: add_region(ID, default_label, dim, false)
 *    - if region is boundary put element into Mesh::bc_elements
 *      else put it into Mesh::element
 *  ---
 *  7) Setup topology - we has to connect Boundary with existing bc_elements, and add the remaining elements,
 *     after we remove support for old bCD files we may skip creation of remaining boundary elements since there will be no way how to set
 *     BC on them.
 *
 */

class RegionDB {
public:
    typedef std::vector<RegionIdx> RegionSet;

    static Region implicit_bulk;
    static Region implicit_boundary;

    static Input::Type::Record region_input_type;
    static Input::Type::Record region_set_input_type;

    TYPEDEF_ERR_INFO( EI_Label, const std::string);
    TYPEDEF_ERR_INFO( EI_ID, unsigned int);
    TYPEDEF_ERR_INFO( EI_IDOfOtherLabel, unsigned int);
    TYPEDEF_ERR_INFO( EI_LabelOfOtherID, const std::string);
    DECLARE_EXCEPTION( ExcAddingIntoClosed, << "Can not add label=" << EI_Label::qval << " into closed MaterialDispatch.\n");
    DECLARE_EXCEPTION( ExcSizeWhileOpen, << "Can not get size of MaterialDispatch yet open.");
    DECLARE_EXCEPTION( ExcInconsistentAdd, << "Inconsistent add of region with id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "other region with same ID but different label: " << EI_LabelOfOtherID::qval << " already exists\n" \
                                             << "OR other region with same label but different ID: " << EI_IDOfOtherLabel::val << " already exists\n" \
                                             << "OR both ID and label match an existing region with different dimension and/or boundary flag.");
    DECLARE_EXCEPTION( ExcCantAdd, << "Can not add new region into DB, id: " << EI_ID::val <<", label: " << EI_Label::qval);


    /// Default constructor
    RegionDB();

    /**
     * Introduce an artificial limit to keep all material indexed arrays
     * of reasonable size.
     */
    static const unsigned int max_n_regions = 64000;

    /**
     * Add new region into database and return its index. This requires full
     * specification of the region that is given in PhysicalNames section of the GMSH MSH format.
     * If the region is already in the DB, check consistency of label and id and return its index.
     *
     * Parameter @p id is any non-negative integer that is unique for the region over all meshes used in the simulation,
     * parameter @p label is unique string identifier of the region, @p dim is dimension of reference elements in the region
     * and @p boundary is true if the region consist of boundary elements (where one can apply boundary condition).
     *
     */
    Region add_region(unsigned int id, const std::string &label, unsigned int dim, bool boundary);

    /**
     * As the previous, but generates automatic label of form 'region_ID' if the region with same ID is not already present. Set bulk region.
     * Meant to be used when reading elements from MSH file. Again, if the region is defined already, we just check consistency.
     */
    Region add_region(unsigned int id, unsigned int dim);

    /**
     * Returns a @p Region with given @p label. If it is not found it returns @p undefined Region.
     */
    Region find_label(const std::string &label);

    /**
     * Returns a @p Region with given @p id. If it is not found it returns @p undefined Region.
     */
    Region find_id(unsigned int id);

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
    unsigned int size();
    /**
     * Returns total number boundary regions.
     */
    unsigned int boundary_size();
    /**
     * Returns total number bulk regions.
     */
    unsigned int bulk_size();

    /**
     * Returns list of boundary regions.
     */
    const RegionSet &boundary_regions();

    /**
     * Returns list of boundary regions.
     */
    const RegionSet &bulk_regions();
    /**
     * Returns list of boundary regions.
     */
    const RegionSet &all_regions();

    /*
     * Add region to given set. Create the set if it does not exist.
     *
     * @param set_name Set to which it is added region
     * @param region Added region
     */
    void add_to_set( const string& set_name, RegionIdx region);

    /**
     * Add a set into map, delete possible previous value.
     *
     * @param set_name Name of added set
     * @param set Added RegionSet
     */
    void add_set( const string& set_name, const RegionSet & set);

    RegionSet union_sets( const string & set_name_1, const string & set_name_2); // sort + std::set_union

    RegionSet intersection( const string & set_name_1, const string & set_name_2);

    RegionSet difference( const string & set_name_1, const string & set_name_2);

    /**
     * Get region set of specified name
     *
     * @param set_name Name of set
     * @return RegionSet of specified name
     */
    const RegionSet & get_region_set(const string & set_name);

    void read_sets_from_input(Input::Record rec); // see structure of RegionDB::region_set_input_type


private:
    /// One item in region database
    struct RegionItem {
        RegionItem(unsigned int index, unsigned int id, const std::string &label, unsigned int dim)
            : index(index), id(id), label(label), dim_(dim) {}

        // unique identifiers
        unsigned int index;
        unsigned int id;
        std::string label;
        // data
        unsigned int dim_;
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
    > RegionTable;

    /// Should be RegionSet that consist from all regions. After RegionSets are implemented.
    RegionTable region_set_;

    /// flag for closed database
    bool closed_;
    /// Number of boundary regions
    unsigned int n_boundary_;
    /// Number of bulk regions
    unsigned int n_bulk_;

    /// Map of region sets
    std::map<std::string, RegionSet > sets_;

    /// Make part of general RegionSet table.
    RegionSet all, bulk, boundary;
};



#endif /* REGION_HH_ */
