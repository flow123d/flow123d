/*
 * region.hh
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 *
 *  TODO:
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
        ASSERT( is_boundary(), "Try to get boundary index of a bulk region with internal idx: %d\n", idx_ );
        return idx_ >> 1; }

    /// Returns index of the region in the bulk set.
    inline unsigned int bulk_idx() const {
        ASSERT( ! is_boundary(), "Try to get bulk index of boundary region with internal idx: %d\n", idx_ );
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

private:
    /**
     * Create accessor from the index. Should be private since implementation specific.
     * We need some way how to iterate over: all regions, boundary regions, bulk regions -
     * solution: have specific RegionSets for these three cases.
     */
    Region(unsigned int index, const RegionDB &db)
    : RegionIdx(index), db_(&db)
    {}

    /// Comparative method of two regions
    static bool comp(const Region &a, const Region &b)
    { return a.idx_ < b.idx_; }

    /// Global variable with information about all regions.
    const RegionDB *db_;

    friend class RegionDB;
};




/**
 * Type representing a set of regions.
 * CAn be used  to set function(field) on more regions at once, possibly across meshes
 *
 * Regions stored in region set are always unique
 */
typedef std::vector<Region> RegionSet;



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
 *
 * typedef std::map<unsigned int, unsigned int> MapElementIDToRegionID;
 * RegionDB::read_regions_from_input(Input::Array region_list, MapElementIDToRegionID &map);
 *
 *
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
    /**
     * Map representing the relevance of elements to regions
     */
    typedef std::map<unsigned int, unsigned int> MapElementIDToRegionID;

    /**
     * Format of input record which defined elements and their affiliation to region sets
     */
    static const Input::Type::Record & get_region_input_type();
    /**
     * Format of input record which defined regions and their affiliation to region sets
     */
    static const Input::Type::Record & get_region_set_input_type();

    TYPEDEF_ERR_INFO( EI_Label, const std::string);
    TYPEDEF_ERR_INFO( EI_ID, unsigned int);
    TYPEDEF_ERR_INFO( EI_IDOfOtherLabel, unsigned int);
    TYPEDEF_ERR_INFO( EI_LabelOfOtherID, const std::string);
    DECLARE_EXCEPTION( ExcAddingIntoClosed, << "Can not add label=" << EI_Label::qval << " into closed MaterialDispatch.\n");
    DECLARE_EXCEPTION( ExcNonuniqueID, << "Non-unique ID during add of region id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "other region with same ID but different label: " << EI_LabelOfOtherID::qval << " already exists\n");
    DECLARE_EXCEPTION( ExcNonuniqueLabel, << "Non-unique label during add of region id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "other region with same label but different ID: " << EI_IDOfOtherLabel::val << " already exists\n");
    DECLARE_EXCEPTION( ExcInconsistentBoundary, << "Inconsistent add of region with id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "both ID and label match an existing region with different boundary flag.");
    DECLARE_EXCEPTION( ExcInconsistentDimension, << "Inconsistent add of region with id: " << EI_ID::val << ", label: " << EI_Label::qval << "\n" \
                                             << "both ID and label match an existing region with different dimension.");

    DECLARE_EXCEPTION( ExcCantAdd, << "Can not add new region into DB, id: " << EI_ID::val <<", label: " << EI_Label::qval);

    DECLARE_EXCEPTION( ExcUnknownSet, << "Operation with unknown region set: " << EI_Label::qval );

    DECLARE_INPUT_EXCEPTION( ExcUnknownSetOperand, << "Operation with unknown region set: " << EI_Label::qval);

    TYPEDEF_ERR_INFO( EI_NumOp, unsigned int);
    DECLARE_INPUT_EXCEPTION( ExcWrongOpNumber, << "Wrong number of operands. Expect 2, given: " << EI_NumOp::val);

    DECLARE_INPUT_EXCEPTION(ExcUniqueRegionId, << "Id of region must be unique, id: " << EI_ID::val );

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
     * @p dim is dimension of reference elements in the region and @p boundary is true if the region consist of boundary elements
     * (where one can apply boundary condition).
     *
     */
    Region add_region(unsigned int id, const std::string &label, unsigned int dim);

    /**
     * As the previous, but set the 'boundary; flag according to the label (labels starting with dot '.' are boundary).
     * Used in read_regions_from_input ( with undefined dimension) to read regions given in 'regions' key of the 'mesh' input record.
     */
    Region add_region(unsigned int id, const std::string &label);

    /**
     * As the previous, but generates automatic label of form 'region_ID' if the region with same ID is not already present. Set bulk region.
     * Meant to be used when reading elements from MSH file. Again, if the region is defined already, we just check consistency.
     */
    Region add_region(unsigned int id, unsigned int dim);

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
     * Get RegionSets of specified names and create their union
     *
     * @param set_name_1 Name of first RegionSet
     * @param set_name_2 Name of second RegionSet
     * @return RegionSet created of union operation
     */
    RegionSet union_sets( const string & set_name_1, const string & set_name_2);

    /**
     * Get RegionSets of specified names and create their intersection.
     * Throws ExcUnknownSet for invalid name.
     *
     * @param set_name_1 Name of first RegionSet
     * @param set_name_2 Name of second RegionSet
     * @return RegionSet created of intersection operation
     */
    RegionSet intersection( const string & set_name_1, const string & set_name_2);

    /**
     * Get RegionSets of specified names and create their difference
     * Throws ExcUnknownSet for invalid name.
     *
     * @param set_name_1 Name of first RegionSet
     * @param set_name_2 Name of second RegionSet
     * @return RegionSet created of difference operation
     */
    RegionSet difference( const string & set_name_1, const string & set_name_2);

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
     * Reads region sets defined by user in input file
     * Format of input record is defined in variable RegionDB::get_region_set_input_type()
     *
     * @param arr Array input records which define region sets
     */
    void read_sets_from_input(Input::Array arr);

    /**
     * Reads elements and their affiliation to region sets defined by user in input file
     * Format of input record is defined in method RegionDB::get_region_input_type()
     *
     * @param region_list Array input records which define region sets and elements
     * @param map Map to which is loaded data
     */
    void read_regions_from_input(Input::Array region_list, MapElementIDToRegionID &map);


private:


    typedef std::pair<unsigned int, unsigned int> DimID;

    /// One item in region database
    struct RegionItem {
        RegionItem(unsigned int index, unsigned int id, const std::string &label, unsigned int dim)
            : index(index), id(dim, id), label(label) {}

        unsigned int get_id() const {return id.second;}
        unsigned int dim() const {return id.first;}

        // unique identifiers
        unsigned int index;
        DimID id;
        std::string label;
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
    RegionTable region_set_;

    /// flag for closed database, no regions can be added, but you can add region sets
    bool closed_;
    /// Number of boundary regions
    unsigned int n_boundary_;
    /// Number of bulk regions
    unsigned int n_bulk_;

    /// Map of region sets
    std::map<std::string, RegionSet > sets_;

    /// Make part of general RegionSet table.
    RegionSet all, bulk, boundary;

    /**
     * Implicit bulk and boundary regions. For GMSH mesh format only implicit_boundary is used for boundary elements
     * that are not explicitly in the mesh file.
     */
    Region implicit_bulk_, implicit_boundary_;

    /**
     * Prepare region sets for union, intersection and difference operation.
     * Get sets of names set_name_1 and set_name_2 and sort them.
     * Throws ExcUnknownSet if the set with given name does not exist.
     */
    void prepare_sets( const string & set_name_1, const string & set_name_2, RegionSet & set_1, RegionSet & set_2);

    /**
     * Read two operands from input array of strings and check if given names
     * are existing sets. Return pair of checked set names.
     */
    pair<string,string> get_and_check_operands(const Input::Array & operands);

    /**
     * Create label of region in format: "region_"+id
     *
     * Use if label is not set.
     */
    void create_label_from_id(const string & label, unsigned int id);

    /**
     * Insert new region into database.
     */
    Region insert_region(unsigned int id, const std::string &label, unsigned int dim, bool boundary);

    /**
     * Replace dimension of existing region with undefined_dim.
     */
    Region replace_region_dim(DimIDIter it_undef_dim, unsigned int dim, bool boundary);

    /**
     * Find existing region given by pair (dim, id).
     */
    Region find_by_dimid(DimIDIter it_id, unsigned int id, const std::string &label, bool boundary);

    /**
     * Return boundary flag for given label. Label of boundary region must start by '.' symbol.
     */
    inline bool is_boundary(const std::string &label) {
    	return (label.size() != 0) && (label[0] == '.');
    }

};



#endif /* REGION_HH_ */
