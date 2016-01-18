/*
 * region_types.hh
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#ifndef REGION_TYPES_HH_
#define REGION_TYPES_HH_

#include "mesh/mesh.h"
#include "mesh/region.hh"
#include "input/accessors.hh"
#include "input/input_type.hh"



/**
 * Base class represented regions.
 *
 * Every descendant must contain:
 *  - constructor what adds region to RegionDB
 *  - static generating function created Input Type Record
 */
class RegionBase : public Region {
public:
	/// Constructor
	RegionBase(unsigned int index, const RegionDB &db)
	: Region(index, db) {}

	/**
     * Returns whole tree of input types for Region with all descendants.
     */
    static Input::Type::Abstract & get_input_type();

    TYPEDEF_ERR_INFO( EI_Region_Label, const std::string);
    DECLARE_EXCEPTION( ExcNonexistingLabel, << "Non-existing label of region: " << EI_Region_Label::qval << "\n" \
                                             << "You must also set ID or use existing label.\n");

    /// Call appropriate add_region methods of RegionDB
    Region add_region(Mesh *mesh, unsigned int id, const std::string &label);
    Region add_region(Mesh *mesh, unsigned int id, const std::string &label, unsigned int dim);
    Region add_region(Mesh *mesh, unsigned int id, unsigned int dim);
    void add_set(Mesh *mesh, const string& set_name, const RegionSet & set);

protected:
    /// Empty constructor
	RegionBase()
	: Region() {}
};


/**
 * Region declared by id and name.
 */
class RegionFromId : public RegionBase {
public:
	/// Constructor
	RegionFromId(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

};


/**
 * Region declared by mesh_label and name
 */
class RegionFromLabel : public RegionBase {
public:
	/// Constructor
	RegionFromLabel(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

};


/**
 * Region declared by name and enum of elements
 */
class RegionFromElements : public RegionBase {
public:
	/// Constructor
	RegionFromElements(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

};


/**
 * Need new implementation, will be solved later.
 */
/*
class RegionBoundary : public RegionBase {
public:
	/// Constructor
	RegionBoundary(const Input::Record &rec, Mesh *mesh);

    /// Returns Input Type Record of Region
    static const Input::Type::Record & get_region_input_type();

};
*/


/**
 * Region defined as union of other regions
 */
class RegionUnion : public RegionBase {
public:
	/// Constructor
	RegionUnion(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

};


/**
 * Region defined as difference of two other regions
 */
class RegionDifference : public RegionBase {
public:
	/// Constructor
	RegionDifference(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

};


/**
 * Region defined as intersection of other regions
 */
class RegionIntersection : public RegionBase {
public:
	/// Constructor
	RegionIntersection(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

};



#endif /* REGION_TYPES_HH_ */
