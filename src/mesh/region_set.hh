/*
 * region_set.hh
 *
 *  Created on: Nov 27, 2012
 *      Author: jb
 */

#ifndef REGION_SET_HH_
#define REGION_SET_HH_

#include "mesh/mesh.h"
#include "mesh/region.hh"
#include "input/accessors.hh"
#include "input/input_type.hh"



/**
 * Base class represented regions.
 *
 * Every descendant must contain:
 *  - constructor what adds region to RegionDB with arguments Input::Record and Mesh
 *  - static generating function created Input Type Record
 *  - static registrar class member for registration class to factory
 */
class RegionSetBase {
public:
	/**
     * Returns whole tree of input types for Region with all descendants.
     */
    static Input::Type::Abstract & get_input_type();

    TYPEDEF_ERR_INFO( EI_Operation_Type, const std::string);
    DECLARE_INPUT_EXCEPTION( ExcEmptyRegionSetResult, << "Empty result of " << EI_Operation_Type::val << " operation.");

protected:
    /// Constructor
    RegionSetBase(Mesh *mesh);
    /// Reference to region_db_ of Mesh
    RegionDB &region_db_;
    /// Reference to map stored relevance of elements to regions.
    RegionDB::MapElementIDToRegionID &el_to_reg_map_;

    unsigned int get_max_region_id() {
    	return region_db_.max_id_+1;
    }
};


/**
 * Region declared by id and name.
 *
 * Allows to create new region with given id and label or specify existing region by id
 * which will be renamed. If existing label is given, it must correspond with appropriate
 * id in RegionDB.
 */
class RegionSetFromId : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
    RegionSetFromId(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

};


/**
 * Region allows to rename existing region specified by mesh_label (e.g. physical volume
 * name in GMSH format).
 */
class RegionSetFromLabel : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
    RegionSetFromLabel(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

};


/**
 * Region declared by name, ID and enum of elements.
 *
 * Allows to get existing region or create new and assign elements to its. Elements are
 * specified by ids. If id of new region is not set, it is generated automatically.
 */
class RegionSetFromElements : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
    RegionSetFromElements(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

};


/**
 * Need new implementation, will be solved later.
 */
/*
class RegionSetBoundary : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
	RegionSetBoundary(const Input::Record &rec, Mesh *mesh);

    /// Returns Input Type Record of Region
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

};
*/


/**
 * Defines region as a union of given two or more other regions.
 *
 * Regions can be given by names or IDs or both ways together. Non-empty set must be the result
 * of the operation.
 */
class RegionSetUnion : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
    RegionSetUnion(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

};


/**
 * Defines region as a difference of given pair of regions.
 *
 * Non-empty set must be the result of the operation.
 */
class RegionSetDifference : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
    RegionSetDifference(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

};


/**
 * Defines region as an intersection of given two or more regions.
 *
 * Non-empty set must be the result of the operation.
 */
class RegionSetIntersection : public RegionSetBase {
public:
    typedef RegionSetBase FactoryBaseType;

	/// Constructor
    RegionSetIntersection(const Input::Record &rec, Mesh *mesh);

	/**
     * Returns Input Type Record of Region
     */
    static const Input::Type::Record & get_region_input_type();

private:
    /// Registrar of class to factory
    static const int registrar;

    /**
     * Get RegionSet of specified name and create its intersection with target RegionSet.
     *
     * @param target_set First RegionSet
     * @param set_name Name of second RegionSet
     * @return RegionSet created of intersection operation
     */
    RegionSet intersection( RegionSet target_set, const string & source_set_name) const;
};



#endif /* REGION_SET_HH_ */
