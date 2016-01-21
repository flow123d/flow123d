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
 *  - constructor what adds region to RegionDB
 *  - static generating function created Input Type Record
 */
class RegionSetBase {
public:
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
    RegionSetBase() {}
};


/**
 * Region declared by id and name.
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
 * Region declared by mesh_label and name
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
 * Region declared by name and enum of elements
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
 * Region defined as union of other regions
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
 * Region defined as difference of two other regions
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
 * Region defined as intersection of other regions
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

};



#endif /* REGION_SET_HH_ */
