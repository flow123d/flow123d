/*
 * type_generic.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_GENERIC_HH_
#define TYPE_GENERIC_HH_



#include <input/type_base.hh>



namespace Input {

namespace Type {


/**
 * Class for representing parametric types in IST.
 *
 * Instances of this class are used only in generic types and during generation
 * of Record are replaced by types of IST (Integer, String, Selection etc.)
 */
class Parameter : public TypeBase {
public:
	Parameter(const string & parameter_name);

	Parameter(const Parameter & other);

    /// Parameter type name getter.
    virtual string type_name() const;

    TypeHash content_hash() const  override;

    virtual bool valid_default(const string &str) const;

protected:
	const string name_;
};



/**
 * Helper class that stores data of generic types.
 */
class Instance : public TypeBase {
public:
	typedef std::pair< std::string, boost::shared_ptr<const TypeBase> > ParameterPair;

	Instance(const TypeBase &generic_type, std::vector<ParameterPair> parameters);

    TypeHash content_hash() const  override;

    virtual bool valid_default(const string &str) const;

protected:
	const TypeBase &generic_type_;

	std::vector<ParameterPair> parameters_;
};



} // closing namespace Type
} // closing namespace Input



#endif /* TYPE_GENERIC_HH_ */
