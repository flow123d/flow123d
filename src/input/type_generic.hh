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
    virtual string type_name() const override;

    TypeHash content_hash() const  override;

    virtual bool valid_default(const string &str) const override;

    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

protected:
	const string name_;
};



/**
 * Helper class that stores data of generic types.
 */
class Instance : public TypeBase {
public:
	Instance(TypeBase &generic_type, std::vector<TypeBase::ParameterPair> parameters);

    TypeHash content_hash() const  override;

    virtual bool valid_default(const string &str) const override;

    /// Used for set Instance to TypeRepository
    const Instance &close() const;

    /// Finish declaration of the Instance type. Call finish of stored @p generic_type_
    bool finish(bool is_generic = false) override;

    // Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

protected:
	TypeBase &generic_type_;

	std::vector<TypeBase::ParameterPair> parameters_;
};



} // closing namespace Type
} // closing namespace Input



#endif /* TYPE_GENERIC_HH_ */
