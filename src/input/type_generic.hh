/*
 * type_generic.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_GENERIC_HH_
#define TYPE_GENERIC_HH_



#include <input/type_base.hh>
#include <input/accessors.hh>



namespace Input {

namespace Type {

// Declaration of exception
TYPEDEF_ERR_INFO(EI_Object, std::string);
TYPEDEF_ERR_INFO(EI_ParameterList, std::string);
DECLARE_INPUT_EXCEPTION(ExcParamaterNotSubsituted,
        << "No input type substitution for input type parameter " << EI_Object::qval
		<< " found during creation of instance with parameter list: " << EI_ParameterList::val << ".");



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

    /// Implements @p TypeBase::content_hash.
    TypeHash content_hash() const  override;

    /// Implements @p TypeBase::valid_default.
    virtual bool valid_default(const string &str) const override;

    /// Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

    /// Implements @p TypeBase::finish.
    bool finish(bool is_generic = false) override;

protected:
    /// name of parameter
	const string name_;
};



/**
 * Helper class that stores data of generic types.
 */
class Instance : public TypeBase {
public:
	Instance(TypeBase &generic_type, std::vector<TypeBase::ParameterPair> parameters);

	/// Implements @p TypeBase::content_hash.
    TypeHash content_hash() const  override;

    /// Implements @p TypeBase::valid_default.
    virtual bool valid_default(const string &str) const override;

    /// Used for set Instance to TypeRepository
    const Instance &close() const;

    /// Finish declaration of the Instance type. Call finish of stored @p generic_type_
    bool finish(bool is_generic = false) override;

    // Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

protected:
    /// Reference to generic types (contains some descendants of type @p Parameter).
	TypeBase &generic_type_;

	/// Stores pairs of (name, Input::Type), that are used for replace of parameters in generic types.
	std::vector<TypeBase::ParameterPair> parameters_;
};



} // closing namespace Type
} // closing namespace Input



#endif /* TYPE_GENERIC_HH_ */
