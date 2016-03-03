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
 * @file    type_generic.hh
 * @brief   
 */

#ifndef TYPE_GENERIC_HH_
#define TYPE_GENERIC_HH_



#include <input/type_base.hh>
#include <input/input_exception.hh>



namespace Input {

namespace Type {

// Declaration of exception
TYPEDEF_ERR_INFO(EI_Object, std::string);
TYPEDEF_ERR_INFO(EI_ParameterList, std::string);
DECLARE_INPUT_EXCEPTION(ExcParamaterNotSubsituted,
        << "No input type substitution for input type parameter " << EI_Object::qval
		<< " found during creation of instance with parameter list: " << EI_ParameterList::val << ".");
DECLARE_INPUT_EXCEPTION(ExcParamaterInIst,
		<< "Parameter " << EI_Object::qval << " appears in the IST. Check where Instance is missing.");
DECLARE_INPUT_EXCEPTION(ExcGenericWithoutInstance,
		<< "Root of generic subtree " << EI_Object::qval << " used without Instance.");



/**
 * @brief Class for representing parametric types in IST.
 *
 * Instances of this class are used only in generic types and during generation
 * of Record are replaced by types of IST (Integer, String, Selection etc.)
 */
class Parameter : public TypeBase {
public:
	/// Constructor.
	Parameter(const string & parameter_name);

	/// Copy constructor.
	Parameter(const Parameter & other);

    /**
     * @brief Implements @p Type::TypeBase::type_name.
     *
     * Parameter type name getter.
     */
    string type_name() const override;
    /// Override @p Type::TypeBase::class_name.
    string class_name() const override { return "Parameter"; }

    /// Implements @p TypeBase::content_hash.
    TypeHash content_hash() const  override;

    /// Implements @p TypeBase::make_instance.
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

    /// Implements @p TypeBase::finish.
    bool finish(bool is_generic = false) override;

protected:
    /// name of parameter
	const string name_;
};



/**
 * @brief Helper class that stores data of generic types.
 */
class Instance : public TypeBase {
public:
	/// Constructor.
	Instance(TypeBase &generic_type, std::vector<TypeBase::ParameterPair> parameters);

	/**
	 * @brief Implements @p TypeBase::content_hash.
	 *
	 * Hash is calculated by hash of generic type and hash of parameters.
	 */
    TypeHash content_hash() const  override;

    /// Used for set Instance to TypeRepository
    const Instance &close() const;

    /// Finish declaration of the Instance type. Call finish of stored @p generic_type_
    bool finish(bool is_generic = false) override;

    /**
     * @brief Implements @p TypeBase::make_instance.
     *
     * In first call creates instance and stores its to @p created_instance_.
     *
     * At each successive call returns this stored type.
     */
    MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

protected:
    /// Reference to generic types (contains some descendants of type @p Parameter).
	TypeBase &generic_type_;

	/// Stores pairs of (name, Input::Type), that are used for replace of parameters in generic types.
	std::vector<TypeBase::ParameterPair> parameters_;

	/**
	 * @brief Stores returned type created in first call of @p make_instance method.
	 *
	 * At each successive call of make_instance returns this stored type.
	 */
	mutable MakeInstanceReturnType created_instance_;
};



} // closing namespace Type
} // closing namespace Input



#endif /* TYPE_GENERIC_HH_ */
