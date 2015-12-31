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
 * @file    type_abstract.hh
 * @brief
 */

#ifndef TYPE_ABSTRACT_HH_
#define TYPE_ABSTRACT_HH_

#include "system/exceptions.hh"

#include "type_base.hh"
#include "type_selection.hh"


namespace Input {
namespace Type {



/**
 * @brief Class for declaration of polymorphic Record.
 *
 * Like the Record class this is only proxy class. It derives all methods from Record, but
 * further  there is method \p no_more_descendants to close adding descendants. After this
 * you can not derive any Record from this Abstract.
 *
 *
 *
 * A static method (typically part of an abstract class) for construction of an AbstractType can look like:
 *
 @code
    static Input::Type::Abstract &SomeAbstractClass::get_input_type() {
        using namespace Input::Type;
        static Abstract rec("Function",
            "Input for generic time-space function.");

        if (! rec.is_finished()) {
            // declaration of keys that should be derived
            rec.declare_key("domain",Domain::get_input_type(),Default::optional(),
                "Possibly restrict domain of the function.");
            // Finish adding keys.
            rec.finish();

            // call construction of derived Records
            FunctionLinear::get_input_type();
            FunctionPolynomial::get_input_type();
            FunctionInterpreted::get_input_type();

            // finish adding descendants.
            rec.finish();
        }

        return rec;
    }
 @endcode
 *
 * @ingroup input_types
 */
class Abstract : public TypeBase {
	friend class OutputBase;
	//friend class Record;
	friend class AdHocAbstract;

protected:

    /**
     * Actual data of the abstract record.
     */
    class ChildData {
    public:
        ChildData(const string &name, const string &description)
        : selection_of_childs( boost::make_shared<Selection> (name + "_TYPE_selection") ),
		  description_(description),
		  type_name_(name),
		  finished_(false),
		  closed_(false),
		  selection_default_(Default::obligatory())
        {}

        /**
         * Selection composed from names of derived Records. Indices are according to
         * the order of derivation (starting from zero).
         */
        boost::shared_ptr< Selection> selection_of_childs;

        /**
         * Vector of derived Records (proxies) in order of derivation.
         */
        vector< Record > list_of_childs;

        /// Description of the whole Abstract type.
        const string description_;

        /// type_name of the whole Abstract type.
        const string type_name_;

        /// Abstract is finished when it has added all descendant records.
        bool finished_;

        /// If Abstract is closed, we do not allow any further declaration calls.
        bool closed_;

        /**
         * Default value of selection_of_childs (used for automatic conversion).
         *
         * If default value isn't set, selection_default_ is set to obligatory.
         */
        Default selection_default_;

        /**
         * Allow store hash of part of generic subtree.
         *
         * This hash can be typically used if descendants of Abstract contains different
         * structure of parameter location.
         * For example we have have Records with key represents generic part of subtree:
         *  - in first descendant this key is of the type Parameter
         *  - in second descendant this key is of the type Array of Parameter
         *  - in third descendant this key is of the type Array of Parameter with fixed size
         *  etc.
         */
        //TypeHash generic_content_hash_;

    };

public:
    /**
     * Public typedef of constant iterator into array of keys.
     */
    typedef std::vector< Record >::const_iterator ChildDataIter;

    /**
     * Default constructor.
     */
    Abstract();

    /**
     * Copy constructor. We check that other is non empty.
     */
    Abstract(const Abstract& other);


    /**
     * Basic constructor. You has to provide \p type_name of the new declared Record type and
     * its \p description.
     */
    Abstract(const string & type_name_in, const string & description);

    /**
     * Implements @p TypeBase::content_hash.
     *
     * Hash is calculated by type name, description and hash of attributes.
     */
    TypeHash content_hash() const   override;

    /**
     * Allows shorter input of the Abstract providing the default value to the "TYPE" key.
     * If the input reader come across the Abstract in the declaration tree and the input
     * is not 'record-like' with specified value for TYPE, it tries to use the descendant Record specified by
     * @p type_default parameter of this method. Further auto conversion of such Record may be possible.
     */
    Abstract &allow_auto_conversion(const string &type_default);

    /**
     *  Can be used to close the Abstract for further declarations of keys.
     */
    Abstract &close();

    /**
     *  Finish declaration of the Abstract type.
     *
     *  Set Abstract as parent of derived Records (for mechanism of
     *  set parent and descendant see \p Record::derive_from)
     */
    bool finish(bool is_generic = false) override;

    /**
     * Returns reference to the inherited Record with given name.
     */
    const Record  &get_descendant(const string& name) const;

    /**
     * Returns default descendant if TYPE key has default value, otherwise returns empty Record.
     */
    const Record * get_default_descendant() const;

    /**
     * Returns reference to Selection type of the implicit key TYPE.
     */
    const Selection &get_type_selection() const;

    /**
     * Returns number of keys in the child_data_.
     */
    unsigned int child_size() const;

    /**
     * Implements @p TypeBase::is_finished.
     */
    virtual bool is_finished() const override;

    /// Returns true if @p data_ is closed.
    virtual bool is_closed() const override;

    /// Abstract type name getter.
    virtual string type_name() const override;
    string class_name() const override { return "Abstract"; }

    /**
     * Container-like access to the data of the Record. Returns iterator to the first data.
     */
    ChildDataIter begin_child_data() const;

    /**
     * Container-like access to the data of the Record. Returns iterator to the last data.
     */
    ChildDataIter end_child_data() const;

    /**
     * Add inherited Record. This method is used primarily in combination with registration
     * variable. @see Input::Factory
     *
     * Example of usage:
	 @code
		 class SomeBase
		 {
		 public:
    		/// the specification of input abstract record
    		static const Input::Type::Abstract & get_input_type();
			...
		 }

		 class SomeDescendant : public SomeBase
		 {
		 public:
    		/// the specification of input record
    		static const Input::Type::Record & get_input_type();
			...
		 private:
			/// registers class to factory
			static const int reg;
		 }

		 /// implementation of registration variable
		 const int SomeDescendant::reg =
			 Input::register_class< SomeDescendant >("SomeDescendant") +
			 SomeBase::get_input_type().add_child(SomeDescendant::get_input_type());
	 @endcode
     */
    int add_child(const Record &subrec);

    // Get default value of selection_of_childs
    Default &get_selection_default() const;

    // Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

    /// Create deep copy of Abstract (copy all data stored in shared pointers etc.)
    Abstract deep_copy() const;

    /// Set flag @p root_of_generic_subtree_ to true
    Abstract &root_of_generic_subtree();

    /// Set @p generic_content_hash_
    //Abstract &set_generic_content_hash(TypeHash generic_content_hash);

protected:
    /**
     * This method intentionally have no implementation to
     * prevents deriving an Abstract form other Abstract.
     * In such a case the linker should report an undefined reference.
     */
    Record &derive_from(Abstract &parent);

    /**
     * Check if type has set value of default descendants.
     */
    bool have_default_descendant() const;

    /// Actual data of the Abstract.
    boost::shared_ptr<ChildData> child_data_;

    friend class Record;
};


/** ******************************************************************************************************************************
 * Class for declaration of polymorphic Record.
 *
 * Abstract extends on list of descendants provided immediately
 * after construction by add_child(). These descendants derive
 * only keys from common AR. AdHocAR has separate instance for every
 * key of this type.
 *
 * @ingroup input_types
 */
class AdHocAbstract : public Abstract {
	friend class OutputBase;
public:
	/**
	 * Constructor
	 */
	AdHocAbstract(const Abstract &ancestor);

	TypeHash content_hash() const   override;

	string class_name() const override { return "AdHocAbstract"; }


    /**
     * Finish declaration of the AdHocAbstract type. Adds descendants of ancestor Abstract,
     * calls close() and complete keys with non-null pointers to lazy types.
     */
    bool finish(bool is_generic = false) override;

    /**
     * Add inherited Record.
     */
    AdHocAbstract &add_child(const Record &subrec);

protected:
    /// Reference to ancestor Abstract
    const Abstract &ancestor_;
};



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_ABSTRACT_HH_ */
