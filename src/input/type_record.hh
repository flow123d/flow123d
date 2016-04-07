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
 * @file    type_record.hh
 * @brief   
 */

#ifndef TYPE_RECORD_HH_
#define TYPE_RECORD_HH_

#include "system/exceptions.hh"

#include "type_base.hh"
#include "storage.hh"


namespace Input {
namespace Type {


/** *********************************************************************************************************************************
 * @brief Class \p Input::Type::Default specifies default value of keys of a \p Input::Type::Record.
 *
 * It contains type of default value and possibly the value itself stored as \p std::string. Currently we distinguish four
 * cases:
 * - \b Default value given \b at \b declaration time, i.e. the default value is part of the \p Input::Type specification.
 *   This should be preferred way to give the default value since it can by documented as part of
 *   Record type specification.
 * - \b Default value given \b at \b read time, i.e. when you ask for the value through the Input::Record accessor. It should be used only
 *   if the default value is not constant, i.e. is taken from a global setting. In this case you should provide textual description
 *   where the default value comes from. Of course it could be difficult if the value is read on several places with different default values.
 * - \b No \b default value is given and input value is \b obligatory. An exception is thrown when the key is missing on the input.
 * - \b No \b default value is given and input value is \b optional.
 *
 *
 *
 * @ingroup input_types
 */
class Default {
	friend class Record;
	friend class OutputBase;
	//friend class OutputJSONTemplate;

private:
    /**
     * Possible types of default values.
     */
    enum DefaultType {
        default_at_declaration,        ///< Default value given at declaration time.
        default_at_read_time,          ///< Some default value will be given when the key is read. The description of this value should be provided.
        no_default_optional_type,      ///< No default value, optional key. This is default type of the Default.
        no_default_obligatory_type     ///< No default value, obligatory key.
    };
public:

    /**
     * Constructor with given default value (at declaration time)
     */
    Default(const std::string & value);

    /**
     * Hash of the Default specification, counted of type_ and value_.
     */
    TypeBase::TypeHash content_hash() const;

    /**
     * Factory function to make an default value that will be specified at the time when a key will be read.
     * You have to provide a string with description of the default value used at the read time., e.g.
     * the key \p time_governer of an equation can specify default value as
     * @code
     *      Default::read_time("By default the global time governor is used.")
     * @endcode
     * To get the value of such key from the input you have to use non-throwing variant of the method
     * Input::Record::key, which returns the value through reference and allows checking presence of the key on the input.
     *
     * Example of usage:
     * @code
     *      some_record.declare_key("time_governor",TimeGovernor(),Default::optional(),"description");
     * @endcode
     */
    static Default read_time(const std::string & description)
    { return Default(default_at_read_time, description ); }

    /**
     * Factory function to make an empty default value which is obligatory.
     * This and following factory functions should be used instead of private constructors.
     *
     * Example of usage:
     * @code
     *      some_record.declare_key("some_key",Integer(),Default::obligatory(),"description");
     * @endcode
     */
    inline static Default obligatory()
    { return Default(no_default_obligatory_type, "OBLIGATORY"); }

    /**
     * Factory function to make an empty default value which is optional.
     * To get the value of such key from the input you have to use non-throwing variant of the method
     * Input::Record::key, which returns the value through reference and allows checking presence of the key on the input.
     *
     * Example of usage:
     * @code
     *      some_record.declare_key("some_key",Integer(),Default::optional(),"description");
     * @endcode
     */
    inline static Default optional()
    { return Default(no_default_optional_type, "OPTIONAL"); }

    /**
     * Returns true if the default value is or will be available when someone tries to read the value.
     */
    inline bool has_value_at_read_time() const
    { return (type_ == default_at_read_time); }

    /**
     * Returns true if the default value is or will be available when someone tries to read the value.
     */
    inline bool has_value_at_declaration() const
    { return (type_ == default_at_declaration); }


    /**
     * Returns true if the key is obligatory and thus must be specified on input. No default value is given.
     */
    inline bool is_obligatory() const
    { return (type_ == no_default_obligatory_type); }

    /**
     * Returns true if the key is optional.
     */
    inline bool is_optional() const
    { return (type_ == no_default_optional_type); }

    /**
     * Returns stored value. Possibly empty string.
     */
    inline const string & value() const
    { return (value_); }

    /**
     * Compares values type_ of two Default objects.
     */
    inline bool has_same_type(const Default &other) const
        {return type_ == other.type_; }

    /**
     * Check validity of @p value_ using the JSON reader
     * if default type is default_at_declaration.
     */
    bool check_validity(boost::shared_ptr<TypeBase> type) const;

    /**
     * Return @p storage_, if storage_ is NULL, call check_validity method
     */
    Input::StorageBase *get_storage(boost::shared_ptr<TypeBase> type) const;

private:
    string value_;                          ///< Stored value.
    enum DefaultType type_;                 ///< Type of the Default.
    mutable Input::StorageBase *storage_;   ///< Storage of default value read by reader

    /**
     * Forbids default constructor.
     */
    Default();

    /**
     * Constructor for other types then 'declaration'.
     */
    Default(enum DefaultType type, const std::string &value);
};


class Abstract;


/** ******************************************************************************************************************************
 * @brief Record type proxy class.
 *
 * To keep consistency, we have to prevent copies of the actual Record data. Therefore this class is just a proxy  that
 * can be freely (and cheaply) copied.
 *
 * @ingroup input_types
 */
class Record : public TypeBase {
	friend class OutputBase;
	friend class Abstract;
	friend class AdHocAbstract;

public:

    /*
     * Exceptions specific to this class.
     */
    TYPEDEF_ERR_INFO( EI_Record, Record );
    TYPEDEF_ERR_INFO( EI_RecordName, const string);
    DECLARE_EXCEPTION( ExcRecordKeyNotFound, << "Key " << EI_KeyName::qval <<" not found in Record:\n" <<  EI_Record::val );

    /**
     *  Structure for description of one key in record.
     *  The members dflt_type_ and default have reasonable meaning only for
     *  type_ == Scalar
     */
    struct Key {
        unsigned int key_index;                     ///< Position inside the record.
        string key_;                                ///< Key identifier.
        string description_;                        ///< Key description in context of particular Record type.
        boost::shared_ptr<TypeBase> type_;          ///< Type of the key.
        Default default_;                      ///< Default, type and possibly value itself.
        bool derived;                               ///< Is true if the key was only derived from the parent Record, but not explicitly declared.
    };

    /**
     * Public typedef of constant iterator into array of keys.
     */
    typedef std::vector<struct Key>::const_iterator KeyIter;

    /**
     * Default constructor. Empty handle.
     */
    Record();

    /**
     * Copy constructor. We allow only copies of non-empty records.
     */
    Record(const Record & other);


    /**
     * Basic constructor. You have to provide \p type_name of the new declared Record type and
     * its \p description.
     */
    Record(const string & type_name_in, const string & description);


    /**
     * Implements @p TypeBase::content_hash.
     *
     * Hash is calculated by type name, description, auto conversion key, hash of keys and attributes.
     */
    virtual TypeHash content_hash() const  override;


    /**
     * Method to derive new Record from an AbstractRecord @p parent. This register the @p parent in the newly
     * created Record. Method creates TYPE key of Record and must be call before declaration of keys.
     *
     * Mechanism of set parent to derived Record and child to parent Abstract is a bit more complicated. For
     * correct finish it must be done in these steps:
     *
     * - in derive_from is set @p parent to derived Record
     * - in \p close is set derived Record to parent (or to all parents for multiple inheritance) and registered
     *   parents in derived Record are erased
     * - in \p AbstractRecord::finish is re-registered parents to descendant (through \p add_parent method)
     */
    virtual Record &derive_from(Abstract &parent);

    /**
     * Copy keys from other record. If @p other record is not yet constructed, we postpone copy to the finish phase.
     */
    Record &copy_keys(const Record &other);

    /**
     * Allows shorter input of the Record providing only value of the \p from_key given as the parameter.
     * All other keys of the Record must have default values specified at declaration. This is checked when the
     * \p finish method is called.
     *
     * If the input reader come across the Record in declaration tree, but there is not 'record-like' input, it
     * save default values into storage tree and tries to match the input with the type of the \p from_key.
     */
    virtual Record &allow_auto_conversion(const string &from_key);

    /**
     * Declares a key of the Record with name given by parameter @p key, the type given by target of pointer @p type,
     * default value by parameter @p default_value, and with given @p description.
     * The parameter @p type points to a descendant of TypeBase.
     */
    Record &declare_key(const string &key, boost::shared_ptr<TypeBase> type,
                            const Default &default_value, const string &description);

    /**
     * Declares a key of the Record with name given by parameter @p key, the type given by parameter @p type,
     * default value by parameter @p default_value, and with given @p description.
     * The parameter @p type has a descendant of TypeBase.
     */
    template <class KeyType>
    Record &declare_key(const string &key, const KeyType &type,
                            const Default &default_value, const string &description);


    /**
     * Same as previous method but without given default value (same as Default() - optional key )
     */
    template <class KeyType>
    Record &declare_key(const string &key, const KeyType &type,
                            const string &description);


    /**
     * Close the Record for further declarations of keys.
     *
     * Add Record to type repository (see @p TypeRepository::add_type) and set Record
     * as descendant of parent if Record is derived (for mechanism of set parent and
     * descendant see \p derive_from)
     */
    const Record &close() const;


    /**
     * Implements @p TypeBase::is_finished.
     */
    bool is_finished() const override;

    /// Returns true if @p data_ is closed.
    bool is_closed() const override;

    /// Record type name getter.
    string type_name() const override;
    virtual string class_name() const override { return "Record"; }

    /// Class comparison and Record type name comparision.
    bool operator==(const TypeBase &other) const;

    /**
     * Interface to mapping key -> index in record. Returns index (in continuous array) for given key.
     *
     * Works also for unfinished Record.
     */
    inline unsigned int key_index(const string& key) const;

    /**
     * Returns iterator to the key struct for given key string.
     */
    inline KeyIter key_iterator(const string& key) const;

    /**
     * Returns iterator to auto-conversion key (see Record::allow_auto_conversion), or end() if the auto conversion
     * is not allowed.
     */
    KeyIter auto_conversion_key_iter() const;

    /**
     * Returns iterator to the key struct for given key string.
     *
     */
    inline bool has_key_iterator(const string& key, KeyIter &it) const;

    /**
     * Container-like access to the keys of the Record. Returns iterator to the first key.
     */
    inline KeyIter begin() const;

    /**
     * Container-like access to the keys of the Record. Returns iterator to the last key.
     */
    inline KeyIter end() const;

    /**
     * Returns true if the Record contains key with given name.
     */
    inline bool has_key(const string& key) const;

    /**
     * Returns number of keys in the Record.
     */
    inline unsigned int size() const;

    /**
     * Finish declaration of the Record type. Calls close() and complete keys with non-null pointers to lazy types.
     */
    bool finish(bool is_generic = false) override;

    /**
     * Add TYPE key as obligatory.
     *
     * This method can't be used for derived record.
     */
    Record &has_obligatory_type_key();

    /// Implements @p TypeBase::make_instance.
    virtual MakeInstanceReturnType make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const override;

    /// Create deep copy of Record (copy all data stored in shared pointers etc.)
    Record deep_copy() const;

    /// Set flag @p root_of_generic_subtree_ to true
    virtual Record &root_of_generic_subtree();


protected:

    /**
     * Assertion for finished Type::Record.
     */
    inline void finished_check() const {
        ASSERT( is_finished(), "Asking for information of unfinished Record type: %s\n", type_name().c_str());
    }

    /// Auxiliary method that actually makes the copy of keys.
    void make_copy_keys(Record &origin);

    /**
     * Declares a TYPE key of the Record.
     */
    Record &declare_type_key();

    /**
     * Set parent Abstract of Record.
     *
     * This method is created for correct functionality of generic types. It must be called
     * in Abstract::finish() and refill @p parent_vec_ vector of correct parents (for complete
     * mechanism of set parent and descendant see \p derive_from)
     */
    const Record &add_parent(Abstract &parent) const;

    /**
     * @brief Set data of Instance of generic type.
     *
     * Called from make_instance method and set data of Record or its descendants.
     */
    void set_instance_data(Record &rec, ParameterMap &parameter_map, std::vector<ParameterPair> vec) const;

    /**
     * Internal data class.
     */
    class RecordData  {
    public:
        /// Constructor
        RecordData(const string & type_name_in, const string & description);

        /**
         * Declares a key and stores its type. The type parameter has to be finished at the call of declare_key().
         * If the parameter @p type_temporary is NULL, the parameter @p type provides pointer to
         * already finished type that will be assigned to the key. On the other hand, if @p type_temporary is not NULL,
         * only this raw pointer is stored and key is fully completed later through TypeBase::lazy_finish().
         */
        void declare_key(const string &key,
                         boost::shared_ptr<TypeBase> type,
                         const Default &default_value, const string &description);


        Record::KeyIter auto_conversion_key_iter() const;

        /// Count hash of RecordData.
        void content_hash(TypeBase::TypeHash &seed) const;

        /// Database of valid keys
        std::map<KeyHash, unsigned int> key_to_index;
        typedef std::map<KeyHash, unsigned int>::const_iterator key_to_index_const_iter;

        /// Keys in order as they where declared.
        std::vector<struct Key> keys;

        /// Description of the whole record type.
        const string description_;
        const string type_name_;

        /// Permanent pointer to parent Abstract, necessary for output.
        std::vector< boost::shared_ptr<Abstract> > parent_vec_;

        /// Record is finished when it is correctly derived (optional) and have correct shared pointers to types in all keys.
        bool finished;

        /// If record is closed, we do not allow any further declare_key calls.
        bool closed_;

        /// True for derived records after make_derived.
        bool derived_;

        /**
         * Initial value is = -1, when allow_auto_conversion is called we set this to 0.
         * Final value can be assigned just after possible inheritance copy of keys from parent Abstract.
         */
        int auto_conversion_key_idx;
        /**
         * Name of key to use for auto conversion.
         */
        std::string auto_conversion_key;

    };

    /// Data handle.
    boost::shared_ptr<RecordData> data_;
};


/*********************************************************
 * Implementation
 */




inline unsigned int Record::key_index(const string& key) const
{
    KeyHash key_h = key_hash(key);
    RecordData::key_to_index_const_iter it = data_->key_to_index.find(key_h);
    if (it != data_->key_to_index.end()) return it->second;
    else
        THROW( ExcRecordKeyNotFound() << EI_KeyName(key) << EI_Record(*this) );

    return size();
}



inline Record::KeyIter Record::key_iterator(const string& key) const
{

    finished_check();
    return begin() + key_index(key);
}



inline bool Record::has_key_iterator(const string& key, KeyIter &it) const
{
    finished_check();
    KeyHash key_h = key_hash(key);
    RecordData::key_to_index_const_iter data_it = data_->key_to_index.find(key_h);
    if (data_it == data_->key_to_index.end()) {
        return false;
    } else {
        it = begin()+data_it->second;
        return true;
    }
}



inline Record::KeyIter Record::begin() const
{
    finished_check();
    return data_->keys.begin();
}



inline Record::KeyIter Record::end() const
{
    finished_check();
    return data_->keys.end();
}



inline bool Record::has_key(const string& key) const
{
    return key_iterator(key) != end();
}



inline unsigned int Record::size() const {
	ASSERT( is_closed(), "Asking for information of unclosed Record type: %s\n", type_name().c_str());
    ASSERT_EQUAL( data_->keys.size(), data_->key_to_index.size());
    return data_->keys.size();
}






} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_RECORD_HH_ */
