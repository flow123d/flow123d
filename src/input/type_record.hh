/*
 * type_record.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_RECORD_HH_
#define TYPE_RECORD_HH_

#include "type_base.hh"
#include "type_selection.hh"

/**
 * Macro to create a key object. The reason for this is twofold:
 * 1) We may use it to implement compile time computed hashes for faster access to the data.
 * 2) We can store line, function, filename where the key where used in order to report more specific error messages.
 */
#define KEY(name) #name


namespace Input {
namespace Type {


class AbstractRecord;


/**
 * @brief Record type proxy class.
 *
 * To keep consistency, we have to prevent copies of the actual Record data. Therefore this class is just a proxy  that
 * can be freely (and cheaply) copied.
 *
 */
class Record : public TypeBase {
public:
    /**
     *  Structure for description of one key in record.
     *  The members dflt_type_ and default have reasonable meaning only for
     *  type_ == Scalar
     */
    struct Key {
        unsigned int key_index;                     ///< Position inside the record.
        string key_;                                ///< Key identifier.
        string description_;                        ///< Key description in context of particular Record type.
        boost::shared_ptr<const TypeBase> type_;    ///< Type of the key.
        DefaultValue default_;                      ///< DefaultValue, type and possibly value itself.
        bool derived;                               ///< Is true if the key was only derived from the parent Record, but not explicitly declared.
    };

    /**
     * Public typedef of constant iterator into array of keys.
     */
    typedef std::vector<struct Key>::const_iterator KeyIter;

    /**
     * Default constructor. Empty handle.
     */
    Record()
    { finished=false; }

    /**
     * Basic constructor. You has to provide @t type_name of the new declared Record type and
     * its @t description.
     */
    Record(const string & type_name_in, const string & description);

    /**
     * Constructor to derive new Record from an AbstractRecord @p parent. This copy all keys from the @p parent and register the newly created Record
     * in the @p parent. You are free to overwrite copied keys, but you can not delete them.
     */
    Record(AbstractRecord parent, const string & type_name_in, const string & description);

    /**
     * Declares a key of the Record with name given by parameter @p key, the type given by parameter @p type, default value by parameter @p default_value, and with given
     * @p description. The parameter @p type has to be any of descendants of TypeBase.
     *
     */
    template <class KeyType>
    void declare_key(const string &key,
                            const KeyType &type,
                            const DefaultValue &default_value, const string &description)
    {
        // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );

        if (data_.use_count() == 0) xprintf(PrgErr, "Can not declare a new key `%s` of an empty Record proxy.\n", key.c_str());
        if (finished) xprintf(PrgErr, "Declaration of key: %s in finished Record type: %s\n", key.c_str(), type_name().c_str());

        // If KeyType is not derived from Scalar, we check emptiness of the default value.
        if (boost::is_base_of<Scalar, KeyType>::value == false && default_value.has_value() )
            xprintf(Err, "Default value for non scalar type in declaration of key: %s in Record type: %s \n", key.c_str(), type_name().c_str() );

        if (! is_valid_identifier(key))
            xprintf(PrgErr, "Invalid key identifier %s in declaration of Record type: %s\n", key.c_str(), type_name().c_str());

        // We do not allow declaration with unfinished type. The only exception is internal "TYPE"
        // key defined by AbstractRecord.
        if ( ! type.is_finished() )
            xprintf(PrgErr, "Unfinished type of declaring key: %s in Record type: %s\n", key.c_str(), type_name().c_str() );

        // make our own copy of type object allocated at heap (could be expensive, but we don't care)
        boost::shared_ptr<const KeyType> type_copy = boost::make_shared<KeyType>(type);

        data_->declare_key(key, type_copy, default_value, description);
    }

    /**
     * Same as previous method but without given default value (same as DefaultValue(DefaultValue::none) )
     */
    template <class KeyType>
    void declare_key(const string &key,
                            const KeyType &type,
                            const string &description)
    {
        declare_key(key,type, DefaultValue(), description);
    }

    /**
     * Finish declaration of the Record type. Now further declarations can be added.
     */
    void finish() { finished = true; }

    /**
     * @brief Implements @p Type:TypeBase::documentation.
     */
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;

    /**
     * Set made_extensive_doc = false for this Record and all its descendants.
     */
    virtual void  reset_doc_flags() const;

    virtual string type_name() const;

    virtual bool operator==(const TypeBase &other) const
    { return  typeid(*this) == typeid(other) &&
                     (type_name() == static_cast<const Record *>(&other)->type_name() );
    }


    /**
     * Interface to mapping key -> index in record. Returns index (in continuous array) for given key.
     */
    inline unsigned int key_index(const string& key) const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name().c_str());
        KeyHash key_h = key_hash(key);
        RecordData::key_to_index_const_iter it = data_->key_to_index.find(key_h);
        if (it != data_->key_to_index.end()) return it->second;
        else
            throw KeyNotFound() << KeyName_EI(key) << RecordName_EI(data_->type_name_);

        return size();
    }

    /**
     * Returns iterator to the key struct for given key string.
     *
     */
    inline KeyIter key_iterator(const string& key) const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name().c_str());
        return begin() + key_index(key);
    }

    inline KeyIter begin() const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name().c_str());
        return data_->keys.begin();
    }

    inline KeyIter end() const
    {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name().c_str());
        return data_->keys.end();
    }

    inline unsigned int size() const {
        ASSERT( finished, "Asking for information of unfinished Record type: %s\n", type_name().c_str());
        ASSERT( data_->keys.size() == data_->key_to_index.size(), "Sizes of Type:Record doesn't match. (map: %d vec: %d)\n", data_->key_to_index.size(), data_->keys.size());
        return data_->keys.size();
    }


protected:

    /**
     * Internal data class.
     */
    class RecordData {
    public:
        RecordData(const string & type_name_in, const string & description);

        void declare_key(const string &key,
                         boost::shared_ptr<const TypeBase> type,
                         const DefaultValue &default_value, const string &description);

        std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;

        void  reset_doc_flags() const;

        /// Database of valid keys
        std::map<KeyHash, unsigned int> key_to_index;
        typedef std::map<KeyHash, unsigned int>::const_iterator key_to_index_const_iter;

        /// Keys in order as they where declared.
        std::vector<struct Key> keys;


        /// Description of the whole record type.
        const string description_;
        const string type_name_;

        /**
         * This flag is set to true when documentation of the Record was called with extensive==true
         * and full description of the Record was produced.
         *
         * This member is marked 'mutable' since it doesn't change structure or description of the type. It only influence the output.
         */
        mutable bool made_extensive_doc;

    };


    /// Data handle.
    boost::shared_ptr<RecordData> data_;



};




class AbstractRecord : public Record {
protected:
    class ChildData {
    public:
        ChildData(const string &name)
        : selection_of_childs( name )
        {}

        Selection<unsigned int> selection_of_childs;
        vector< Record > list_of_childs;
    };

public:
    /**
     * Basic constructor. You has to provide @t type_name of the new declared Record type and
     * its @t description.
     */
    AbstractRecord(const string & type_name_in, const string & description)
    : Record(type_name_in, description),
      child_data_( boost::make_shared<ChildData>( type_name_in + "_selection" ) )
    {
        // make our own copy of type object allocated at heap (could be expensive, but we don't care)
        boost::shared_ptr< Selection<unsigned int> > type_copy = boost::make_shared< Selection<unsigned int> >(child_data_->selection_of_childs);

        data_->declare_key("TYPE", type_copy, DefaultValue(DefaultValue::obligatory),
                     "Sub-record selection.");

        finished=false;
    }

    void add_descendant(const Record &subrec)
    {
        child_data_->selection_of_childs.add_value(child_data_->list_of_childs.size(), subrec.type_name() );
        child_data_->list_of_childs.push_back(subrec);
    }

    const Record  &get_descendant(const string& name) const {
        unsigned int idx;

        idx = child_data_->selection_of_childs.name_to_value(name);
        ASSERT( idx < child_data_->list_of_childs.size() , "Size mismatch.\n");
        return child_data_->list_of_childs[idx];
    }

protected:
    boost::shared_ptr<ChildData> child_data_;

};

} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_RECORD_HH_ */
