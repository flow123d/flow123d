/*
 * type_record.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_RECORD_HH_
#define TYPE_RECORD_HH_

#include "system/exceptions.hh"

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
     * Exceptions specific to this class.
     */
    TYPEDEF_ERR_INFO( EI_Record, Record );
    TYPEDEF_ERR_INFO( EI_RecordName, const string);
    DECLARE_EXCEPTION( ExcRecordKeyNotFound, << "Key " << EI_KeyName::qval <<" not found in Record:\n" <<  EI_Record::val );
    DECLARE_EXCEPTION( ExcDeriveNonEmpty, << "Can not derive from Record " << EI_RecordName::qval << " into"
            "non-empty Record:\n" << EI_Record::val );

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
    Record() {}

    /**
     * Basic constructor. You has to provide @t type_name of the new declared Record type and
     * its @t description.
     */
    Record(const string & type_name_in, const string & description);

    /**
     * Method to derive new Record from an AbstractRecord @p parent. This copy all keys from the @p parent and register the newly created Record
     * in the @p parent. You are free to overwrite copied keys, but you can not delete them.
     */
    void derive_from(AbstractRecord parent);

    /**
     * Declares a key of the Record with name given by parameter @p key, the type given by parameter @p type, default value by parameter @p default_value, and with given
     * @p description. The parameter @p type has to be any of descendants of TypeBase.
     *
     */
    template <class KeyType>
    void declare_key(const string &key,
                            const KeyType &type,
                            const Default &default_value, const string &description)
    {
        // ASSERT MESSAGE: The type of declared keys has to be a class derived from TypeBase.
        BOOST_STATIC_ASSERT( (boost::is_base_of<TypeBase, KeyType>::value) );

        empty_check();
        if (is_finished() ) xprintf(PrgErr, "Declaration of key: %s in finished Record type: %s\n", key.c_str(), type_name().c_str());

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
     * Same as previous method but without given default value (same as Default(Default::none) )
     */
    template <class KeyType>
    void declare_key(const string &key,
                            const KeyType &type,
                            const string &description)
    {
        declare_key(key,type, Default::optional(), description);
    }

    /**
     * Finish declaration of the Record type. Now further declarations can be added.
     */
    void finish() {
        empty_check();
        data_->finished = true;
    }

    virtual bool is_finished() const {
        empty_check();
        return data_->finished;
    }

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
        finished_check();
        KeyHash key_h = key_hash(key);
        RecordData::key_to_index_const_iter it = data_->key_to_index.find(key_h);
        if (it != data_->key_to_index.end()) return it->second;
        else
            THROW( ExcRecordKeyNotFound() << EI_KeyName(key) << EI_Record(*this) );

        return size();
    }


    /**
     * Returns iterator to the key struct for given key string.
     *
     */
    inline KeyIter key_iterator(const string& key) const
    {
        finished_check();
        return begin() + key_index(key);
    }

    /**
     * Returns iterator to the key struct for given key string.
     *
     */
    inline bool has_key_iterator(const string& key, KeyIter &it) const
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


    inline KeyIter begin() const
    {
        finished_check();
        return data_->keys.begin();
    }

    inline KeyIter end() const
    {
        finished_check();
        return data_->keys.end();
    }

    inline bool has_key(const string& key) const
    {
        return key_iterator(key) != end();
    }

    inline unsigned int size() const {
        finished_check();
        ASSERT( data_->keys.size() == data_->key_to_index.size(), "Sizes of Type:Record doesn't match. (map: %d vec: %d)\n", data_->key_to_index.size(), data_->keys.size());
        return data_->keys.size();
    }


protected:

    inline void empty_check() const {
        ASSERT( data_.use_count() != 0, "Empty Record handle.\n");
    }

    inline void finished_check() const {
        ASSERT( is_finished(), "Asking for information of unfinished Record type: %s\n", type_name().c_str());
    }

    /**
     * Internal data class.
     */
    class RecordData {
    public:
        RecordData(const string & type_name_in, const string & description);

        void declare_key(const string &key,
                         boost::shared_ptr<const TypeBase> type,
                         const Default &default_value, const string &description);

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

        bool finished;

    };


    /// Data handle.
    boost::shared_ptr<RecordData> data_;



};




class AbstractRecord : public Record {
protected:
    class ChildData {
    public:
        ChildData(const string &name)
        : selection_of_childs( boost::make_shared<Selection<unsigned int> > (name) )
        {}

        boost::shared_ptr< Selection<unsigned int> > selection_of_childs;
        vector< Record > list_of_childs;
    };

public:
    /**
     * Basic constructor. You has to provide @t type_name of the new declared Record type and
     * its @t description.
     */
    AbstractRecord(const string & type_name_in, const string & description);

    /**
     * This method close an AbstractRecord for any descendants (since they modify the parent). Maybe we should not use
     * a Selection for list of descendants, since current interface do not expose this Selection. Then this method
     * could be removed.
     */
    void no_more_descendants();

    void add_descendant(const Record &subrec);

    /**
     * @brief Implements @p Type:TypeBase::documentation.
     */
    virtual std::ostream& documentation(std::ostream& stream, bool extensive=false, unsigned int pad=0) const;

    /**
     * Set made_extensive_doc = false for this Record and all its descendants.
     */
    virtual void  reset_doc_flags() const;

    const Record  &get_descendant(const string& name) const;

    const Record  &get_descendant(unsigned int idx) const;


protected:
    boost::shared_ptr<ChildData> child_data_;

};


/*********************************************************
 * Implementation
 */



} // closing namespace Type
} // closing namespace Input




#endif /* TYPE_RECORD_HH_ */
