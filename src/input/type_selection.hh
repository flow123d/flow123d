/*
 * type_selection.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_SELECTION_HH_
#define TYPE_SELECTION_HH_

#include "system/exceptions.hh"

#include "system/system.hh"
#include "type_base.hh"

namespace Input {
namespace Type {

using std::string;



/**
 * @brief Template for classes storing finite set of named values.
 *
 * The primary purpose of this class is initialization of enum variables. Since C++ provides no reflection,
 * in particular no access to enum identifiers as strings, you has to construct the Selection object consistent with an enum you want to initialize.
 *
 * Similarly to Type::Record and Type::AbstractRecord the Selection class is only proxy to the actual data.
 *
 * Usage:
 @code
    enum Colors { blue, white };

    Selection colors("Colors");
    colors.add_value(blue, "blue");
    colors.add_value(white,"white","White color"); // with optional item description
    colors.finish();
 @endcode
 *
 *
 * TODO: We can not guarantee full compatibility of the Selection with corresponding Enum type
 *       the Selection can have fewer values since we can not get number of values in the Enum.
 *       Therefore we either have to move under C++11, where enum classes may provide elementary
 *       reflection or have Selection of simple ints.
 *
 * @ingroup input_types
 */

class Selection : public Scalar {
	friend class OutputBase;

public:
    /*
     * Exceptions specific to this class.
     */
    TYPEDEF_ERR_INFO( EI_Selection, const Selection );
    TYPEDEF_ERR_INFO( EI_Value, const int);
    DECLARE_EXCEPTION( ExcSelectionKeyNotFound,
            << "Key " << EI_KeyName::qval <<" not found in Selection:\n" <<  EI_Selection::val );
    DECLARE_EXCEPTION( ExcSelectionValueNotFound,
                << "Value " << EI_Value::val <<" not found in Selection:\n" <<  EI_Selection::val );

    /**
     * Structure for description of one key in selection
     */
    struct Key {
        unsigned int key_index;
        string key_;
        string description_;
        int value;
    };

    /**
     * Public typedef of constant iterator into array of keys
     */
    typedef std::vector<struct Key>::const_iterator keys_const_iterator;

    /**
     * Default constructor. Empty selection.
     */
    Selection();


    /**
     * Copy constructor.
     */
    Selection(const Selection& other);

    /**
     * Creates a handle pointing to the new SelectionData.
     */
    Selection(const string &name, const std::string &description = "");

    /**
     * Adds one new @p value with name given by @p key to the Selection. The @p description of meaning of the value could be provided.
     */
    Selection &add_value(const int value, const std::string &key, const std::string &description = "");


    /**
     * Close the Selection, no more values can be added.
     */
    const Selection &close() const;


    TypeHash content_hash() const   override;


    /// Implements \p TypeBase::is_finished
    virtual bool is_finished() const;

    /// Implements \p TypeBase::type_name
    virtual string type_name() const;

    /// Implements \p TypeBase::full_type_name
    virtual string full_type_name() const;

    /// Implements \p TypeBase::operator==  compare also Selection names.
    virtual bool operator==(const TypeBase &other) const;

    /**
     * Container-like access to the keys of the Record. Returns iterator to the first key.
     */
    inline keys_const_iterator begin() const;

    /**
     * Container-like access to the keys of the Record. Returns iterator to the last key.
     */
    inline keys_const_iterator end() const;

    /**
     * Returns iterator to the key struct for given key string.
     */
    inline keys_const_iterator key_iterator(const string& key) const;

    /**
     * Converts given value name \p key to the value. Throws exception if the value name does not exist.
     */
    int name_to_int(const string &key) const;

    /**
     * Returns value name for the given \p value. Throws exception if the value does not exist.
     */
    string int_to_name(const int &value) const;

    Selection &copy_values(const Selection &sel);

    /**
     * Same as \p Selection::name_to_int, but throws different exception, when string comes from default value.
     */
    int from_default(const string &str) const;

    /// Implements  @p Type::TypeBase::valid_defaults.
    virtual bool valid_default(const string &str) const;

    /**
     * Just check if there is a particular name in the Selection.
     */
    inline bool has_name(const string &key) const;

    /**
     *  Check if there is a particular value in the Selection.
     */
    inline bool has_value(const int &val) const;

    /**
     * Returns number of values in the Selection.
     */
    inline unsigned int size() const;


    bool finish()
        { ASSERT(data_->closed_, "Finished Selection '%s' must be closed!", this->type_name().c_str()); return true; }


    /// Implements \p TypeBase::is_closed
    virtual bool is_closed() const override;


    // Implements @p TypeBase::make_instance.
    virtual boost::shared_ptr<TypeBase> make_instance(std::vector<ParameterPair> vec = std::vector<ParameterPair>()) const;
private:

    /**
     * Assertion for finished Selection (methods are called in correct order).
     */
    void finished_check() const;

    /**
     * Used in error messaged, where we can not use desc(), which can lead to infinite loop due to TYPE selection of AbstractRecord.
     */
    string key_list() const;

    /**
     * Actual Selection data.
     */
    class SelectionData  {
    public:

        SelectionData(const string &name)
        : type_name_(name), closed_(false)
        {}

        void add_value(const int value, const std::string &key, const std::string &description);


        /// Name of the Selection.
        string type_name_;

        /// Map : valid value name to index.
        std::map<KeyHash, unsigned int> key_to_index_;
        typedef std::map<KeyHash, unsigned int>::const_iterator key_to_index_const_iter;

        /// Map : valid value to index.
        std::map<int, unsigned int> value_to_index_;
        typedef std::map<int, unsigned int>::const_iterator value_to_index_const_iter;

        /// Vector of values of the Selection
        std::vector<Key> keys_;

        /// Text description of the whole Selection object.
        std::string description_;

        /// Indicator of closed Selection.
        mutable bool closed_;
    };

    /// Handle to actual Selection data.
    boost::shared_ptr<SelectionData> data_;

};





/******************************************************************************
 * Implementation of inline methods.
 */

inline bool Selection::has_name(const string &key) const {
    finished_check();
    KeyHash key_h = key_hash(key);
    return (data_->key_to_index_.find(key_h) != data_->key_to_index_.end());
}



inline bool Selection::has_value(const int &val) const {
    finished_check();
    return (data_->value_to_index_.find(val) != data_->value_to_index_.end());
}



inline unsigned int Selection::size() const {
    finished_check();
    ASSERT( data_->keys_.size() == data_->key_to_index_.size(),
            "Sizes of Type:Selection doesn't match. (map: %ld vec: %ld)\n", data_->key_to_index_.size(), data_->keys_.size());
    ASSERT_EQUAL( data_->keys_.size(), data_->key_to_index_.size());
    return data_->keys_.size();
}




inline void Selection::finished_check() const {
    ASSERT(data_->closed_, "Accessing unfinished Selection '%s'\n", type_name().c_str() );
}



inline Selection::keys_const_iterator Selection::begin() const
{
    finished_check();
    return data_->keys_.begin();
}



inline Selection::keys_const_iterator Selection::end() const
{
    finished_check();
    return data_->keys_.end();
}


inline Selection::keys_const_iterator Selection::key_iterator(const string& key) const
{
    finished_check();
    return begin() + name_to_int(key);
}


} // closing namespace Type
} // closing namespace Input

#endif /* TYPE_SELECTION_HH_ */
