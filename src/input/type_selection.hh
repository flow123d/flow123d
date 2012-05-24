/*
 * type_selection.hh
 *
 *  Created on: May 1, 2012
 *      Author: jb
 */

#ifndef TYPE_SELECTION_HH_
#define TYPE_SELECTION_HH_

#include "system/exceptions.hh"

#include "system.hh"
#include "type_base.hh"

namespace Input {
namespace Type {

using namespace std;

/**
 * Base of Selection templates.
 */
class SelectionBase: public Scalar {
public:

    /**
     * Exceptions specific to this class.
     */
    TYPEDEF_ERR_INFO( EI_Selection, const SelectionBase * );
    DECLARE_EXCEPTION( ExcSelectionKeyNotFound, << "Key " << EI_KeyName::qval <<" not found in Selection:\n" <<  *EI_Selection::ref(_exc) );


    virtual int name_to_int(const string& name) const =0;
};

/**
 * @brief Template for classes storing finite set of named values.
 *
 * The primary purpose of this class is initialization of enum variables. In such a case
 * the template parameter @p Enum is the particular enum type. However, any type
 * that can be used as a key in std::map (comparable, copy and default constructable) can be used here.
 * In order to by type safe this class has to be templated by actual enum to which
 * the input value could be converted.
 *
 * This template assumes that Enum is an enum type or enum class type (in C++11).
 *
 * It would be nice to have specialization that could be constructed from enum types produced by CppEnumMacro.
 * Then we can drop add_value method.
 *
 * Future shows, if this is not too restrictive. Maybe, it is more practical drop the template and use just plain ints for values.
 *
 * TODO: We can not guarantee full compatibility of the Selection with corresponding Enum type
 *       the Selection can have fewer values since we can not get number of values in the Enum.
 *       Therefore we either have to move under C++11, where enum classes may provide elementary
 *       reflection or have Selection of simple ints.
 */
template<class Enum>
class Selection: public SelectionBase {
private:
    BOOST_STATIC_ASSERT( boost::is_integral<Enum>::value || boost::is_enum<Enum>::value );
public:



    /**
     * Default constructor. Empty handle.
     */
    Selection()
    {}

    /**
     * Creates a handle pointing to the new SelectionData.
     */
    Selection(const string &name) :
            data_(boost::make_shared<SelectionData>(name))
    {}

    /**
     * Adds one new @p value with name given by @p key to the Selection. The @p description of meaning of the value could be provided.
     */
    void add_value(const Enum value, const std::string &key, const std::string &description = "") {
        empty_check();
        if (is_finished())
            xprintf(PrgErr, "Declaration of new name: %s in finished Selection type: %s\n", key.c_str(), type_name().c_str());

        data_->add_value(value, key, description);
    }

    void finish() {
        empty_check();
        data_->finished = true;
    }

    virtual bool is_finished() const {
        empty_check();
        return data_->finished;
    }

    virtual std::ostream& documentation(std::ostream& stream, bool extensive = false, unsigned int pad = 0) const {
        if (! is_finished()) xprintf(Warn, "Printing documentation of unfinished Input::Type::Selection!\n");
        return data_->documentation(stream, extensive, pad);
    }

    virtual void reset_doc_flags() const {
        if (data_.use_count() != 0) data_->made_extensive_doc = false;
    }

    virtual string type_name() const {
        if (data_.use_count() == 0) return "empty_selection_handle";
        else return data_->type_name_;
    }

    virtual bool operator==(const TypeBase &other) const
    { return  typeid(*this) == typeid(other) &&
                     (type_name() == static_cast<const Selection *>(&other)->type_name() );
    }


    /***
     * Converts name (on input) to the value. Throws if the name do not exist.
     */
    virtual int name_to_int(const string &key) const {
        finished_check();
        KeyHash key_h = key_hash(key);
        typename SelectionData::key_to_index_const_iter it = data_->key_to_index_.find(key_h);
        if (it != data_->key_to_index_.end())
            return (data_->keys_[it->second].value);
        else
            throw ExcSelectionKeyNotFound() << EI_KeyName(key) << EI_Selection(this);
    }

    /**
     * Just check if there is a particular name in the Selection.
     */
    inline bool has_name(const string &key) const {
        finished_check();
        KeyHash key_h = key_hash(key);
        return (data_->key_to_index_.find(key_h) != data_->key_to_index_.end());
    }

    /***
     *  Check if there is a particular value in the Selection.
     */
    inline bool has_value(const Enum &val) const {
        finished_check();
        return (data_->value_to_index_.find(val) != data_->value_to_index_.end());
    }

    inline unsigned int size() const {
        finished_check();
        ASSERT( data_->keys_.size() == data_->key_to_index_.size(), "Sizes of Type:Selection doesn't match. (map: %d vec: %d)\n", data_->key_to_index_.size(), data_->keys_.size());
        return data_->keys_.size();
    }

private:

    inline void empty_check() const {
        ASSERT( data_.use_count() != 0, "Empty Selection handle.\n");
    }

    inline void finished_check() const {
        ASSERT( is_finished(), "Asking for information of unfinished Seleciton type: %s\n", type_name().c_str());
    }

    class SelectionData {
    public:
        SelectionData(const string &name)
        : type_name_(name), made_extensive_doc(false), finished(false)
        {}

        void add_value(const Enum value, const std::string &key, const std::string &description) {
            F_ENTRY;

            KeyHash key_h = TypeBase::key_hash(key);
            if (key_to_index_.find(key_h) != key_to_index_.end()) {
                xprintf(PrgErr, "Name '%s' already exists in Selection: %s\n", key.c_str(), type_name_.c_str());
                return;
            }
            value_to_index_const_iter it = value_to_index_.find(value);
            if (it != value_to_index_.end()) {
                xprintf(
                        PrgErr, "Value %d of new name '%s' conflicts with value %d of previous name '%s' in Selection: '%s'.\n", value, key.c_str(), keys_[it->second].value, keys_[it->second].key_.c_str(), type_name_.c_str());
                return;
            }

            unsigned int new_idx = key_to_index_.size();
            key_to_index_.insert(std::make_pair(key_h, new_idx));
            value_to_index_.insert(std::make_pair(value, new_idx));

            Key tmp_key = { new_idx, key, description, value };
            keys_.push_back(tmp_key);
        }

        std::ostream& documentation(std::ostream& stream, bool extensive, unsigned int pad) const
        {

            if (!extensive) {
                stream << "Selection '" << type_name_ << "' of " << keys_.size() << " values.";
            }
            if (extensive && !made_extensive_doc) {
                made_extensive_doc = true;

                stream << endl << "Selection '" << type_name_ << "' of " << keys_.size() << " values." << endl;
                stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << endl;
                // keys
                for (keys_const_iterator it = keys_.begin(); it != keys_.end(); ++it) {
                    stream << setw(4) << "" << it->key_ << " = " << it->value;
                    if (it->description_ != "")
                        stream << " (" << it->description_ << ")";
                    stream << endl;
                }
                stream << setw(pad) << "" << std::setfill('-') << setw(10) << "" << std::setfill(' ') << " " << type_name_ << endl;
            }
            return stream;
        }




        string type_name_;

        struct Key {
            unsigned int key_index;
            string key_;
            string description_;
            Enum value;
        };
        /// Map of valid keys to index.
        std::map<KeyHash, unsigned int> key_to_index_;
        typedef std::map<KeyHash, unsigned int>::const_iterator key_to_index_const_iter;

        /// Map of valid values to index.
        typename std::map<Enum, unsigned int> value_to_index_;
        typedef typename std::map<Enum, unsigned int>::const_iterator value_to_index_const_iter;

        std::vector<Key> keys_;
        typedef typename std::vector<struct Key>::const_iterator keys_const_iterator;

        /**
         * This flag is set to true when documentation of the Record was called with extensive==true
         * and full description of the Record was produced.
         *
         * This member is marked 'mutable' since it doesn't change structure or description of the type. It only influence the output.
         */
        mutable bool made_extensive_doc;

        bool finished;
    };

    boost::shared_ptr<SelectionData> data_;
};

} // closing namespace Type
} // closing namespace Input

#endif /* TYPE_SELECTION_HH_ */
