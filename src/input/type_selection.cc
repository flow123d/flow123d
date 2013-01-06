/*
 * type_selection.cc
 *
 *  Created on: Aug 1, 2012
 *      Author: jb
 */

#include "input/type_selection.hh"

namespace Input {
namespace Type {

using std::string;

Selection &Selection::add_value(const int value, const std::string &key, const std::string &description) {
    empty_check();
    if (is_finished())
        xprintf(PrgErr, "Declaration of new name: %s in finished Selection type: %s\n", key.c_str(), type_name().c_str());

    data_->add_value(value, key, description);

    return *this;
}



void Selection::finish() {
    empty_check();
    data_->finished = true;
}



bool Selection::is_finished() const {
    empty_check();
    return data_->finished;
}



std::ostream& Selection::documentation(std::ostream& stream, DocType extensive, unsigned int pad) const {
    if (!is_finished())
        xprintf(Warn, "Printing documentation of unfinished Input::Type::Selection!\n");
    return data_->documentation(stream, extensive, pad);
}



void Selection::reset_doc_flags() const {
    if (data_.use_count() != 0)
        data_->made_extensive_doc = false;
}



string Selection::type_name() const {
    if (data_.use_count() == 0)
        return "empty_selection_handle";
    else
        return data_->type_name_;
}



bool Selection::operator==(const TypeBase &other) const {
    return typeid(*this) == typeid(other) && (type_name() == static_cast<const Selection *>(&other)->type_name());
}



int Selection::name_to_int(const string &key) const {
    finished_check();
    KeyHash key_h = key_hash(key);
    SelectionData::key_to_index_const_iter it = data_->key_to_index_.find(key_h);
    if (it != data_->key_to_index_.end())
        return (data_->keys_[it->second].value);
    else
        throw ExcSelectionKeyNotFound() << EI_KeyName(key) << EI_Selection(*this);
}



int Selection::from_default(const string &str) const {
    try {
        return name_to_int(str);
    } catch (ExcSelectionKeyNotFound &e) {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName(type_name()));
    }
}



void Selection::SelectionData::add_value(const int value, const std::string &key, const std::string &description) {
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



std::ostream& Selection::SelectionData::documentation(std::ostream& stream, DocType extensive, unsigned int pad) const
{

    if (extensive == record_key) {
        stream << "Selection '" << type_name_ << "' of " << keys_.size() << " values.";
    } else
    if (! made_extensive_doc) {
        made_extensive_doc = true;
        pad=0;

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



} // closing namespace Type
} // close namespace Input
