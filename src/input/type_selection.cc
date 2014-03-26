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

Selection::Selection()
: data_(boost::make_shared<SelectionData>("EmptySelection"))
{
    close();
}



Selection::Selection(const Selection& other)
: Scalar(other), data_(other.data_)
{
    ASSERT( TypeBase::was_constructed(&other), "Trying to copy non-constructed Selection.\n");
}



Selection::Selection(const string &name, const string &desc)
: data_(boost::make_shared<SelectionData>(name))
{
    TypeBase::lazy_type_list().push_back( boost::make_shared<Selection>( *this) );
    data_->description_=desc;
}



Selection &Selection::add_value(const int value, const std::string &key, const std::string &description) {
    if (is_finished())
        xprintf(PrgErr, "Declaration of new name: %s in finished Selection type: %s\n", key.c_str(), type_name().c_str());

    data_->add_value(value, key, description);
    return *this;
}


const Selection & Selection::close() const {
    data_->finished=true;
    return *this;
}



bool Selection::valid_default(const string &str) const {
    if (! has_name(str))
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName( type_name() + " with values: "+key_list() ));
    return true;
}


bool Selection::is_finished() const {
    return data_->finished;
}


string Selection::type_name() const {
   return data_->type_name_;
}


string Selection::full_type_name() const {
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


string Selection::int_to_name(const int &val) const {
	finished_check();
	auto it = data_->value_to_index_.find(val);
	if (it != data_->value_to_index_.end())
		return data_->keys_[it->second].key_;
	else
		throw ExcSelectionValueNotFound() << EI_Value(val) << EI_Selection(*this);
}


Selection &Selection::copy_values(const Selection &sel)
{
	for (auto it = sel.begin(); it != sel.end(); ++it)
		add_value(it->value, it->key_, it->description_);

	return *this;
}



int Selection::from_default(const string &str) const {
    try {
        return name_to_int(str);
    } catch (ExcSelectionKeyNotFound &e) {
        THROW( ExcWrongDefault() << EI_DefaultStr( str ) << EI_TypeName( type_name() + " with values: "+key_list() ));
    }
    return -1;
}


string Selection::key_list() const {
    ostringstream os;
    for(unsigned int i=0; i<size(); i++) os << "'" <<data_->keys_[i].key_ << "' ";
    return os.str();
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





} // closing namespace Type
} // close namespace Input
