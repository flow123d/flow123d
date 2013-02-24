/*
 * accessors_impl.hh
 *
 *  Created on: Aug 1, 2012
 *      Author: jb
 */

#ifndef ACCESSORS_IMPL_HH_
#define ACCESSORS_IMPL_HH_


namespace Input {

using std::string;

/******************************************************************************************
 * Implementation of Input::Record
 */
template <class Ret>
inline const Ret Record::val(const string &key) const {
    try {
        Type::Record::KeyIter key_it = record_type_.key_iterator(key);

        ASSERT( key_it->default_.is_obligatory() || key_it->default_.has_value_at_declaration(),
                "The key '%s' is declared as optional or with default value at read time,"
                " you have to use Record::find instead.\n", key.c_str());

        Iterator<Ret> it = Iterator<Ret>( *(key_it->type_), address_, key_it->key_index);
        return *it;
    }
    // we catch all possible exceptions
    catch (Type::Record::ExcRecordKeyNotFound & e) {
        throw;
    }
    catch (ExcTypeMismatch & e) {
        e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
        throw;
    }
    catch (ExcStorageTypeMismatch &e) {
        throw;
    }
    catch (ExcAccessorForNullStorage &e) {
        throw;
    }
}



template <class Ret>
inline const Ret Record::val(const string &key, const Ret default_val ) const {
    try {
        Type::Record::KeyIter key_it = record_type_.key_iterator(key);

        ASSERT( key_it->default_.has_value_at_read_time(),
                "The key %s is not declared with default value at read time,"
                " you have to use Record::val or Record::find instead.\n", key.c_str());

        Iterator<Ret> it = Iterator<Ret>( *(key_it->type_), address_, key_it->key_index);
        if (it)
            return *it;
        else
            return default_val;
    }
    // we catch all possible exceptions
    catch (Type::Record::ExcRecordKeyNotFound & e) {
        throw;
    }
    catch (ExcTypeMismatch & e) {
        e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
        throw;
    }
    catch (ExcStorageTypeMismatch &e) {
        throw;
    }
    catch (ExcAccessorForNullStorage &e) {
        throw;
    }
}



template <class Ret>
inline Iterator<Ret> Record::find(const string &key) const {
    try {
        Type::Record::KeyIter key_it = record_type_.key_iterator(key);
        return Iterator<Ret>( *(key_it->type_), address_, key_it->key_index);
    }
    // we catch all possible exceptions
    catch (Type::Record::ExcRecordKeyNotFound & e) {
        throw;
    }
    catch (ExcTypeMismatch & e) {
        e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
        throw;
    }
}

template <class Ret>
inline bool Record::opt_val(const string &key, Ret &value) const {
    try {
        Type::Record::KeyIter key_it = record_type_.key_iterator(key);
        Iterator<Ret> it=Iterator<Ret>( *(key_it->type_), address_, key_it->key_index);
        if (it) {
            value = *it;
        } else {
            return false;
        }
    }
    // we catch all possible exceptions
    catch (Type::Record::ExcRecordKeyNotFound & e) {
        throw;
    }
    catch (ExcTypeMismatch & e) {
        e << EI_CPPRequiredType(typeid(Ret).name()) << EI_KeyName(key);
        throw;
    }

    return true;
}


/******************************************************************************************
 * Implementation of Input::Array
 */

template <class ValueType>
inline Iterator<ValueType> Array::begin() const {
    try {
        return Iterator<ValueType>(array_type_.get_sub_type(), address_, 0);
    }
    catch (ExcTypeMismatch & e) {
        e << EI_CPPRequiredType(typeid(ValueType).name()) << EI_KeyName("begin()");
        throw e;
    }
}



inline IteratorBase Array::end() const {
    return IteratorBase(address_, address_->storage_head()->get_array_size());
}



inline unsigned int Array::size() const {
    return address_->storage_head()->get_array_size();
}



template <class Container>
void Array::copy_to(Container &out) const {
    out.clear();
    Iterator<typename Container::value_type> it = begin<typename Container::value_type>();

    for(;it != end(); ++ it) {
        out.push_back(*it);
    }
}


/******************************************************************************************
 * Implementation of Input::IteratorBase
 */

inline bool IteratorBase::operator == (const IteratorBase &that) const
        { return ( address_->storage_head()  == that.address_->storage_head()  && index_ == that.index_); }



inline bool IteratorBase::operator != (const IteratorBase &that) const
        { return ! ( *this == that ); }



inline IteratorBase::operator bool() const {
    const StorageBase *s = address_->storage_head()->get_item(index_);
    return ( s && ! s->is_null() );
}



inline unsigned int IteratorBase::idx() const {
    return index_;
}


/******************************************************************************************
 * Implementation of Input::Iterator<Type>
 */


template<class T>
inline Iterator<T> & Iterator<T>::operator ++() {
    index_++;
    return *this;
}

template<class T>
inline typename Iterator<T>::OutputType Iterator<T>::operator *() const {

	const Address *a = address_->down(index_);

    ASSERT(a->storage_head(), "NULL pointer in storage!!! \n");

    return internal::TypeDispatch < DispatchType > ::value(a, type_);
}

template<class T>
inline typename Iterator<T>::OutputType * Iterator<T>::operator ->() const {
    BOOST_STATIC_ASSERT(
            (boost::is_same < Record, OutputType > ::value || boost::is_same < AbstractRecord, OutputType > ::value
                    || boost::is_same < Array, OutputType > ::value));

    // we have to make save temporary
    temporary_value_ = this->operator*();
    return &(temporary_value_);

}


template<class T>
typename Iterator<T>::InputType Iterator<T>::type_check_and_convert(const Input::Type::TypeBase &type) {
    if (typeid(type) == typeid(InputType)) {
        return static_cast<const InputType &>(type);
    } else {
        THROW(ExcTypeMismatch() << EI_InputType(type.type_name()) << EI_RequiredType(typeid(InputType).name()));
    }
}


} // namespace Input

#endif /* ACCESSORS_IMPL_HH_ */
