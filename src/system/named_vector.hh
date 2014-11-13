/*
 * named_vector.hh
 *
 *  Created on: 16.10., 2014
 *      Author: pe
 */

#ifndef NAMED_VECTOR_H_
#define NAMED_VECTOR_H_

#include <system/exceptions.hh>
#include "input/accessors.hh"


template<class T>
class NamedVector
{
public:
    TYPEDEF_ERR_INFO(EI_Member, std::string);
    TYPEDEF_ERR_INFO(EI_Index, unsigned int);
    DECLARE_EXCEPTION(ExcUnknownMember, << "NamedVector has no member with the name: " << NamedVector::EI_Member::qval);
    DECLARE_EXCEPTION(ExcIndexOutOfBounds, << "Index exceeded NamedVector size: " << NamedVector::EI_Index::val);
    
    /// Iterator over the vector.
    typedef typename std::vector<T>::iterator iterator;
    
    /// @name Operators.
    //@{ 
    /// Add a new member.
    NamedVector<T> &operator += (T &new_member);

    /// Add new members passed in another vector @p another.
    NamedVector<T> &operator += (const NamedVector<T> &another);
    
    /// Returns reference to the member given by @p member_name.
    /**
     * Throws exception if the member with given name is not found.
     */
    T &operator[](const std::string &member_name) const;
    
    /// Returns reference to the member given by @p index in the vector.
    /**
     * It is faster then looking for the member by name.
     * Throws exception if the index exceeds the vector length.
     */
    T &operator[](unsigned int index) const;
    
    //@}
    
    /// Returns pointer to the member given by name @p member_name.
    /**
     * Returns nullptr if not found.
     */
    T *member(const std::string &member_name) const;
    
    /// Looks for the member with name @p member_name.
    /** 
     * Returns the index in the vector of the member with name @p member_name.
     * The found member is returned by a pointer @p member_out.
     */
    unsigned int idx(const std::string &member_name, T* member_out) const;
    
    /// Create a new NamedVector as a subset of *this. 
    /**
     * The new NamedVector contains members with names given by the vector @p names.
     * If @p strict is true (default), then it will throw exception when a member with one of the names 
     * is not found. If @p strict is set false, then the names, which were not found, are omitted.
     */
    NamedVector<T> subset(std::vector<std::string> names, bool strict=true) const;
    
    /// Returns the size of the vector.
    inline unsigned int size() const { return members_.size(); }
        
    ///Provides read access to the data.
    const std::vector<T> &data()
    {return members_;}
    
    inline iterator begin() {return members_.begin();}
    inline iterator end() {return members_.end();}
    
    inline bool empty() {return begin() == end();}
    
    inline void clear(){ members_.clear();}
    
protected:
    std::vector<T *> members_;
    //std::vector<std::string> names_;
};

namespace NamedVectorTools {

template<class T>
static void initialize(NamedVector<T> &vector, const Input::Array& in_array)
{
    vector.clear();
    for (auto it = in_array.begin<Input::Record>(); it != in_array.end(); ++it)
    {
        vector += *(new T(*it));
    }
}

} //namespace NamedVectorTools


////////////////////////////////////////////////////////// IMPLEMENTATION

#include "system/global_defs.h"

template<class T>
T* NamedVector<T>::member(const std::string& member_name) const
{
    for(auto member : members_)
    {
        if (member->name() == member_name) return member;
    }
    return nullptr;
}

template<class T>
T& NamedVector<T>::operator[](const std::string& member_name) const
{
    T *found_member = member(member_name);
    if (found_member) return *found_member;

    THROW(ExcUnknownMember() << NamedVector<T>::EI_Member(member_name));
    return *(members_[0]); // formal to prevent compiler warning
}

template<class T>
T& NamedVector<T>::operator[](unsigned int index) const
{
    //ASSERT(index < members_.size(),"Member index is out of bounds.");
    if(index < members_.size())
        return members_[index]; 
    
    THROW(ExcIndexOutOfBounds() << NamedVector<T>::EI_Index(index));
    return *(members_[0]); // formal to prevent compiler warning
}


template<class T>
NamedVector< T >& NamedVector<T>::operator+=(T& new_member)
{
    T *found_member = member(new_member.name());
    if (found_member) {
        ASSERT(&new_member == found_member, "Another member of the same name exists when adding member: %s\n",
                new_member.name().c_str());
    } else {
        members_.push_back(&new_member);
    }
    return *this;
}

template<class T>
NamedVector< T >& NamedVector<T>::operator+=(const NamedVector< T >& another_vector)
{
    for(auto member : another_vector.members_) this->operator +=(*member);
    return *this;
}

template<class T>
NamedVector< T > NamedVector<T>::subset(vector< string > names, bool strict) const
{
    NamedVector<T> subset;
    if(strict)  //can throw exception if not found
        for(auto name : names) subset += (*this)[name];
    else        //this will not throw exception if member is not found
    {
        T* found_member;
        for(auto name : names) 
        {
            found_member = member(name);
            if(found_member)
                subset += *found_member;
        }
    }
    return subset;
}

template<class T>
unsigned int NamedVector<T>::idx(const string& member_name, T* member_out) const
{
    unsigned int idx = 0;
    for(auto member : members_)
    {
        if (member->name() == member_name) 
        {   member_out = member;
            return idx;
        }
        idx++;
    }
    member_out = nullptr;
    return -1;
}


#endif /* NAMED_VECTOR_H_ */