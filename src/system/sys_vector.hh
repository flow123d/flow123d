/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Vector classes to support both Iterator, index and Id access and creating co-located vectors.
 *
 * @todo Implement co-located vector class that takes reference of the original Vector and
 *       is of the same size. FullIterator should have method for accessing elements in the original vector.
 */
#ifndef SYS_VECTOR_HH_
#define SYS_VECTOR_HH_

#include <vector>
#include <map>
#include "global_defs.h"
#include "system.hh"

// OPEN NAME SPACE "flow"
namespace flow {

/**
 * @brief Iterator that keeps also reference to its container. Safer and provides indexes.
 *
 *  Full iterator into a Vector<T> container. It makes
 *  an envelop around std:vector<T>::iterator, but stores also reference of
 *  the container it points into. Consequently FullIter can return index.
 *  The container has to be provided when FullIterator is created and this container can not be changed
 *  during the life of the variable. This is intentional since more safe.
 *
 *  This class is larger then normal pointer i.e. Iter type and therefore is questionable
 *  to use it in the large arrays. To this and you can use Iter and later convert it to FullIter
 *  providing the container It points into.
 *
 */

template <class Cont>
class FullIterator
{
public:
    /// Constructor to make uninitialized FullIterator for a container.
    FullIterator(Cont &cont_par)
            : cont(cont_par),
              iter(cont_par.storage.end()) {}

    /// Constructor to make FullIterator from container and its std::vector iterator.
    /// Mainly for internal use.
    FullIterator(const Cont &cont_par,const  typename  Cont::iterator it)
        : cont(cont_par),
          iter(it) {}

    /// Constructor to make FullIterator from container and Iter i.e. T*
    FullIterator(Cont &cont_par,const typename Cont::Iter it)
        : cont(cont_par)
    {
        if (it == NULL) iter = cont_par.storage.end();
        else iter = cont_par.storage.begin() + cont.index(it);
    }

    /// Check validity of FullIterator.
    inline bool is_valid()
        { return this->iter != this->cont.storage.end(); }

    /**
     *  Get index of FullIterator in its container.
     *  Return invalid index -1 for undefined iterator.
     */
    inline int index() const
    { if (this->iter != this->cont.storage.end()) return this->iter - cont.storage.begin();
      else  return (-1);
    }

    /// Get reference to the container of the FullIterator.
    inline Cont &container()
    { return cont; }

    /// Assign operator. In debugging version check that containers of both FullIterator match.
    inline FullIterator & operator =(const FullIterator &orig)
    {
      //  container.storage.begin() == orig.container.storage.begin()
      ASSERT( (&(this->cont) == &(orig.cont)),"Can not change container of FulIter.\n");
      this->iter=orig.iter;
      return (*this);
    }

    /// Type cast to Iter i.e. standard pointer.
    inline operator typename Cont::Iter () const
            { return &(*(this->iter)); }

    ///  * dereference operator
    inline typename Cont::ElementType & operator *() const
            { return *(this->iter); }

    /// -> dereference operator
    inline typename Cont::ElementType * operator ->() const
            { return &(*(this->iter)); }

    inline FullIterator &operator += (const int shift)
            { this->iter+=shift; return (*this); }

    inline FullIterator &operator -= (const int shift)
            { this->iter-=shift; return (*this); }

    /// Comparison of two FullIterator.
    inline bool operator == (const FullIterator &it) const
            { return ( &(this->cont) == &(it.cont)) && (this->iter == it.iter); }

    inline bool operator != (const FullIterator &it) const
            { return ( &(this->cont) != &(it.cont)) || (this->iter != it.iter) ; }

    /// Prefix. Advance operator.
    inline FullIterator &operator ++ ()
    {
        ASSERT( iter != cont.storage.end(), "Can not advance iterator at the end.\n");
        ++(this->iter); return *this;
    }

    /// Postfix. Should not be used since involves iterator copy.
    inline FullIterator operator ++ (int);
    //{
    //    //xprinf(Warn, "Postfix advance opeartor should")
    //    ASSERT( iter==storage.end(), "Can not advance iterator at the end.\n");
    //    FullIterator x(*this); this->iter++; return x;
    //}

    /// Prefix. Advance to previous operator.
    inline FullIterator &operator -- ()
    {
        ASSERT( iter != cont.storage.begin(), "Can not advance iterator to previous of begin().\n");
        this->iter--; return *this;
    }

    /// Postfix. Should not be used since involves iterator copy.
    inline FullIterator operator -- (int);
    //{
    //    //xprinf(Warn, "Postfix advance opeartor should")
    //    ASSERT( iter==storage.begin(), "Can not advance iterator to previous of begin().\n");
    //    FullIterator x(*this); this->iter--; return x;
    //}

    // + - opeartors is better to define outside if we ever allow them.
protected:
    /// We have here reference to the container of the iterator. So, an instance of this class can not change the container.
    const Cont &cont;
    // Conatainer iterator.
    typename Cont :: iterator iter;
};

/**
 * Small extension of the class FullIterator<Cont>. Provides id() method to get id number of the iterator.
 * In contrast to FullIterator<> this is templated only by the element type of the container which is always VectorId<T>.
 */
template <class T>
class VectorId;

template <class T>
class FullIteratorId : public FullIterator< VectorId<T> > {
public:
        /// Constructor. Just call base class.
        FullIteratorId(VectorId<T> &cont)
                : FullIterator<VectorId<T> >(cont) {}
        /// Constructor. Just call base class.
        FullIteratorId(const VectorId<T> &cont, typename std::vector<T>::iterator it)
            : FullIterator< VectorId<T> >(cont, it) {}
        /// Constructor. Just call base class.
        FullIteratorId(VectorId<T> &cont, const typename VectorId<T>::Iter it)
            : FullIterator< VectorId<T> >(cont, it) {}
        /// Returns id of the iterator. See VectorId documentation.
        inline int id()
        { return this->cont.id_storage[ this->index() ]; }
};



/**
 *  @brief Envelop over std::vector, use enhanced iterators.
 *
 *  Differences compared to std::vector.
 *
 *  -# Implement only subset of std::vector methods.
 *
 *  -# iterator Iter - which is nothing else then T *. This we use since it can be initializeted to NULL as universal "non valid value".
 *     for iterators one has to comapre to std::vector.end(), but the container is not always at hand.
 *     This should be used in complex mesh structures to save space.
 *
 *  -# iterator FullIter - this is standard iterator which contains also reference to the containre it points to. Then it has more
 *     save "non valid value" method. Can be constructer from Iter and a container. You can get the index of this iterator.
 *     This is meant to be use as local variable.
 *     Using this instead of standard iterators It provides function to get index of an iterator
 *
 *  -# Provides method to get index of a pointer i.e. Iter.
 *
 *  -# add_item can be used to add a new element and get reference to it so that it can be initialized
 *     Useful even for T be class since when creating a new element the default constructor is called
 *     which can not fill it with values.
 *
 *  <b> Developer note: </b>
 * 
 *  It appears very dangerous to combine reallocating std::vectors with iterators implemented as pointers.
 *  Indeed, when vector array is reallocated old pointers become invalid. One possibility is to strictly
 *  distinguish crating of the array and later creating references into it. Another is to implement iterators by indexes
 *  which means also to use only FullIterators.
 */

template <class T>
class Vector
{
public:
    /// We have to export template type for FullIteraror.
    typedef T ElementType;
    /// Default iterator for this class.
    typedef T * Iter;
    /// For compatibility with std::algorithms we provide also standard iterators.
    /// Theoretically it should work also with Iter. Test that if you need it.
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    /// Type of FullIterator that should be used to iterate through this container.
    typedef FullIterator<Vector<T> > FullIter;


    /// Construct with reserve the space. Includes default constructor.
    Vector(unsigned int size = 0) : storage(0)
        {this->reserve(size);}

    /**
     * Add a new element into the container. Return iterator to it.
     * Provides a way to add new element to the container and get an iterator to set its values.
     * This way we need not copy a new element.
     */
    inline FullIter add_item()
    {
        storage.resize( storage.size() + 1 );
        return FullIter( *this, storage.begin() + storage.size() - 1 );
    }

    /**
     * For given pointer returns the index of the element in the Vector. The first element has zero index.
     */
    inline unsigned int index(Iter pointer) const
        {
          ASSERT( pointer >= &(storage.front()) && pointer <= &(storage.back()),
                  "Wrong pointer %d to obtain its index (%d, %d).\n",pointer, &(storage.front()), &(storage.back()));
          return ( pointer - &(storage.front()) );
        }

     /// Returns FullFullIterer of the first element.
     inline FullIter begin()
         { return FullIter( *this, storage.begin() ); }

     /// Returns FullFullIterer of the fictions past the end element.
     inline FullIter end()
         { return FullIter( *this, storage.end() ); }

     /// Returns size of the container. This is independent of the allocated space.
     inline unsigned int size() const
         { return storage.size(); }

     /// Returns size of the container. This is independent of the allocated space.
     inline void resize(const unsigned int size)
         { storage.resize(size); }


     /// Gets reference to the element specified by index.
     inline T & operator[](unsigned int idx)
     {
         ASSERT( idx < this->size(), "Index %d outside of Vector of size %d\n",idx, this->size());
         return storage[idx];
     }

     /// Gets reference to the element specified by index.
     inline const T & operator[](unsigned int idx) const
     {
         ASSERT( idx < this->size(), "Index %d outside of Vector of size %d\n",idx, this->size());
         return storage[idx];
     }


     /// Gets iterator of the element specified by index.
     inline FullIter operator()(unsigned int idx)
     {
         ASSERT( idx < this->size(), "Index %d outside of Vector of size %d\n",idx, this->size());
         return FullIter( *this, begin()+idx );
     }


     /// Reallocates the container space.
     inline void reserve(unsigned int size)
      {
         ASSERT( size >= this->size(), "Vector can not be reallocated into space %d smaller then its size %d\n",size,this->size());
         storage.reserve(size);
      }

     /// Alternative way to make FullIterer form Iter. This way you need not write
     /// the full FullIter type to call the constructor. Then the result can be assigned to
     /// suitable local FullIter variable.
     inline FullIter full_iter(Iter it)
         { return FullIter(*this, it); }

     /// Since storage is direct member it is deallocated by its own destructor.
     virtual ~Vector() {}

protected:
     /// Underlying vector container.
     std::vector<T> storage;

     friend  class FullIterator< Vector<T> >;
};


/**
 * @brief Small extension of Vector<T> container with support to Id numbers.
 *
 * This stores Id numbers of added elements
 * and can latter access elements by id. An Id can be any integer number that only has to be unique
 * for every element in the container. In particular Id numbers are need not to be continuous neither keep the ordering of the elements.
 * Main application is to keep Id numbers from input, which are only necessary to amke correct references during input and possibly to
 * generate consistent output. Id numbers are useless for calculations.
 *
 * <b> Developer Note: </b>
 * 
 *  I have tried to make one common base class VectorBase which
 * should have storage member and implements basic operations with it.
 * But I can not make it to return FullIterator which has to be initialized
 * with reference to the container so not VectorBase but some of its descendants.
 * Only way seems to construct FullIterator by any pointer which could be type cast
 * to the type of the container which is template parameter of FullIter.
 * But such constructor should be used only by VectorXXX classes, but we can not make
 * VectorBase<T, It> friend for any T and any It.
 * However all this seems to be too complicated and ugly so I and up with lot of repetitive code
 * in implementation of VectorXXX classes.
 *
 */



template <class T>
class VectorId
{
public:
    /// We have to export template type for FullIteraror.
    typedef T ElementType;
    /// Default iterator for this class.
    typedef T * Iter;
    /// For compatibility with std::algorithms we provide also standard iterators.
    /// theoreticaly it should work also with Iter. Test that if you need it.
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    /// Type of FullIterator that should be used to iterate through this container.
    typedef FullIteratorId< T > FullIter;

    /// Construct with reserve the space. Includes default constructor.
    VectorId(unsigned int size = 0) : storage(0), id_storage(0)
        {this->reserve(size);}

    /**
     * Adds a new element into the container. Return iterator to it to set vales of the new element.
     * This way we need not copy a new element. You are forced to set the id when adding a new element.
     */
    FullIter add_item(int id)
    {
        ASSERT( id_map.find(id) == id_map.end(), "Can not add item with id number %d since it already exists.", id);
        id_storage.push_back(id);
        //DBGMSG("Push id: %d\n",id);
        id_map[id]=this->size();

        this->storage.resize( this->storage.size() + 1 );
        return FullIter( *this, this->storage.begin() + this->storage.size() - 1 );
    }


    /**
     * For given pointer returns the index of the element in the Vector. The first element has zero index.
     */
    inline unsigned int index(const T * pointer) const
        {
          ASSERT( pointer >= &(storage.front()) && pointer <= &(storage.back()),
                "Wrong pointer %p to obtain its index (%p, %p).\n",pointer, &(storage.front()), &(storage.back()));
          return ( pointer - &(storage.front()) );
        }

     /// Returns FullFullIterer of the first element.
     /// Do not try to set these as const methods. Such has to produce const iterators
     /// which requires const anlaogy of FullIter - really ugly
     inline FullIter begin()
         { return FullIter( *this, storage.begin() ); }

     /// Returns FullFullIterer of the fictions past the end element.
     inline FullIter end()
         { return FullIter( *this, storage.end() ); }

     /// Returns size of the container. This is independent of the allocated space.
     inline unsigned int size() const
         { return storage.size(); }

     /// Gets reference to the element specified by index.
     inline T & operator[](unsigned int idx)
     {
         ASSERT( idx < this->size(), "Index %d outside of Vector of size %d\n",idx, this->size());
         return storage[idx];
     }


     /// Gets reference to the element specified by index.
     inline const T & operator[](unsigned int idx) const
     {
         ASSERT( idx < this->size(), "Index %d outside of Vector of size %d\n",idx, this->size());
         return storage[idx];
     }


     /// Gets iterator of the element specified by index.
     inline FullIter operator()(unsigned int idx)
     {
         ASSERT( idx < this->size(), "Index %d outside of Vector of size %d\n",idx, this->size());
         return FullIter( *this, begin()+idx );
     }


     /// Alternative way to make FullFullIterer form FullIterer. This way you need not write
     /// the full FullFullIterer typename to call the constructor. Then the result can be assigned to
     /// suitable local FullFullIterer variable.
     inline FullIter full_iter(Iter it)
         { return FullIter(*this, it); }

     /**
      * Returns FullIter of the element with given id. This is implemented by std::map which use balanced trees and ensure access in O(nlog n) time.
      */
     inline FullIter find_id(const int id)
     {
         map<int,unsigned int>::iterator iter = id_map.find(id);
         if ( iter == id_map.end() ) {
             return FullIter(*this, this->storage.end());
         } else
             return FullIter(*this, this->storage.begin() + iter->second);
     }

     /**
      * Returns FullIter of the element with given id. This is implemented by std::map which use balanced trees and ensure access in O(nlog n) time.
      */
     inline const T * find_id(const int id) const
     {
         map<int,unsigned int>::const_iterator iter = id_map.find(id);
         if ( iter == id_map.end() ) {
             return &(*(this->storage.end()));
         } else
             return &(*(this->storage.begin() + iter->second));
     }

    /** Returns Id of the element given by pointer i.e. Iter. FullIter i.e. FullIteratorId<T>
     * provides its own method for the same.
     */
    inline int get_id(Iter it) const
    {
        return *(id_storage.begin() + this->index(it));
    }

    inline int get_id(int idx) const
    {
        return *(id_storage.begin() + idx);
    }

    /// Reallocates the container space.
    inline void reserve(unsigned int size)
     {
        ASSERT( size >= this->size(), "Vector can not be reallocated into space %d smaller then its size %d\n",size,this->size());
        storage.reserve(size);
        id_storage.reserve(size);
     }

    /// Since storage is direct member it is deallocated by its own destructor.
    virtual ~VectorId() {}


private:
    /// Underlying vector container.
    std::vector<T> storage;

    std::vector<int> id_storage;         ///< Space to save id numbers.
    std::map<int, unsigned int> id_map;  ///< Maps id numbers to indexes into storage vector.

    friend class FullIterator< VectorId<T> >;
    friend class FullIteratorId<T>;
};



// CLOSE NAME SPACE "flow"
}


#endif /* SYS_VECTOR_HH_ */
