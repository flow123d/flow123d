/*
 * storage.hh
 *
 *  Created on: May 5, 2012
 *      Author: jb
 */

#ifndef STORAGE_HH_
#define STORAGE_HH_

/**
 * Iterator - like intermediate class for access to the stored values. The storage is a tree where
 * one node can hold:
 * - int, double, bool
 * - pointer to string
 * - pointer to array of storages
 * - special state: NULL (no data)
 *   INCLUDE (have to read another file to provide the value, this may be possible only through particular readers)
 *   ...
 *
 * json_spirit use boost::variant for storing JSON tree and arrays are stored as vectors of these variants
 * not pointers to variants. This is probably more effective, but do not allow effective modifications of the
 * tree and also construction not involving copies is not very intuitive. Therefore we use our own storage and
 * copy the json_spirit tree into it.
 */

#include <iostream>
#include <string>
#include <vector>
#include "system/exceptions.hh"


namespace Input {

TYPEDEF_ERR_INFO( EI_RequestedType, const std::string);
TYPEDEF_ERR_INFO( EI_StoredType, const std::string);
DECLARE_EXCEPTION(ExcStorageTypeMismatch, << "Storage type mismatch. You want value of type "
                                          << EI_RequestedType::qval << " but stored is value of type "
                                          << EI_StoredType::qval);


/**
 * @brief Base class for nodes of a data storage tree.
 *
 * This class as well as its descendants is meant for internal usage only as part of the implementation of the input interface.
 *
 * The leave nodes of the data storage tree can be of types \p StorageBool, \p StorageInt, \p StorageDouble, and StorageNull.
 * The branching nodes of the tree are of type StorageArray. The data storage tree serves to store data with structure described
 * by Input::Type classes. Therefore it provides no way to ask for the type of stored data and an exception \p ExcStorageTypeMismatch
 * is thrown if you use
 * getter that do not match actual type of the node. Moreover, the tree can be only created using bottom-up approach and than can
 * not be modified ( this can change if we want to use same storage for buffered reading of large data). However, we provide a method for
 * deep copy of any subtree.
 *
 * @ingroup input
 */
class StorageBase {
public:
    virtual int get_int() const;
    virtual double get_double() const;
    virtual bool get_bool() const;
    virtual const std::string &get_string() const;
    virtual const StorageBase * get_item(const unsigned int index) const;
    virtual bool is_null() const =0;
    virtual unsigned int get_array_size() const;

    virtual StorageBase *deep_copy()=0;
    virtual void print(std::ostream &stream, int pad=0) const =0;

    virtual ~StorageBase();

};

/**
 * Simple array of heterogeneous values. The values are inserted as pointers
 * (no copies) that is possibly dangerous, but don't care as the Storage is meant for internal usage only.
 */
class StorageArray : public StorageBase {
public:
    StorageArray(unsigned int size);
    StorageArray(const StorageArray &); // deep copy for test purpose
    void new_item(unsigned int index, StorageBase* item);
    virtual const StorageBase * get_item(const unsigned int index) const;
    virtual unsigned int get_array_size() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(std::ostream &stream, int pad=0) const;
    virtual ~StorageArray();
private:
    /// Forbids default constructor to have array set to NULL.
    StorageArray();
    std::vector<StorageBase *> array_;
};


class StorageBool : public StorageBase {
public:
    StorageBool(bool value);
    virtual bool get_bool() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(std::ostream &stream, int pad=0) const;
    virtual ~StorageBool();
private:
    bool value_;
};

class StorageInt : public StorageBase {
public:
    StorageInt(int value);
    virtual int get_int() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(std::ostream &stream, int pad=0) const;
    virtual ~StorageInt();

private:
    int value_;
};

class StorageDouble : public StorageBase {
public:
    StorageDouble(double value);
    virtual double get_double() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(std::ostream &stream, int pad=0) const;
    virtual ~StorageDouble();

private:
    double value_;
};

class StorageString : public StorageBase {
public:
    StorageString(const std::string & value);
    virtual const std::string & get_string() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(std::ostream &stream, int pad=0) const;
    virtual ~StorageString();

private:
    std::string value_;
};

class StorageNull : public StorageBase {
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(std::ostream &stream, int pad=0) const;
    virtual ~StorageNull();

};

} // namespace Input

#endif /* STORAGE_HH_ */


