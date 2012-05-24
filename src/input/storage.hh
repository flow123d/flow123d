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

#include "system.hh"
#include "system/exceptions.hh"

namespace Input {

TYPEDEF_ERR_INFO( EI_RequestedType, const string);
TYPEDEF_ERR_INFO( EI_StoredType, const string);
DECLARE_EXCEPTION(ExcStorageTypeMismatch, << "Storage type mismatch. You want value of type "
                                          << EI_RequestedType::qval << " but stored is value of type "
                                          << EI_StoredType::qval);


/**
 * This class specifies interface to a data storage (and currently implements it by json_spirit library).
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
    virtual void print(ostream &stream, int pad=0) const =0;

    virtual ~StorageBase();

};

/**
 * Simple array of heterogeneous values. The values are inserted as pointers
 * (no copies) that is possibly dangerous, but this is only internal class.
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
    virtual void print(ostream &stream, int pad=0) const;
    virtual ~StorageArray();
private:
    /// Forbids default constructor to have array set to NULL.
    StorageArray();
    vector<StorageBase *> array_;
};


class StorageBool : public StorageBase {
public:
    StorageBool(bool value);
    virtual bool get_bool() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(ostream &stream, int pad=0) const;
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
    virtual void print(ostream &stream, int pad=0) const;
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
    virtual void print(ostream &stream, int pad=0) const;
    virtual ~StorageDouble();

private:
    double value_;
};

class StorageString : public StorageBase {
public:
    StorageString(const string & value);
    virtual const string & get_string() const;
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(ostream &stream, int pad=0) const;
    virtual ~StorageString();

private:
    string value_;
};

class StorageNull : public StorageBase {
    virtual bool is_null() const;
    virtual StorageBase *deep_copy();
    virtual void print(ostream &stream, int pad=0) const;
    virtual ~StorageNull();

};

} // namespace Input

#endif /* STORAGE_HH_ */


