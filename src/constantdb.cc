/**
 * File:   constantdb.cc
 * Author: dalibor
 * Created on October 1, 2010, 11:07 PM
 *
 * ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 * ! !    DON'T USE THIS CLASS - it's just for my private using.     ! !
 * ! !    After finishing my work,                                   ! !
 * ! !     this class will be canceled without refund.               ! !
 * ! !                                                      Dalibor  ! !
 * ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 * 
 */

#include "constantdb.h"
#include "global_defs.h"
#include "system.hh"

ConstantDB* ConstantDB::instance = new ConstantDB();

/**
 *
 */
bool ConstantDB::isKeyChar(const char* key) {
    std::map<const char*, const char*>::iterator it;
    it = charMap.find(key);
    if (it == charMap.end()) {
        return false;
    } else {
        return true;
    }
}

void ConstantDB::setChar(const char* key, const char* value) {
    ASSERT(!isKeyChar(key), "Second insert of constant char '%s'.\n", key);
//    bool val = charMap[key];
    charMap[key] = value;
}

const char* ConstantDB::getChar(const char* key) {
    ASSERT(isKeyChar(key), "Constant char '%s' not found.\n", key);
    const char* value = charMap[key];
    return value;
}
//------------------------------------------------------------------------------

/**
 *
 */
bool ConstantDB::isKeyInt(const char* key) {
    std::map<const char*, int>::iterator it;
    it = intMap.find(key);
    if (it == intMap.end()) {
        return false;
    } else {
        return true;
    }
}

void ConstantDB::setInt(const char* key, int value) {
    ASSERT(!isKeyInt(key), "Second insert of constant int '%s'.\n", key);
    int val = intMap[key];
    intMap[key] = value;
}

int ConstantDB::getInt(const char* key) {
    ASSERT(isKeyInt(key), "Constant int '%s' not found.\n", key);
    int value = intMap[key];
    return value;
}
//------------------------------------------------------------------------------

/**
 *
 */
bool ConstantDB::isKeyDouble(const char* key) {
    std::map<const char*, double>::iterator it;
    it = doubleMap.find(key);
    if (it == doubleMap.end()) {
        return false;
    } else {
        return true;
    }
}

void ConstantDB::setDouble(const char* key, double value) {
    ASSERT(!isKeyDouble(key), "Second insert of constant double '%s'.\n", key);
    double val = doubleMap[key];
    doubleMap[key] = value;
}

double ConstantDB::getDouble(const char* key) {
    ASSERT(isKeyDouble(key), "Constant double '%s' not found.\n", key);
    double value = doubleMap[key];
    return value;
}
//------------------------------------------------------------------------------

/**
 *
 */
bool ConstantDB::isKeyObject(const char* key) {
    std::map<const char*, void*>::iterator it;
    it = objectMap.find(key);
    if (it == objectMap.end()) {
        return false;
    } else {
        return true;
    }
}

void ConstantDB::setObject(const char* key, void* value) {
    ASSERT(!isKeyObject(key), "Second insert of constant instance '%s'.\n", key);
    bool val = objectMap[key];
    objectMap[key] = value;
}

void* ConstantDB::getObject(const char* key) {
    ASSERT(isKeyObject(key), "Constant instance '%s' not found.\n", key);
    void* value = objectMap[key];
    return value;
}
