/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    storage_transpose.cc
 * @brief   
 */

#include <cstdint>

#include "input/storage_transpose.hh"
#include "input/input_type.hh"


namespace Input {

StorageTranspose::StorageTranspose(const Type::TypeBase *target_type, const Type::TypeBase *source_type,
		const StorageBase *source_storage, unsigned int vec_size)
: target_type_(target_type), source_type_(source_type), source_storage_(source_storage), vec_size_(vec_size)
{}


StorageBase * StorageTranspose::get_item(unsigned int index) {
	return modify_storage(target_type_, source_type_, source_storage_, index);
}


StorageBase * StorageTranspose::modify_storage(const Type::TypeBase *target_type, const Type::TypeBase *source_type,
		const StorageBase *source_storage, unsigned int index) {

	if (typeid(*source_storage) == typeid(StorageNull)) {
		return new StorageNull();
	}

	// dispatch types
    if (typeid(*source_type) == typeid(Type::Record)) {
        return modify_storage(target_type, static_cast<const Type::Record *>(source_type),
        		source_storage, index);
    } else
    if (typeid(*source_type) == typeid(Type::Array)) {
        return modify_storage(target_type, static_cast<const Type::Array *>(source_type),
        		source_storage, index);
    } else
    if (typeid(*source_type) == typeid(Type::Integer)) {
    	OLD_ASSERT( typeid(*target_type) == typeid(Type::Integer), "Incompatible type of target. Must be Type::Integer!\n");
    	OLD_ASSERT( typeid(*source_storage) == typeid(StorageInt), "Incompatible type of storage. Must be Integer!\n");
    	return source_storage->deep_copy();
    } else
    if (typeid(*source_type) == typeid(Type::Double)) {
    	OLD_ASSERT( typeid(*target_type) == typeid(Type::Double), "Incompatible type of target. Must be Type::Double!\n");
    	OLD_ASSERT( typeid(*source_storage) == typeid(StorageDouble), "Incompatible type of storage. Must be Double!\n");
    	return source_storage->deep_copy();
    } else
    if (typeid(*source_type) == typeid(Type::Bool)) {
    	OLD_ASSERT( typeid(*target_type) == typeid(Type::Bool), "Incompatible type of target. Must be Type::Boolean!\n");
    	OLD_ASSERT( typeid(*source_storage) == typeid(StorageBool), "Incompatible type of storage. Must be Boolean!\n");
    	return source_storage->deep_copy();
    } else
    if (typeid(*source_type) == typeid(Type::Selection)) {
    	OLD_ASSERT( typeid(*target_type) == typeid(Type::Selection), "Incompatible type of target. Must be Type::Selection!\n");
    	OLD_ASSERT( typeid(*source_storage) == typeid(StorageInt), "Incompatible type of storage. For selection must be Integer!\n");
    	return source_storage->deep_copy();
    } else {
    	const Type::Abstract * abstract_record_type = dynamic_cast<const Type::Abstract *>(source_type);
    	if (abstract_record_type != NULL ) {
            return modify_storage(target_type, abstract_record_type, source_storage, index);
    	}

        const Type::String * source_string = dynamic_cast<const Type::String *>(source_type);
        if (source_string != NULL ) {
        	const Type::String * target_string = dynamic_cast<const Type::String *>(target_type);
        	OLD_ASSERT( target_string != NULL, "Incompatible type of target. Must be Type::String!\n");
        	OLD_ASSERT( typeid(*source_storage) == typeid(StorageString), "Incompatible type of storage. Must be String!\n");
        	return source_storage->deep_copy();
        }

        // default -> error
        xprintf(Err,"Unknown descendant of TypeBase class, name: %s\n", typeid(source_type).name());
    }

	return new StorageNull();
}

StorageBase * StorageTranspose::modify_storage(const Type::TypeBase *target_type, const Type::Record *source_type,
		const StorageBase *source_storage, unsigned int index) {

	OLD_ASSERT( typeid(*target_type) == typeid(Type::Record), "Incompatible type of target type. Must be Record!\n");
	OLD_ASSERT( typeid(*source_storage) == typeid(StorageArray), "Incompatible type of storage. Must be Array!\n");
	ASSERT_EQUAL(source_type->size(), (static_cast<const Type::Record *>(target_type))->size());

	Type::Record::KeyIter target_it= (static_cast<const Type::Record *>(target_type))->begin();
	StorageArray *storage_array = new StorageArray(source_type->size());

	for( Type::Record::KeyIter source_it= source_type->begin();
			source_it != source_type->end();
			++source_it, ++target_it) {

		target_it->default_.has_same_type(source_it->default_);
		OLD_ASSERT(target_it->default_.has_same_type(source_it->default_), "Incompatible default value of source type and target type!\n");

		StorageBase * sb = modify_storage( target_it->type_.get(), source_it->type_.get(),
				source_storage->get_item(source_it->key_index), index );
		storage_array->new_item(source_it->key_index, sb);
	}

	return storage_array;
}


StorageBase * StorageTranspose::modify_storage(const Type::TypeBase *target_type, const Type::Abstract *source_type,
		const StorageBase *source_storage, unsigned int index) {

	const Type::Abstract *target_arec = dynamic_cast<const Type::Abstract *>(target_type);
	OLD_ASSERT( target_arec != NULL, "Incompatible type of target type. Must be Abstract!\n");
	OLD_ASSERT( typeid(*source_storage) == typeid(StorageArray), "Incompatible type of storage. Must be Array!\n");

	string descendant_name = source_storage->get_item(0)->get_string();

	return modify_storage( &( target_arec->get_descendant(descendant_name) ),
	                	   &( source_type->get_descendant(descendant_name) ),
	                	   source_storage, index );
}


StorageBase * StorageTranspose::modify_storage(const Type::TypeBase *target_type, const Type::Array *source_type,
		const StorageBase *source_storage, unsigned int index) {
	OLD_ASSERT( typeid(*source_storage) == typeid(StorageArray), "Incompatible type of storage. Must be Array!\n");

	if ( typeid(*target_type) == typeid(Type::Array) ) { // copy array
		const Type::TypeBase *source_array_type = &( source_type->get_sub_type() );
		const Type::TypeBase *target_array_type = &( static_cast<const Type::Array *>(target_type)->get_sub_type() );
		unsigned int array_size = source_storage->get_array_size();
		StorageArray *storage_array = new StorageArray(array_size);

		for (unsigned int i=0; i<array_size; i++) {
			StorageBase * sb = modify_storage( target_array_type, source_array_type,
					source_storage->get_item(i), index );
			storage_array->new_item(i, sb);
		}

		return storage_array;
	}

	if ( *target_type == source_type->get_sub_type()) { // get member at index position
		OLD_ASSERT(index < vec_size_, "Index of storage descendant is out of range.\n");
		ASSERT_EQUAL(source_storage->get_array_size(), vec_size_);

		return source_storage->get_item(index)->deep_copy();
	}

	OLD_ASSERT( false, "Incompatible type of target type. Must be Array or same type as subtype of source!\n");
	return new StorageNull();
}


} /* namespace Input */
