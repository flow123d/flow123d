/*
 * storage_modifier.cc
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#include "input/storage_modifier.hh"


namespace Input {

const StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::TypeBase *source_type,
		const StorageBase *source_storage, unsigned int index, unsigned int vec_size) {

	// dispatch types
    if (typeid(*source_type) == typeid(Type::Record)) {
        return modify_storage(target_type, static_cast<const Type::Record *>(source_type),
        		source_storage, index, vec_size);
    } else
    if (typeid(*source_type) == typeid(Type::Array)) {
        return modify_storage(target_type, static_cast<const Type::Array *>(source_type),
        		source_storage, index, vec_size);
    } else
    if (typeid(*source_type) == typeid(Type::Integer)) {
        ASSERT( typeid(*source_storage) == typeid(StorageInt), "Incompatible type of storage. Must be Integer!\n");
    	return source_storage;
    } else
    if (typeid(*source_type) == typeid(Type::Double)) {
        ASSERT( typeid(*source_storage) == typeid(StorageDouble), "Incompatible type of storage. Must be Double!\n");
    	return source_storage;
    } else
    if (typeid(*source_type) == typeid(Type::Bool)) {
        ASSERT( typeid(*source_storage) == typeid(StorageBool), "Incompatible type of storage. Must be Boolean!\n");
    	return source_storage;
    } else
    if (typeid(*source_type) == typeid(Type::Selection)) {
        ASSERT( typeid(*source_storage) == typeid(StorageInt), "Incompatible type of storage. For selection must be Integer!\n");
    	return source_storage;
    } else {
    	const Type::AbstractRecord * abstract_record_type = dynamic_cast<const Type::AbstractRecord *>(source_type);
    	if (abstract_record_type != NULL ) {
            return modify_storage(target_type, abstract_record_type, source_storage, index, vec_size);
    	}

        const Type::String * string_type = dynamic_cast<const Type::String *>(source_type);
        if (string_type != NULL ) {
            ASSERT( typeid(*source_storage) == typeid(StorageString), "Incompatible type of storage. Must be String!\n");
        	return source_storage;
        }

        // default -> error
        xprintf(Err,"Unknown descendant of TypeBase class, name: %s\n", typeid(source_type).name());
    }

	return new StorageNull();
}

const StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::Record *source_type,
		const StorageBase *source_storage, unsigned int index, unsigned int vec_size) {

	ASSERT( typeid(*target_type) == typeid(Type::Record), "Incompatible type of target type. Must be Record!\n");

	Type::Record::KeyIter target_it= (static_cast<const Type::Record *>(target_type))->begin();
	StorageArray *storage_array = new StorageArray(source_type->size());

	for( Type::Record::KeyIter source_it= source_type->begin();
			source_it != source_type->end();
			++source_it, ++target_it) {

		const StorageBase * sb = modify_storage( target_it->type_.get(), source_it->type_.get(),
				source_storage->get_item(source_it->key_index), index, vec_size );
		storage_array->new_item(source_it->key_index, const_cast<StorageBase*>(sb));
	}

	return storage_array;
}


const StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::AbstractRecord *source_type,
		const StorageBase *source_storage, unsigned int index, unsigned int vec_size) {

	ASSERT( typeid(*target_type) == typeid(Type::AbstractRecord), "Incompatible type of target type. Must be AbstractRecord!\n");

	int descendant_index = source_storage->get_item(0)->get_int();
	const Type::AbstractRecord *target_arec = static_cast<const Type::AbstractRecord *>(target_type);

	return modify_storage( &( target_arec->get_descendant(descendant_index) ),
	                	   &( source_type->get_descendant(descendant_index) ),
	                	   source_storage, index, vec_size );
}


const StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::Array *source_type,
		const StorageBase *source_storage, unsigned int index, unsigned int vec_size) {

	if ( typeid(*target_type) == typeid(Type::Array) ) { // copy array
		const Type::TypeBase *source_array_type = &( source_type->get_sub_type() );
		const Type::TypeBase *target_array_type = &( static_cast<const Type::Array *>(target_type)->get_sub_type() );
		unsigned int array_size = source_storage->get_array_size();
		StorageArray *storage_array = new StorageArray(array_size);

		for (unsigned int i=0; i<array_size; i++) {
			const StorageBase * sb = modify_storage( target_array_type, source_array_type,
					source_storage->get_item(i), index, vec_size );
			storage_array->new_item(i, const_cast<StorageBase*>(sb));
		}

		return storage_array;
	}

	if ( *target_type == source_type->get_sub_type()) { // get member at index position
		ASSERT(index < vec_size, "Index of storage descendant is out of range.\n");
		ASSERT_EQUAL(source_storage->get_array_size(), vec_size);

		return modify_storage( target_type, &( source_type->get_sub_type() ),
				source_storage->get_item(index), index, vec_size );
	}

	ASSERT( false, "Incompatible type of target type. Must be Array or same type as subtype of source!\n");
	return new StorageNull();
}


} /* namespace Input */
