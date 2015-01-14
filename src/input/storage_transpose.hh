/*
 * storage_modifier.hh
 *
 *  Created on: Feb 13, 2014
 *      Author: jb
 */

#ifndef STORAGE_TRANSPOSE_HH_
#define STORAGE_TRANSPOSE_HH_

#include "input/storage.hh"
#include "input/type_base.hh"
#include "input/type_record.hh"

namespace Input {

class StorageTranspose {
public:
	/**
	 * Constructor
	 */
    StorageTranspose(const Type::TypeBase *target_type, const Type::TypeBase *source_type,
    		StorageBase *source_storage, unsigned int vec_size);

    StorageBase *get_item(unsigned int index);

private:
    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::TypeBase *source_type,
    		StorageBase *source_storage, unsigned int index);

    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::Record *source_type,
    		StorageBase *source_storage, unsigned int index);
    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::AbstractRecord *source_type,
    		StorageBase *source_storage, unsigned int index);
    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::Array *source_type,
    		StorageBase *source_storage, unsigned int index);

    const Type::TypeBase *target_type_;
    const Type::TypeBase *source_type_;
    StorageBase *source_storage_;
    unsigned int vec_size_;
};


} /* namespace Input */

#endif /* STORAGE_TRANSPOSE_HH_ */
