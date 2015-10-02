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
 * @file    storage_transpose.hh
 * @brief   
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
    		StorageBase const *source_storage, unsigned int vec_size);

    StorageBase *get_item(unsigned int index);

private:
    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::TypeBase *source_type,
    		StorageBase const *source_storage, unsigned int index);

    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::Record *source_type,
    		StorageBase const *source_storage, unsigned int index);
    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::AbstractRecord *source_type,
    		StorageBase const *source_storage, unsigned int index);
    StorageBase * modify_storage(const Type::TypeBase *target_type, const Type::Array *source_type,
    		StorageBase const *source_storage, unsigned int index);

    const Type::TypeBase *target_type_;
    const Type::TypeBase *source_type_;
    const StorageBase *source_storage_;
    unsigned int vec_size_;
};


} /* namespace Input */

#endif /* STORAGE_TRANSPOSE_HH_ */
