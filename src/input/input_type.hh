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
 * @file    input_type.hh
 * @brief   
 * @todo
 *  - explicit instantiation of templates in accessors - not so easy
 *
 *  - TYPE is obligatory key of descendants of an Abstract, for consistent documentation it should be reported as
 *    a obligatory key in these Records, however it is quite complicated to have it there , and it is not necessary for
 *    check of the input since we shall look for it explicitely.
 *
 *  - Documentation of Abstract should contain TYPE and common keys, descendants should report only nonderived keys.
 *
 *  - better Doxygen doc
 *
 *  - detailni popis pouziti deklaraci typu na systemu trid (v TypeBase)
 *
 *  - Implementovat selection tak, aby typem Enum byla templatovana jen vnitrni datova struktura ( a pochopitelne access metody, ktere potrebuji
 *    presny typ, tj
 *    ...
 *
 *  - have global list of Record and selection names and guarantee the they are unique, otherwise == can be incorrect.
 *
 *  - when creating a "unique instance" of a lazy type we should check that its name is unique (in derived records we should
 *    distinguish short_name used in Abstract TYPE selection, and full_name that includes name of the parent Abstract.
 *    This is important to prevent Record derive from different local instances of Abstract.
 *
 *  When C++11 specification become more supported, we can introduce class Key  that should be constructed form
 *  constant string during compilation, in particular it should check validity of the key string and compute the hash.
 *  This can provide some speedup for reading if it will be needed (probably not).
 *
 *  - under C++11 we can also construct Selection s from the enum type
 *
 *  - Question: allow implicit value of the TYPE key (and selection the descendant) depending on the keys that appears on the input.
 *  - Question: Support keys with multiple possible types:
 *    we declare more keys like: region_set_INTEGER, region_set_STRING, region_set_ARRAY, but user specify just
 *    region_set=..., and appropriate key is used.
 *
 *  - non-polymorph inheritance of Records, in fact just copy of keys of one record to another + posibility to extend them
 *
 *  Notes:
 *  - copy constructors and usage of Pimpl idiom for more complex types is necessary due to usage of shared pointers - to create them we need copies
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include "type_base.hh"
#include "type_selection.hh"
#include "type_record.hh"
#include "type_abstract.hh"
#include "type_generic.hh"

#endif /* INPUTTYPE_HH_ */
