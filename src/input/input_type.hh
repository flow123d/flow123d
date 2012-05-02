/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 *
 *  todo:
 *
 *  - testy novych metod ..
 *
 *  - Doxygen doc
 *
 *  - gettery pro AbstractRecord
 *  - testovani metod recordu a abstract recordu
 *  - match metody pro Scalarni typy
 *
 *  - zavest priznak pro record, ktery muze byt inicializovan z hodnoty jednoho klice (ten je treba zadat)
 *    a podobne pro pole.
 *    ...
 *
 *  - ?? predelat DefaultValue na hierarchii trid (opet problem s error hlaskou) tkato by se to hlasilo pri kompilaci.
 *    ( snad lepe factory staticke funkce defalut("..."), default(), default_obligatory() )
 *
 *  - Allow non integral values to selection ??
 *
 *  - have global list of Record and selection names and guarantee the they are unique, otherwise == can be incorrect.
 *
 *  When C++11 specification become more supported, we can introduce class Key  that should be constructed form
 *  constant string during compilation, i particular it should check validity of the key string and compute the hash.
 *  This can provide some speedup for reading if it will be needed (probably not).
 *
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include "type_base.hh"
#include "type_selection.hh"
#include "type_record.hh"

#endif /* INPUTTYPE_HH_ */
