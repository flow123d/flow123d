/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 *
 *  todo:

 *  - If we require the types with complex construction are instantiated as static members,
 *    + we store pointers to them in the decalration phase, we may remove pimpl idiom, and disallow copies
 *    (careful! surly there were reason for it, try to forbid copy constr. first)
 *  - this way we can recognize keys with constructed record type in the declare_key, since they
 *    are already registered in LazyTypes static list, in this case we can use lazy initialization only if we need it
 *    Workaround: use reference (to static variable) when returning Record type so that declare_key get correct pointer
 *
 *  - type_use_case_test, + demonstrace vystupu dokumentace
 *  - match metody pro Scalarni typy + converze z default hodnot
 *  - operator-> for iterators to accessors
 *
 *  - TYPE is obligatory key of descendants of an AbstractRecord, for consistent documentation it should be reported as
 *    a obligatory key in these Records, however it is quite complicated to have it there , and it is not necessary for
 *    check of the input since we shall look for it explicitely.
 *  - Documentation of AbstractRecord should contain TYPE and common keys, descendants should report only nonderived keys.
 *  - Descendant of an AbstractRecord should not be used directly !!
 *
 *  - Doxygen doc
 *
 *
 *  - detailni popis pouziti deklaraci typu na systemu trid (v TypeBase)
 *
 *  - Implementovat selection tak, aby typem Enum byla templatovana jen vnitrni datova struktura ( a pochopitelne access metody, ktere potrebuji
 *    presny typ, tj
 *  - zavest priznak pro record, ktery muze byt inicializovan z hodnoty jednoho klice (ten je treba zadat)
 *    a podobne pro pole.
 *    ...
 *
 *  - ?? predelat Default na hierarchii trid (opet problem s error hlaskou) tkato by se to hlasilo pri kompilaci.
 *    ( snad lepe factory staticke funkce defalut("..."), default(), default_obligatory() )
 *
 *  - Allow non integral values to selection ??
 *
 *  - declare_key should take TypeBase as parameter to be more flexible, but is not necessary
 *
 *  - have global list of Record and selection names and guarantee the they are unique, otherwise == can be incorrect.
 *
 *  When C++11 specification become more supported, we can introduce class Key  that should be constructed form
 *  constant string during compilation, in particular it should check validity of the key string and compute the hash.
 *  This can provide some speedup for reading if it will be needed (probably not).
 *
 *
 *   - nevytvaret deklaraci Recordu a dalsich typu runtime, ale pri kompilaci
 *    tedy jako skutecnou hierarchii trid, to by umoznilo statickou kontrolu
 *    kompatibility typu v konstruktoru Iterator<T>
 *    ??? is it possible?
 *
 *
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include "type_base.hh"
#include "type_selection.hh"
#include "type_record.hh"

#endif /* INPUTTYPE_HH_ */
