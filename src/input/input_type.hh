/*
 * InputType.hh
 *
 *  Created on: Mar 28, 2012
 *      Author: jb
 *
 *  todo:
 *
 *  - detailni osetreni existence data_
 *  - korektni destructory
 *  - testy novych metod ..
 *
 *  - Doxygen doc
 *
 *  - projit vsechny tridy a zkontrolovat, zda poskytuji dostatecny pristup ke svemu obsahu.
 *
 *  - zavest dedeni Recordu (neco jako copy constructor), ale s tim, ze rodicovsky Record by si pamatoval sve potomky
 *    a tim bychom se zbavili AbstractRecordu -
 *    v kopii recordu by se mohly prepisovat klice (jak to zaridit, nemel by se z chyby udelat warning?
 *
 *    1) Record si pamatuje sveho rodice. (neni nutne)
 *    2) Pri konstrukci si zkopiruje klice od rodice, potreba poznacit ty zkopirovane a umoznit jejich prepis.
 *    3) Pri konstrukci se potomek ohlasi vhodnou metodou predkovi.
 *    4) Base Record ma specialni klic TYPE typu Selection ten prirazuje jmenum potomky. Zaroven drzi pole ukazatelu na potomky.
 *       Tato Selection se muze menit i po uzavreni Recordu. Ovsem je to proti pravidlu, ze do Recordu s emohou pridavat
 *       jen uzavrane typy
 *
 *    behem cteni vstupu potrebuju prejit ke konkretnimu Recordu podle TYPE
 *    metody:
 *    Record(name, desc, inherit_from shared_ptr<AbstractRecord>) - inheriting constructor
 *    this copy
 *
 *    AbstractRecord - navic:
 *    protected + friend Record ? add_descendent( shared_ptr<Record> ) - ! can not add AbstractRecord, only one level inheritance
 *    shared_ptr<Record> get_descendent(string name)
 *
 *    nebo udelat potomka Record -> RecordBase OK.
 *
 *   - implementovat dedici konstruktor
 *   - testy AbstractRecordu
 *
 *  - zavest priznak pro record, ktery muze byt inicializovan z hodnoty jednoho klice (ten je treba zadat)
 *    a podobne pro pole.
 *    ...
 *
 *  - ?? predelat DefaultValue na hierarchii trid (opet problem s error hlaskou) tkato by se to hlasilo pri kompilaci.
 *
 *
 *  When C++11 specification become more supported, we can introduce class Key  that should be constructed form
 *  constant string during compilation, i particular it should check validity of the key string and compute the hash.
 *  This can provide some speedup for reading if it will be needed (probably not).
 */

#ifndef INPUTTYPE_HH_
#define INPUTTYPE_HH_

#include "type_base.hh"
#include "type_selection.hh"
#include "type_record.hh"

#endif /* INPUTTYPE_HH_ */
