/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 * $Id: test.cpp 1055 2011-04-21 13:43:54Z jan.brezina $
 * $Revision: 1055 $
 * $LastChangedBy: jan.brezina $
 * $LastChangedDate: 2011-04-21 15:43:54 +0200 (Thu, 21 Apr 2011) $
 *
 */

#include <cppunit/config/SourcePrefix.h>
#include <cppunit/extensions/HelperMacros.h>

#include "system/system.hh"
#include <system/sys_vector.hh>

using namespace flow;

// Class for Vector and VectorId template parameter.
struct TestElementClass {
    int i;
    double a;
};

typedef Vector<TestElementClass> Cont;

#define CHECK(expr,msg) CPPUNIT_ASSERT_MESSAGE(msg,expr)

class SysVectorTest: public CPPUNIT_NS::TestFixture {

    CPPUNIT_TEST_SUITE( SysVectorTest);
    CPPUNIT_TEST( sys_vector_classes);
    CPPUNIT_TEST( sys_vector_id_classes);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp() {}
    void tearDown() {}

protected:
    void sys_vector_classes() {// TEST CONTAINER Vector

        // create a container
        // test constructors and types
        Cont vec_def;
        Cont vec(1);
        Cont::Iter small_iter;
        Cont::iterator true_iter;
        // test FullIter default constructor
        Cont::FullIter iter(vec);

        // test add_values
        // test dereference of FullIter
        iter = vec.add_item();
        iter->i = 1;
        iter->a = 3.14;
        iter = vec.add_item();
        iter->i = 2;
        iter->a = 6.28;

        // test assign and type cast of FullIter
        small_iter = iter;
        //test index of Vector
        CHECK(vec.index(small_iter) == 1, "");
        //test begin, end FullIter::index
        iter = vec.begin();
        // test FullIter prefix ++
        ++iter;
        CHECK(iter.index() == 1, "");
        iter = vec.end();
        // test FullIter prefix --
        --iter;
        CHECK(iter.index() == 1, "");
        CHECK(vec.size() == 2, "");
        // access items by [] and (), typecast form FullIter to Iter
        CHECK(vec[0].i == 1, "");
        Vector<TestElementClass>::Iter small_it;
        small_it = vec(1);
        CHECK(small_it->a == 6.28, "");
        // Do not know how to test Vector::reserve ??

        // test full_iter
        CHECK(vec.full_iter(small_iter)->i == 2, "");

        // TEST remaining FullIter methods
        //
        Cont::FullIter full_iter_1(vec, small_iter);
        CHECK(full_iter_1->i == 2, "");
        // test dereference *
        CHECK((*full_iter_1).i == 2, "");
        // test operator ==
        iter = vec.begin();
        CHECK(iter == vec.begin(), "");
        //test += and -=
        iter += 2;
        CHECK(iter == vec.end(), "");
        iter -= 1;
        CHECK(iter->i == 2, "");
    }

    void sys_vector_id_classes() {
        // TEST VectorId
        typedef VectorId<TestElementClass> VecId;
        // test constructor
        VecId vec_id(2);
        // test inhereted types
        VecId::Iter small_iter;
        VecId::FullIter iter(vec_id);
        // test add_item
        iter = vec_id.add_item(30);
        iter->i = 1;
        iter->a = 3.14;
        iter = vec_id.add_item(40);
        iter->i = 2;
        iter->a = 6.28;
        // test find_id
        CHECK(vec_id.find_id(40)->i == 2, "");
        // get_id, VectorId::FullIter::id
        CHECK(iter.id() == 40, "");
        small_iter = iter;
        CHECK(vec_id.get_id(small_iter) == 40, "");
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( SysVectorTest );
