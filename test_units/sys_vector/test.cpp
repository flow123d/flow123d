#include "system.hh"
#include <sys_vector.hh>

using namespace flow;

struct TestElementClass {
  int i;
  double a;
};

#define CHECK(x,...) INPUT_CHECK(x,__VA_ARGS__)

typedef Vector<TestElementClass> Cont;

int main(int argc, char * argv[]) {
  
  {// TEST CONTAINER Vector
  
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
  small_iter=iter;
  //test index of Vector
  CHECK( vec.index(small_iter) == 1, "");
  //test begin, end FullIter::index
  iter=vec.begin();
  // test FullIter prefix ++
  ++iter;
  CHECK(iter.index() == 1,"");
  iter=vec.end();
  // test FullIter prefix --
  --iter;
  CHECK(iter.index() == 1,"");
  CHECK(vec.size() == 2, "");
  // access items by [] and (), typecast form FullIter to Iter
  CHECK( vec[0].i == 1, "" );
  Vector<TestElementClass>::Iter small_it;
  small_it = vec(1);
  CHECK( small_it->a == 6.28, "" );
  // Do not know how to test Vector::reserve ??
 
  // test full_iter 
  CHECK( vec.full_iter(small_iter)->i == 2, "");
  
  // TEST remaining FullIter methods
  //
  Cont::FullIter full_iter_1(vec, small_iter);
  CHECK(full_iter_1->i == 2, "");
  // test dereference *
  CHECK( (*full_iter_1).i == 2, "");
  // test operator ==
  iter=vec.begin();
  CHECK( iter == vec.begin(), "");
  //test += and -=
  iter+=2;
  CHECK( iter == vec.end(), "");
  iter-=1;
  CHECK( iter->i == 2, "");
  }
  
  {// TEST VectorId
  typedef VectorId<TestElementClass> VecId;
  // test constructor
  VecId vec_id(2);
  // test inhereted types
  VecId::Iter small_iter;
  VecId::FullIter iter(vec_id);
  // test add_item
  iter=vec_id.add_item(30);
  iter->i = 1;
  iter->a = 3.14;
  iter = vec_id.add_item(40);
  iter->i = 2;
  iter->a = 6.28;
  // test find_id
  CHECK( vec_id.find_id(40)->i == 2,""); 
  // get_id, VectorId::FullIter::id
  CHECK(iter.id() == 40, "");
  small_iter=iter;
  CHECK( vec_id.get_id(small_iter) == 40, "");
  }
  
}  