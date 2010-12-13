#include "system.hh"

int main(int argc, char * argv[]) {
  
  int *i = xnew int (2);
  INPUT_CHECK( *i == 2 , "new operator ... don't pass\n");
  delete i;
  
  int *ii = xnew int [10];
  INPUT_CHECK (ii[0] == 0, "new[] operator ... don't pass\n");
  delete ii;
  
  int (*iii)[3];
  iii=xnew (int[3]);
}  