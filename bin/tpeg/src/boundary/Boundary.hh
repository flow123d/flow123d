#ifndef F123DTPEG_BOUNDARY_HH
#define F123DTPEG_BOUNDARY_HH

using namespace std;

#define CONDITION_ON_SIDE 2

#define DIRICHLET_BC      1
#define NEUMANN_BC        2
#define NEWTON_BC         3

class Boundary {
private:
   void init(int, int, double, int, int, int, int);

   int     id;
   int     type;
   double  value;
   double  sigma;
   int     where;
   int     elmId;
   int     sidId;
   int     tag;

public:
   Boundary(int, int, double, int, int, int, int);
   Boundary(int, int, double, double, int, int, int, int);
   ~Boundary();

   int     getId();
   int     getType();
   double  getValue();
   double  getSigma();
   int     getWhere();
   int     getElmId();
   int     getSidId();
   int     getTag();
};

#endif
