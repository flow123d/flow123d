#ifndef F123DTPEG_STRUCT_HH
#define F123DTPEG_STRUCT_HH

using namespace std;

class AbstractBoundary {
   virtual int     getId();
   virtual int     getType();
   virtual double  getValue();
   virtual double  getSigma();
   virtual short   getWhere();
   virtual int     getElmId();
   virtual short   getSidId();
   virtual short   getTag();
};

#endif
