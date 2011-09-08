/* Automatic differentiation subroutines - Forward mode.
   Jukka Toivanen, Jan Stebel
*/

#include <iostream>
#include <math.h>
#include "ad.h"
#include "sparse.h"

using namespace std;

int CVar::count = 0;

void CVar::init(int size=MAX_AD)
{
  ind = -1;
  asize = size;
  der=new double[asize]; 
  dermap=new int[asize];
  nder=0;
}

CVar::CVar()
{
  init();
  val = 0.0;
}

CVar::CVar(double ival)
{
  init();
  val = ival;
}

CVar::CVar(double ival, int size)
{
  init(size);
  val = ival;
}


CVar::CVar(const CVar& var)
{
  ind = -1;
  asize = var.asize;
  val = var.val;
  der = new double[asize];
  dermap = new int[asize];
  nder = var.nder;
  for (int i=0; i<var.nder; i++)
  {
    dermap[i]=var.dermap[i];
    der[i]=var.der[i];
  }
}


CVar::~CVar()
{
  delete[] dermap;
  delete[] der;
}



int CVar::getNder() const
{
  return nder;
}



void CVar::setDer(int ind, double val)
{
  int i;
  for (i=nder-1; i>=0; i--)
  {
    if (dermap[i]==ind)
    {
      der[i]=val;
      return;
    }
    if (dermap[i]<ind) break;
  }
  i++;

  for (int j=nder; j>i; j--)
  {
    dermap[j]=dermap[j-1];
    der[j]=der[j-1];
  }
  dermap[i]=ind;
  der[i]=val;

  if (nder==asize)
  {
    cout << "dermap full " << endl;
    exit(1);
  }

  nder++;
}

double CVar::getDer(int ind) const
{
  return der[ind];
}

void CVar::getDer(int low, int ncomp, double *d) const
{
  for (int i=0; i<ncomp; i++) d[i] = 0.0;

  for (int i=0; i<nder; i++)
    if (dermap[i]>=low && dermap[i]<low+ncomp) 
      d[dermap[i]-low]=der[i];
}  

double CVar::getDerOfVar(int ind) const {
  for (int i=0; i<nder; i++)
    if (dermap[i]==ind) return der[i];
  return 0;
}  

double CVar::getDerOfVar(const CVar& var) const {
  for (int i=0; i<nder; i++)
    if (dermap[i]==var.ind) return der[i];
  return 0;
}

int CVar::getDerVar(int ind) const {
  return dermap[ind];
}

int CVar::getInd() const {
  return ind;
}

inline cder CVar::getFirstCder(const CVar& param) {
  fi=0;
  gi=0;
  return getNextCder(param);
}

inline cder CVar::getNextCder(const CVar& param) {
  static cder res;

  if ((fi>=nder) && (gi>=param.nder)) { 
    res.ind=-1; return res; 
  }
  if (fi>=nder) {
    res.ind=param.dermap[gi];
    res.gd=param.der[gi];
    res.fd=0;
    gi++;
    return res;
  }
  if (gi>=param.nder) {
    res.ind=dermap[fi];
    res.fd=der[fi];
    res.gd=0;
    fi++;
    return res;
  }
  if (dermap[fi]<param.dermap[gi]) {
    res.ind=dermap[fi];
    res.fd=der[fi];
    res.gd=0;
    fi++;
    return res;
  } 
  if (dermap[fi]>param.dermap[gi]) {
    res.ind=param.dermap[gi];
    res.gd=param.der[gi];
    res.fd=0;
    gi++;
    return res;
  }

  res.ind=dermap[fi];
  res.fd=der[fi];
  res.gd=param.der[gi];
  fi++;
  gi++;
  return res;

}
    
void CVar::setIndependent()
{
  ind=count++;
  setDer(ind, 1.0);
}

void CVar::resetDer()
{
  nder = 0;
}


CVar& CVar::operator = (double nval)
{
  val = nval;
  nder = 0;
  return *this;
}  

CVar& CVar::operator = (const CVar& param)
{
  val = param.val;
  if (asize<param.asize)
  {
    delete[] dermap;
    delete[] der;
    asize=param.asize;
    der = new double[asize];
    dermap = new int[asize];
  }
  nder = param.nder;
  for (int i=0; i<nder; i++)
  {
    dermap[i]=param.dermap[i];
    der[i]=param.der[i];
  }
  return *this;
}

CVar CVar::operator+ (double d)
{
  CVar tmp(*this);
  tmp.val += d;
  return tmp;
}

CVar CVar::operator+ (const CVar& param)
{
  CVar tmp(val+param.val, max(asize, param.asize));
  cder der = getFirstCder(param);
  while (der.ind>=0)
  {
    tmp.setDer(der.ind, der.fd+der.gd);
    der=getNextCder(param);
  }
  return tmp;
}

CVar& CVar::operator+= (const CVar& param)
{
  val += param.val;
  assert(asize>=param.asize);
  cder der;
  der = getFirstCder(param);
  while (der.ind>=0) {
    setDer(der.ind, der.fd+der.gd);
    der=getNextCder(param);
  }
  return *this;
}

CVar CVar::operator- (double d)
{
  CVar tmp(*this);
  tmp.val -= d;
  return tmp;
}

CVar CVar::operator- (const CVar& param)
{
  CVar tmp(val-param.val, max(asize, param.asize));
  cder der;

  der = getFirstCder(param);

  while(der.ind>=0) {
    tmp.setDer(der.ind, der.fd-der.gd);
    der = getNextCder(param);
  }
  return tmp;
}

CVar CVar::operator - ()
{
  CVar tmp(0);
  return tmp-(*this);
}

CVar CVar::operator* (double d)
{
  CVar tmp(val*d, asize);
  for (int i=0;i<nder; i++) tmp.setDer(dermap[i], der[i]*d);
  return tmp;
}

CVar CVar::operator* (const CVar& param)
{
  CVar tmp(val*param.val, max(asize, param.asize));
  cder der;
  
  der = getFirstCder(param);
  while(der.ind>=0)
  {
    tmp.setDer(der.ind, der.fd*param.val+val*der.gd);
    der = getNextCder(param);
  }
  return tmp;
}

CVar CVar::operator / (const CVar& param)
{
  CVar tmp(val/param.val, max(asize, param.asize));
  cder der;

  der = getFirstCder(param);
  
  while(der.ind>=0)
  {
    tmp.setDer(der.ind, (der.fd*param.val-val*der.gd)/(param.val*param.val));
    der = getNextCder(param);
  }
  return tmp;
}

double sgn(double x)
{
  if (x > 0) { return 1.0; }
  else if (x < 0) { return -1.0; }
  else return 0.0;
}

CVar abs(CVar var)
{
  CVar tmp(fabs(var.val), var.asize);
  
  double s=sgn(var.val);
  for (int i=0; i<var.nder; i++) tmp.setDer(var.dermap[i], s*var.der[i]);

  return tmp;
}

CVar xabsx(CVar var)
{
  CVar tmp(fabs(var.val)*var.val, var.asize);
  
  double s=2.0*fabs(var.val);
  for (int i=0; i<var.nder; i++) tmp.setDer(var.dermap[i], s*var.der[i]);

  return tmp;
}

CVar sqrt(CVar param)
{
  CVar tmp(sqrt(param.val));
  for (int i=0; i<param.nder; i++) tmp.setDer(param.dermap[i], 0.5*param.der[i]/tmp.val);
  return tmp;
}

CVar sin(CVar var)
{
  double c;
  CVar tmp(sin(var.val));
  
  c = cos(var.val);
  for (int i=0; i<var.nder; i++) tmp.setDer(var.dermap[i], c*var.der[i]);
  return tmp;
}

CVar cos(CVar var)
{
  double c;
  CVar tmp(cos(var.val));
  c = -sin(var.val);
  for (int i=0; i<var.nder; i++) tmp.setDer(var.dermap[i], c*var.der[i]);
  return tmp;
}

CVar tanh(CVar var)
{
  double c;
  CVar tmp(tanh(var.val));
  
  c = 1.0 / (cosh(var.val) * cosh(var.val));
  for (int i=0; i<var.nder; i++) tmp.setDer(var.dermap[i], c*var.der[i]);

  return tmp;
}

CVar pow(CVar base, double exponent)
{
  CVar tmp(pow(base.val, exponent));
  double d = pow(base.val, exponent-1.0)*exponent;
  for (int i=0; i<base.nder; i++) tmp.setDer(base.dermap[i], d*base.der[i]);
  return tmp;
}

CVar pow(double base, CVar exponent)
{
  CVar tmp(pow(base, exponent.val));
  double d = tmp.val*log(base);
  for (int i=0; i<exponent.nder; i++) tmp.setDer(exponent.dermap[i], d*exponent.der[i]);
  return tmp;
}

CVec::CVec()
{
  asize=0;
  vec = new CVar[asize];
}

CVec::CVec(int size)
{
  asize=size;
  vec = new CVar[asize];
}

CVec::CVec(const CVec &old)
{
  asize = old.asize;
  vec = new CVar[asize];
  for (int i=0; i<asize; i++) vec[i] = old.vec[i];
}

CVec::~CVec()
{
  delete[] vec;
}

void CVec::getJac(int low, int ncomp, CSparseMat &jac) const
{
   int nnz=0, var, rows = asize, cols = ncomp, high = low+ncomp-1, tot;
   double d;

   for (int i=0; i<rows; i++)
     for (int j=0; j<vec[i].getNder(); j++)
     {
       var=vec[i].getDerVar(j);
       if (var>=low && var<=high) nnz++;
     }
     
   jac.resize(rows, cols, nnz);

   jac.ia[0]=0;
   tot = 0;
   for (int i=0; i<rows; i++)
   {
     jac.ia[i+1]=jac.ia[i];
     for (int j=0; j<vec[i].getNder(); j++)
     {
       var=vec[i].getDerVar(j);
       d = vec[i].getDer(j);
       if (var>=low && var<=high)
       {
         jac.v[tot]=d;
         jac.ja[tot]=var-low;
         jac.ia[i+1]++;
         tot++;
       }
     }
   }
}

/**
    Return the matrix
     d(this_i)
     ---------, i=M,...,M+m, j=low,...,low+ncomp
     dx_j
*/
void CVec::getSubJac(int M, int m, int low, int ncomp, CSparseMat &jac) const
{
   int nnz=0, var, rows = m, cols = ncomp, high = low+ncomp-1, tot;
   double d;

   for (int i=M; i<M+rows; i++)
     for (int j=0; j<vec[i].getNder(); j++)
     {
       var=vec[i].getDerVar(j);
       if (var>=low && var<=high) nnz++;
     }
     
   jac.resize(rows, cols, nnz);

   jac.ia[0]=0;
   tot = 0;
   for (int i=0; i<rows; i++)
   {
     jac.ia[i+1]=jac.ia[i];
     for (int j=0; j<vec[M+i].getNder(); j++)
     {
       var=vec[M+i].getDerVar(j);
       d = vec[M+i].getDer(j);
       if (var>=low && var<=high)
       {
         jac.v[tot]=d;
         jac.ja[tot]=var-low;
         jac.ia[i+1]++;
         tot++;
       }
     }
   }
}

void CVec::getVal(double *val) const
{
  for (int i=0; i<asize; i++)
    val[i]=vec[i].getVal();
}

double CVec::norm2()
{
  double nor = 0.0;
  for (int i=0; i<asize; i++) nor += vec[i].getVal()*vec[i].getVal();
  return sqrt(nor);
}

double CVec::min()
{
  double m = 1.0e20, v;
  for (int i=0; i<asize; i++)
  {
    v=vec[i].getVal();
    if (v<m) m=v;
  }
  return m;
}

double CVec::max()
{
  double m = -1.0e20, v;
  for (int i=0; i<asize; i++)
  {
    v=vec[i].getVal();
    if (v>m) m=v;
  }
  return m;
}

void CVec::setIndependent()
{
  for (int i=0; i<asize; i++) vec[i].setIndependent();
}

void CVec::resetDer()
{
  for (int i=0; i<asize; i++) vec[i].resetDer();
}


void CVec::resize(int size)
{
  if (asize!=size)
  {
    if (asize > 0) delete[] vec;
    asize=size;
    vec = new CVar[asize];
  }
}

int CVec::size()
{
  return asize;
}

CVec& CVec::operator = (const CVec& param)
{
  delete[] vec;
  asize = param.asize;
  vec = new CVar[asize];
  for (int i=0; i<asize; i++) vec[i]=param.vec[i];

  return *this;
}

CVec& CVec::operator = (const double& param)
{
  for (int i=0; i<asize; i++) vec[i]=param;

  return *this;
}


CVar& CVec::operator() (int loc)
{
  if (loc < 0 || loc >= asize)
  {
    cerr << "Array index out of range" << endl;
    throw 1;
    exit(1);
  }
  return vec[loc];
}

CVar CVec::operator | (const CVec& param)
{
  CVar tmp(0,asize);
  if (asize != param.asize)
  {
    cerr << "Vector size mismatch in dot product" << endl;
    exit(1); 
  }
  for (int i=0; i<asize; i++) tmp=tmp+vec[i]*param.vec[i];
  return tmp;
} 

