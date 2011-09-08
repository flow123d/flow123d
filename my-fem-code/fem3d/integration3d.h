/* Integration points (quadratures)
   Jan Stebel
*/

#ifndef _INTEGRATION3D_H__
#define _INTEGRATION3D_H__

#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include "mesh3d.h"

/** evaluation of function only, or gradient w.r.t. q or gradient w.r.t. a */
typedef enum
{
  EVAL_F, EVAL_FG, EVAL_FA
} eval_t;

// Definition of quadrature formulas

/**
  Structure for storing the quadrature formulas
*/
struct CQuadrature
{
  int n;
  const double *x, *y, *z, *w;
};

// 3D Gaussian quadrature formulas for tetrahedra (exact for 4-th order polynomials)
const double Q3d11px[] = { 0.2500000, 0.7857143, 0.07142857, 0.07142857, 0.07142857, 0.1005964, 0.1005964, 0.1005964, 0.3994034, 0.3994034, 0.3994034, 0, 1, 0, 0, 0.5, 0.5, 0,   0,   0,   0.5 };
const double Q3d11py[] = { 0.2500000, 0.07142857, 0.7857143, 0.07142857, 0.07142857, 0.1005964, 0.3994034, 0.3994034, 0.1005964, 0.1005964, 0.3994034, 0, 0, 1, 0, 0,   0.5, 0.5, 0,   0.5, 0   };
const double Q3d11pz[] = { 0.2500000, 0.07142857, 0.07142857, 0.7857143, 0.07142857, 0.3994034, 0.1005964, 0.3994034, 0.1005964, 0.3994034, 0.1005964, 0, 0, 0, 1, 0,   0,   0,   0.5, 0.5, 0.5 };
const double Q3d11pw[] = { -0.01315556, 0.007622222, 0.007622222, 0.007622222, 0.007622222, 0.02488889, 0.02488889, 0.02488889, 0.02488889, 0.02488889, 0.02488889, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
const CQuadrature Q3d11p = { 11, Q3d11px, Q3d11py, Q3d11pz, Q3d11pw };

const double Q3d15px[] = { 0.2500000, 0.0000000, 0.3333333, 0.3333333, 0.3333333, 0.7272727, 0.09090909, 0.09090909, 0.09090909, 0.4334498, 0.4334498, 0.4334498, 0.06655015, 0.06655015, 0.06655015, 0, 1, 0, 0, 0.5, 0.5, 0,   0,   0,   0.5  };
const double Q3d15py[] = { 0.2500000, 0.3333333, 0.0000000, 0.3333333, 0.3333333, 0.09090909, 0.7272727, 0.09090909, 0.09090909, 0.4334498, 0.06655015, 0.06655015, 0.4334498, 0.4334498, 0.06655015, 0, 0, 1, 0, 0,   0.5, 0.5, 0,   0.5, 0    };
const double Q3d15pz[] = { 0.2500000, 0.3333333, 0.3333333, 0.0000000, 0.3333333, 0.09090909, 0.09090909, 0.7272727, 0.09090909, 0.06655015, 0.4334498, 0.06655015, 0.4334498, 0.06655015, 0.4334498, 0, 0, 0, 1, 0,   0,   0,   0.5, 0.5, 0.5  };
const double Q3d15pw[] = { 0.03028368, 0.006026786, 0.006026786, 0.006026786, 0.006026786, 0.01164525, 0.01164525, 0.01164525, 0.01164525, 0.01094914, 0.01094914, 0.01094914, 0.01094914, 0.01094914, 0.01094914, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
const CQuadrature Q3d15p = { 15, Q3d15px, Q3d15py, Q3d15pz, Q3d15pw };

const double Q3d45px[] = { 0.2500000,  0.6175872,  0.1274709,  0.1274709,  0.1274709,
                           0.9037635,  0.03207883, 0.03207883, 0.03207883, 0.4502229,
			   0.4502229,  0.4502229,  0.04977710, 0.04977710, 0.04977710,
			   0.3162696,  0.3162696,  0.3162696,  0.1837304,  0.1837304,
			   0.1837304,  0.5132800,  0.5132800,  0.5132800,  0.02291779,
			   0.02291779, 0.02291779, 0.2319011,  0.2319011,  0.2319011,
			   0.2319011,  0.2319011,  0.2319011,  0.1937465,  0.1937465,
			   0.1937465,  0.7303134,  0.7303134,  0.7303134,  0.03797005,
			   0.03797005, 0.03797005, 0.03797005, 0.03797005, 0.03797005,
			   0, 1, 0, 0, 0.5, 0.5, 0,   0,   0,   0.5 };
const double Q3d45py[] = { 0.2500000,  0.1274709,  0.6175872,  0.1274709,  0.1274709,
                           0.03207883, 0.9037635,  0.03207883, 0.03207883, 0.4502229,
			   0.04977710, 0.04977710, 0.4502229,  0.4502229,  0.04977710,
			   0.3162696,  0.1837304,  0.1837304,  0.3162696,  0.3162696,
			   0.1837304,  0.02291779, 0.2319011,  0.2319011,  0.5132800,
			   0.2319011,  0.2319011,  0.5132800,  0.5132800,  0.02291779,
			   0.02291779, 0.2319011,  0.2319011,  0.7303134,  0.03797005,
			   0.03797005, 0.1937465,  0.03797005, 0.03797005, 0.1937465,
			   0.1937465,  0.7303134,  0.7303134,  0.03797005, 0.03797005,
			   0, 0, 1, 0, 0,   0.5, 0.5, 0,   0.5, 0 };
const double Q3d45pz[] = { 0.2500000,  0.1274709,  0.1274709,  0.6175872,  0.1274709,
                           0.03207883, 0.03207883, 0.9037635,  0.03207883, 0.04977710,
			   0.4502229,  0.04977710, 0.4502229,  0.04977710, 0.4502229,
			   0.1837304,  0.3162696,  0.1837304,  0.3162696,  0.1837304,
			   0.3162696,  0.2319011,  0.02291779, 0.2319011,  0.2319011,
			   0.5132800,  0.2319011,  0.02291779, 0.2319011,  0.5132800, 
			   0.2319011,  0.5132800,  0.02291779, 0.03797005, 0.7303134,
			   0.03797005, 0.03797005, 0.1937465,  0.03797005, 0.7303134,
			   0.03797005, 0.1937465,  0.03797005, 0.1937465,  0.7303134,
			   0, 0, 0, 1, 0,   0,   0,   0.5, 0.5, 0.5 };
const double Q3d45pw[] = { -0.03932701,   0.004081316,  0.004081316,  0.004081316,  0.004081316,
                            0.0006580868, 0.0006580868, 0.0006580868, 0.0006580868, 0.004384259,
                            0.004384259,  0.004384259,  0.004384259,  0.004384259,  0.004384259,
                            0.01383006,   0.01383006,   0.01383006,   0.01383006,   0.01383006,
                            0.01383006,   0.004240437,  0.004240437,  0.004240437,  0.004240437,
			    0.004240437,  0.004240437,  0.004240437,  0.004240437,  0.004240437,
			    0.004240437,  0.004240437,  0.004240437,  0.002238740,  0.002238740,
			    0.002238740,  0.002238740,  0.002238740,  0.002238740,  0.002238740,
			    0.002238740,  0.002238740,  0.002238740,  0.002238740,  0.002238740,
			    1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
const CQuadrature Q3d45p = { 45, Q3d45px, Q3d45py, Q3d45pz, Q3d45pw };

// 2D Gaussian quadrature formulas for triangles
const double Q2d4px[] = { 1./3,    0.2,   0.2,   0.6,   0, 1, 0, 0.5, 0.5, 0 };
const double Q2d4py[] = { 1./3,    0.6,   0.2,   0.2,   0, 0, 1, 0,   0.5, 0.5 };
const double Q2d4pw[] = { -27./96, 25./96, 25./96, 25./96, 1, 1, 1, 1,   1,   1 };
const CQuadrature Q2d4p = { 4, Q2d4px, Q2d4py, 0, Q2d4pw };

const double Q2d6px[] = { 0.8168476,  0.09157621, 0.09157621, 0.1081030, 0.4459485, 0.4459485, 0, 1, 0, 0.5, 0.5, 0 };
const double Q2d6py[] = { 0.09157621, 0.8168476 , 0.09157621, 0.4459485, 0.1081030, 0.4459485, 0, 0, 1, 0,   0.5, 0.5 };
const double Q2d6pw[] = { 0.05497587, 0.05497587, 0.05497587, 0.1116908, 0.1116908, 0.1116908, 1, 1, 1, 1, 1, 1 };
const CQuadrature Q2d6p = { 6, Q2d6px, Q2d6py, 0, Q2d6pw };

// 1D Gaussian quadrature formulas
const double Q1d2px[] = { -sqrt(3)/3, sqrt(3)/3, 0., 1., 0.5 };
const double Q1d2pw[] = { 1., 1., 1., 1., 1. };
const CQuadrature Q1d2p = { 2, Q1d2px, 0, 0, Q1d2pw };

const double Q1d3px[] = { 0., -sqrt(15)/5, sqrt(15)/5, 0., 1., 0.5 };
const double Q1d3pw[] = { 8./9, 5./9, 5./9, 1., 1., 1. };
const CQuadrature Q1d3p = { 3, Q1d3px, 0, 0, Q1d3pw };

const double Q1d4px[] = { -sqrt(525.-70.*sqrt(30.))/35., sqrt(525.-70.*sqrt(30.))/35., -sqrt(525.+70.*sqrt(30.))/35., sqrt(525.+70.*sqrt(30.))/35., 0., 1., 0.5 };
const double Q1d4pw[] = { (18.+sqrt(30.))/36., (18.+sqrt(30.))/36., (18.-sqrt(30.))/36., (18.-sqrt(30.))/36., 1., 1., 1. };
const CQuadrature Q1d4p = { 4, Q1d4px, 0, 0, Q1d4pw };




using namespace std;

/**
  Definition of Integration (interpolation) points
*/
class CIntPoint
{
 protected:
  C3DMesh *imesh;       /**< Associated mesh */
  const CQuadrature *Q; /**< Associated quadrature formula */
  CMat<CVar> J_a,       /**< Jacobian in CVar variables */
             IJ_a,      /**< Inverted jacobian in CVar variables */
	     x_a;       /**< Nodal coordinates in CVar variables */
  CVar IdJ_a,           /**< Determinant of inverted jacobian in CVar variables */
       dJ_a;            /**< Determinant of jacobian in CVar variables */
  CMat<double> SF,      /**< Shape functions */
               GSF,     /**< Gradient of shape functions. */
	       J,       /**< Jacobian */
	       IJ,      /**< Inverted jacobian */
	       x;       /**< Nodal coordinates */
  map<int, map<int, map<int,double> > > HSF; /**< hessian of shape functions */
  double IdJ,           /**< Determinant of inverted jacobian */
         dJ;            /**< Determinant of jacobian */
  CVec dvec;            /**< Auxiliary vector of d.o.f. */
  int *dintvec;         /**< Array of indices of associated d.o.f. */
  int indof;            /**< No. of associated d.o.f. */
  int dim;              /**< Space dimension of integration point */
  
  /**
    Calculate jacobian.
    @param E Type of evaluation
  */
  void calcJ(eval_t E);
  
  /**
    Evaluate function value and gradient at quadrature point.
    
    @param dofvec Vector of d.o.f.
    @param dof    Array of indices of d.o.f.
    @param val    Returned function value
    @param grad   Returned spatial gradient
    @param E      Flag indicating type of evaluation
  */
  virtual void evaluate(CVec &dofvec, int *dof, CVar &val, CVec &grad, eval_t E);
  
  /** Set gradient of shape functions. */
  virtual void setGSF() = 0;
  
  /** Set shape functions. */
  virtual void setSF()  = 0;
  
  /** True if integration point is linear (i.e. gradient of s.f. is constant) */
  virtual bool isLinear() = 0;
  
  
 public:
  double ksi, eta, zeta; /**< Local coordinates */
  CVar dx_a;             /**< Integration weight (in CVar type) */
  double dx;             /**< Integration weight */
  
  /** Constructor (not to be used) */
  CIntPoint() {};
  
  /**
    Default constructor.
  
    Initializes the integration point and associates it with mesh and quadrature.
    @param mesh Associated mesh.
    @param q    Associated quadrature
  */
  CIntPoint(C3DMesh *mesh, CQuadrature *q) { init(mesh, q); };
  
  /** Destructor */
  virtual ~CIntPoint() { delete[] dintvec; };
  
  /** Quadrature getter */
  const CQuadrature *getQ() { return Q; };
  
  /** Getter for dJ */
  double getDJ() { return dJ; };
  
  /**
    Initialization of data structures.
    @param mesh Associated mesh
    @param q    Associated quadrature
    @param ndof Number of d.o.f. per element
    @param dm  Space dimension of integration point
  */
  void         init(C3DMesh *mesh, const CQuadrature *q, int ndof, int dm);
  
  /**
    Virtual method for initialization of particular integration point types.
    @param mesh Associated mesh.
    @param q    Associated quadrature
   */
  virtual void init(C3DMesh *mesh, const CQuadrature *q) {};
  
  /** Return number of d.o.f. per element. */
  int          getNdof() { return indof; }
  
  /**
    Return number of d.o.f. per element in case that mesh
    consists of different element types.
    @param elem Element number
  */
  virtual int  getNdof(int elem) { return indof; }
  
  /** Return number of mesh nodes. */
  virtual int  getMeshNdof() { return imesh->getNnodes(); }
  
  /**
    Virtual method for calculation of element data.
    @param elem Element number
    @param E    Type of evaluation
  */
  virtual void calcE(int elem, eval_t E) = 0;
  
  /**
    Virtual method for calculation of quadrature point data.
    @param qp Quadrature point
    @param E  Type of evaluation
  */
  virtual void calcQP(int qp, eval_t E)  = 0;
  
  /**
    Virtual method for calculation of quadrature point data at arbitrary point.
    @param l1 Local (barycentric) coordinate
    @param l2 Local (barycentric) coordinate
    @param l3 Local (barycentric) coordinate
  */
  virtual void calcQP(double l1, double l2, double l3) {};
  
  /**
    Evaluate base function and spatial gradient in quadrature point.
    @param basis  Number of base function
    @param elem   Element number
    @param offset Offset of current variable within global vector of d.o.f.
    @param val    Returned function value
    @param grad   Returned spatial gradient
    @param E      Type of evaluation
  */
  void         evalBasis(int basis, CVar &val, CVec &grad, eval_t E);
  
  /**
    Virtual method for evaluation of function and gradient in quadrature point.
    @param dof    Global vector of d.o.f.
    @param elem   Element number
    @param offset Offset of current variable within global vector of d.o.f.
    @param val    Returned function value
    @param grad   Returned spatial gradient
    @param E      Type of evaluation
  */
  virtual void eval  (CVec &dof, int elem, int offset, CVar &val, CVec &grad, eval_t E);
  
  /** Return global x-coordinate of actual quadrature point */
  CVar getX();
  /** Return global y-coordinate of actual quadrature point */
  CVar getY();
  /** Return global z-coordinate of actual quadrature point */
  CVar getZ();
  
  /**
    Return global number of d.o.f. given by element and local number.
    @param elem Element number
    @param dof  Local number of d.o.f.
  */
  inline virtual int mapDof (int elem, int dof) { return imesh->getVolumeNode(elem, dof); }; // map local to global index of degree of freedom
};

/**
  Common class for 3D integration points.
*/
class C3DIntPoint : public CIntPoint
{
 protected:
  /**
    Set the coordinates of d.o.f.
    @param elem Element number
    @param E    Type of evaluation
  */
  void setX(int elem, eval_t E);
  
  /**
    Calculate element data.
    @param elem Element number
    @param E    Type of evaluation
  */
  void calcE(int elem, eval_t E);
  
  /**
    Calculate quadrature point data.
    @param qp Quadrature point
    @param E  Type of evaluation
  */
  void calcQP(int qp, eval_t E);
};

/**
  Common class for 2D integration points.
*/
class C2DIntPoint : public CIntPoint
{
 protected:
  /**
    Set the coordinates of d.o.f.
    @param elem Element number
    @param E    Type of evaluation
  */
  void setX(int elem, eval_t E);
  
  /**
    Calculate element data.
    @param elem Element number
    @param E    Type of evaluation
  */
  void calcE(int elem, eval_t E);
  
  /**
    Calculate quadrature point data.
    @param qp Quadrature point
    @param E  Type of evaluation
  */
  void calcQP(int qp, eval_t E);
  
  virtual void setHSF() {};
  void evalHess(CVec &dofvec, int elem, int offset, CMat<CVar> &h, eval_t E);
};

/**
  Common class for 1D integration points.
*/
class C1DIntPoint : public CIntPoint
{
 protected:
  /**
    Set the coordinates of d.o.f.
    @param elem Element number
    @param E    Type of evaluation
  */
  void setX(int elem, eval_t E);
  
  /**
    Calculate element data.
    @param elem Element number
    @param E    Type of evaluation
  */
  void calcE(int elem, eval_t E);
  
  /**
    Calculate quadrature point data.
    @param qp Quadrature point
    @param E  Type of evaluation
  */
  void calcQP(int qp, eval_t E);
};

/** 3D Linear interpolation */
class C3DLinPoint : public C3DIntPoint {
 protected:
  bool isLinear() { return true; };
 public:
  C3DLinPoint() {};
  C3DLinPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 4, 3); };
  void setGSF();
  void setSF();
  void calcQP(double l1, double l2, double l3) {};
  inline int getMeshNdof() { return imesh->getNlnodes(); }
};

/** 3D Linear interpolation on quadratic mesh */
class C3DLinPointOnQuadMesh : public C3DLinPoint
{
 public:
  C3DLinPointOnQuadMesh(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 4, 3); };
  inline int mapDof(int e, int n)
  {
    return ((C3DQuadMesh *)imesh)->getNoOfLnode(imesh->getVolumeNode(e, n));
  };
  inline int getMeshNdof() { return imesh->getNlnodes(); }
};

/** 3D Linear discontinuous interpolation */
class C3DLinDiscPoint : public C3DLinPoint {
 public:
  C3DLinDiscPoint() {};
  C3DLinDiscPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 4, 3); };
  virtual int mapDof(int e, int n) { return e*4+n; };
};

/** 3D Linear + cubic bubble interpolation */
class C3DLinBubblePoint : public C3DIntPoint {
 protected:
  bool isLinear() { return false; };
 public:
  C3DLinBubblePoint() {};
  C3DLinBubblePoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 5, 3); };
  void setGSF();
  void setSF();
  virtual int mapDof(int e, int n);
};

/** 3D Quadratic interpolation */
class C3DQuadPoint : public C3DIntPoint {
 protected:
  bool isLinear() { return false; };
 public:
  C3DQuadPoint() {};
  C3DQuadPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 10, 3); };
  void setGSF();
  void setHSF();
  void setSF();
  void calcQP(double l1, double l2, double l3);
};

/** 2D Linear interpolation */
class C2DLinPoint : public C2DIntPoint {
 protected:
  bool isLinear() { return true; };
 public:
  C2DLinPoint() {};
  C2DLinPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 3, 2); };
  void setGSF();
  void setSF();
  void calcQP(double l1, double l2, double l3) {};
  virtual int mapDof (int elem, int dof) { return imesh->getFaceNode(elem, dof); };
};

/** 2D Linear interpolation on quadratic mesh */
class C2DLinPointOnQuadMesh : public C2DLinPoint
{
 public:
  C2DLinPointOnQuadMesh(C3DMesh *mesh, const CQuadrature *q) : C2DLinPoint(mesh, q) {};
  inline int mapDof(int e, int n)
  {
    return ((C3DQuadMesh *)imesh)->getNoOfLnode(imesh->getFaceNode(e, n));
  };
  int getMeshNdof() { return imesh->getNlnodes(); }
    
};

/** 2D Linear discontinuous interpolation */
class C2DLinDiscPoint : public C2DLinPoint
{
 protected:
  bool isLinear() { return true; };
 public:
  C2DLinDiscPoint() {};
  C2DLinDiscPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 3, 2); };
  virtual int mapDof(int n, int e) { return e*indof+n; }; // map local to global index
};

/** 2D Linear + cubic bubble interpolation */
class C2DLinBubblePoint : public C2DIntPoint {
 protected:
  bool isLinear() { return false; };
 public:
  C2DLinBubblePoint() {};
  C2DLinBubblePoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 4, 2); };
  void setGSF();
  void setSF();
  virtual int mapDof(int n, int e); // map local to global index
  void evalHess(CVec &dof, int elem, int offset, CMat<CVar> &hess);
  void evalHessBasis(int basis, CMat<CVar> &hess);
};

/** 2D non-conforming linear interpolation */
class C2DNCLinPoint : public C2DIntPoint {
 protected:
  bool isLinear() { return true; };
 public:
  C2DNCLinPoint() {};
  C2DNCLinPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 3, 2); };
  void setGSF();
  void setSF();
};

/** 2D Quadratic interpolation */
class C2DQuadPoint: public C2DIntPoint {
 protected:
  bool isLinear() { return false; };
 public:
  C2DQuadPoint() {};
  C2DQuadPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 6, 2); };
  void setGSF();
  void setHSF();
  void setSF();
  void evalHess(CVec &dof, int elem, int offset, CMat<CVar> &hess);
  void evalHessBasis(int basis, CMat<CVar> &hess);
  virtual int mapDof (int elem, int dof) { return imesh->getFaceNode(elem, dof); };
};

/** 2D Quadratic+cubic bubble interpolation */
class C2DQuadBubblePoint : public C2DIntPoint {
 protected:
  bool isLinear() { return false; };
 public:
  C2DQuadBubblePoint() {};
  C2DQuadBubblePoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 7, 2); };
  void setGSF();
  void setSF();
  virtual int mapDof(int n, int e); // map local to global index
};

/** 1D Linear interpolation */
class C1DLinPoint : public C1DIntPoint {
 protected:
  bool isLinear() { return true; };
 public:
  C1DLinPoint() {};
  C1DLinPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 2, 1); };
  void setGSF();
  void setSF();
  void calcQP(double l1, double l2, double l3) {};
  virtual int mapDof(int n, int e); // map local to global index
};

/** 1D Linear interpolation on quadratic mesh */
class C1DLinPointOnQuadMesh : public C1DLinPoint
{
 public:
  C1DLinPointOnQuadMesh(C3DMesh *mesh, const CQuadrature *q) : C1DLinPoint(mesh, q) {};
  inline int mapDof(int e, int n)
  {
    return ((C3DQuadMesh *)imesh)->getNoOfLnode(imesh->getEdgeNode(e, n));
  };
  int getMeshNdof() { return imesh->getNlnodes(); }    
};

/** 1D Quadratic interpolation */
class C1DQuadPoint : public C1DIntPoint
{
 protected:
  bool isLinear() { return false; };
 public:
  C1DQuadPoint() {};
  C1DQuadPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 3, 1); };
  void setGSF();
  void setSF();
  void calcQP(double l1, double l2, double l3) {};
  virtual int mapDof(int n, int e); // map local to global index
};

/** 1D non-conforming linear interpolation */
class C1DNCLinPoint : public C1DIntPoint
{
 protected:
  bool isLinear() { return true; };
 public:
  C1DNCLinPoint() {};
  C1DNCLinPoint(C3DMesh *mesh, const CQuadrature *q) { CIntPoint::init(mesh, q, 3, 1); indof = 1; };
  void setGSF();
  void setSF();
  virtual int mapDof(int n, int e); // map local to global index
};




#endif
