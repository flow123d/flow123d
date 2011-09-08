#ifndef _STATEPROBLEM_H__
#define _STATEPROBLEM_H__

// unknown orders
#define UORDER3DLIN  1
#define UORDER3DQUAD 2

#define UT_SCALAR 1
#define UT_VECTOR 3
#define UT_TENSOR 9

#include <stdio.h>
#include <vector>
#include <map>
#include <ncurses.h>

#include "ad.h"
#include "integration3d.h"
#include "muParser.h"
#include "wmanager.h"

/** Start measuring time. */
void start_clock(clock_t &t);
/** Stop measuring time. */
double stop_clock(clock_t &t);

using namespace mu;

/** Type for derivatives of residual vector w.r.t. degrees of freedom. */
typedef map<int,double> t_drq;

/** Type for map of integration points. */
typedef map<int,CIntPoint*> t_IP;

/** Auxiliary function for calculation of orthogonal vectors. */
void calculateTangentsFromNormal(CVec &n, CVec &t1, CVec &t2);

/** Class for information about components of vector of unknowns (e.g. scalars, vector components). */
class CComponent
{
 public:
  CVar val,     /**< Value at quadrature point.            */
       o_val,   /**< Value from previous timestep.         */
       oo_val;  /**< Value from pre-previous timestep.     */
  CVec grad,    /**< Spatial gradient at quadrature point. */
       o_grad,  /**< Gradient from previous timestep.      */
       oo_grad; /**< Gradient from pre-previous timestep.  */
  CIntPoint *ip3d, *ip2d, *ip1d; /**< Integration points. */
  bool   set_mean_value; /**< Whether to set meanvalue of this component. */
  int    mv_segment;     /**< Physical entity where the meanvalue will be fixed. */
  double mean_value,     /**< Prescribed meanvalue. */
         volume;         /**< Volume (measure) of physical entity mv_segment. */
  double *mean;          /**< Volume (measure) of particular elements. */
  map<string,Parser*> p; /**< Expression parser for Dirichlet boundary condition. */
  double *node_value,    /**< Nodal values (for postprocessing). */
         *node_grad;     /**< Nodal gradients (for postprocessing). */
  
  CComponent();
  ~CComponent();
};

/** Class for information about unknowns (e.g. scalars, vectors, tensors). */
class CUnknown
{
 public:
  string name; /**< Real name of unknown. */
  int type;    /**< Type of unknown (scalar, vector, tensor). */
  int index;   /**< Index of first component within global set of components. */
  map<string,Parser*> normalbc,  /**< Expression parser for Dirichlet b.c. for vector normal components. */
                      tangentbc; /**< Expression parser for Dirichlet b.c. for vector tangent components. */
  
  CUnknown();
  ~CUnknown();
};


/** Class for definition of state problem = variational formulation of system of PDE and FEM approximation. */
class CStateProblem
{
 protected:
  C3DMesh *mesh;                      /**< Associated mesh. */
  const CQuadrature *Q3d, *Q2d, *Q1d; /**< Associated quadrature points. */
  vector<CComponent *> U;             /**< Vector of components. */
  vector<CUnknown *> Uinfo;           /**< Vector of unknowns. */
  t_IP IP3D, IP2D, IP1D;              /**< Integration (interpolation) points. */
  CWManager *wm;                      /**< Window manager. */
  
  CVec q, oq, ooq; /**< Vector of dofs, oq,ooq=from previous time steps */
  double *r;       /**< Residual vector. */
  t_drq *drq;      /**< Derivatives of residual vector w.r.t. state vector (degrees of freedom). */
  
  CVar x, y, z;         /**< Global coordinates. */
  double x_d, y_d, z_d, /**< Global coordinates in CVar type. */
         tau,           /**< Timestep size. */
	 tmax;          /**< Stopping time. */
  CVar phi;             /**< Test function. */
  CVec gphi;            /**< Spatial gradient of test function. */
  
  map<string,string> params; /**< Parameters parsed from input file. */
  bool loaded_solution,      /**< True if solution was loaded from file. */
       calculated_solution;  /**< True if solution has been calculated (for postprocessing). */
  
  FILE *logfile;  /**< File where logs will be saved. */
  map<string, Parser *> PPscalars; /**< Expression parsers for arbitrary postprocessing fields to be saved. */
  
  /**
    Load mesh from file.
    @param meshtype Type of mesh (linear, quadratic).
    @param file     File name.
  */
  void loadMesh(int meshtype, const char *file);
  
  /**
    Add single component (e.g. scalar, vector component).
    @param utype       Order of approximation (linear, quadratic).
    @param _set_mv     Whether to fix meanvalue.
    @param _mv_segment Physical entity where the meanvalue will be fixed.
    @param _mv         Level of meanvalue.
  */
  void addComponent(int utype, bool _set_mv = false, int _mv_segment = 0, double _mv = 0);
  
  /**
    Log progress during computation.
    @param window Window for output.
    @param pct    Percent to be shown.
  */
  void logPct(int window, int pct);
  
  /**
    Add residual and its derivatives.
    @param r_loc Value and derivatives to be added.
    @param ind   Index of residual component.
  */
  void addResidual(CVar &r_loc, int ind);
  
  /**
    Set residual component and its derivatives.
    @param r_loc Value and derivatives to be set.
    @param ind   Index of residual component.
  */
  void setResidual(CVar &r_loc, int ind);
  
  /** Assemble volume integrals */
  void assembleVolumeIntegral();
  
  /** Assemble face integrals. */
  void assembleFaceIntegral();
  
  /** Assemble edge integrals. */
  void assembleEdgeIntegral();
  
  /** Set normal and tangent boundary conditions on faces. */
  void assembleFaceNormalTangentBC();
  
  /** Sert normal and tangent boundary conditions on edges. */
  void assembleEdgeNormalTangentBC();
  
  /** Set Dirichlet b.c. - in the beginning of assembly process. */
  void setDirichletBC();
  
  /** Set residuals for Dirichlet b.c. - at the end of assembly. */
  void setResidualsForDirichletBC();
  
  /** Fix meanvalues. */
  void setMeanValues();
  
  /** Virtual method called after config is readed. */
  virtual void afterReadConfig() {};
  
  /** Virtual method called after data structures are prepared. */
  virtual void afterPrepare() {};
  
  /** Virtual method called after volume components are evaluated. */
  virtual void afterEvaluateVolumeComponents(int elem, int qp) {};
  
  /** Virtual method called after face components are evaluated. */
  virtual void afterEvaluateFaceComponents(int elem, int qp) {};
  
  /** Virtual method called after edge components are evaluated. */
  virtual void afterEvaluateEdgeComponents(int elem, int qp) {};
  
  /** Whether volume integrais will be calculated. */
  virtual bool useVolumeIntegral() = 0;
  
  /** Whether face integrais will be calculated. */
  virtual bool useFaceIntegral()   = 0;
  
  /** Whether edge integrais will be calculated. */
  virtual bool useEdgeIntegral()   = 0;
  
  /**
    Abstract method for definition of volume integrand.
    @param comp Component number.
    @param elem Element number.
    @param qp   Quadrature point.
  */
  virtual CVar volumeIntegral(int comp, int elem, int qp) { return 0; }
  
  /**
    Abstract method for definition of face integrand.
    @param comp Component number.
    @param elem Element number.
  */
  virtual CVar faceIntegral(int comp, int elem) { return 0; }
  
  /**
    Abstract method for definition of edge integrand.
    @param comp Component number.
    @param elem Element number.
  */
  virtual CVar edgeIntegral(int comp, int elem) { return 0; }
  
  /** Evaluate solution at mesh points (for postprocessing). */
  void calculateSolution();
 
 public:
 
  CStateProblem();
  virtual ~CStateProblem();
  
  /**
    Read config file.
    @param file File name.
  */
  void readConfig(char *file);
  
  /**
    Define quadratures.
    @param _q3d 3D quadrature.
    @param _q2d 2D quadrature.
    @param _q1d 1D quadrature.
  */
  void registerQuadratures(const CQuadrature *_q3d, const CQuadrature *_q2d, const CQuadrature *_q1d);
  
  /**
    Define scalar unknown.
    @param name        Real name (e.g. pressure).
    @param utype       Order of approximation (linear, quadratic).
    @param _set_mv     Whether to fix meanvalue.
    @param _mv_segment Physical entity where the meanvalue will be fixed.
    @param _mv         Level of meanvalue.
  */
  void addScalar(string name, int utype, bool _set_mv = false, int _mv_segment = 0, double _mv = 0);
  
  /**
    Define vector unknown.
    @param name        Real name (e.g. velocity).
    @param utype       Order of approximation (linear, quadratic).
    @param _set_mv     Whether to fix meanvalue (of all components).
    @param _mv_segment Physical entity where the meanvalue will be fixed.
    @param _mv         Level of meanvalue.
  */
  void addVector(string name, int utype, bool _set_mv = false, int _mv_segment = 0, double _mv = 0);
  
  /**
    Define tensor unknown.
    @param name        Real name (e.g. stress).
    @param utype       Order of approximation (linear, quadratic).
    @param _set_mv     Whether to fix meanvalue (of all components).
    @param _mv_segment Physical entity where the meanvalue will be fixed.
    @param _mv         Level of meanvalue.
  */
  void addTensor(string name, int utype, bool _set_mv = false, int _mv_segment = 0, double _mv = 0);
  
  /** Calculate size of vector of d.o.f. */
  int  getProblemNdof();
  
  /** Allocate data structures. */
  void prepare();
  
  /** Assemble residual and problem matrix. */
  void assemble();
  
  /**
    Solve linearized system.
    @param relax Relaxation parameter.
  */
  int  solve(double relax = 1);
  
  /**
    Solve linearized time-dependent system.
    @param relax Relaxation parameter.
  */
  void tSolve(double relax = 1);
  
  /**
    Solve nonlinear problem.
    @param relax Relaxation parameter.
  */
  void nlSolve(double relax);
  
  /** Whether problem is time-dependent. */
  bool isTimeDependent() { return tau>0; }
  
  /**
    Save matrix to file.
    @param file File name.
  */
  void saveMatrix(const char *file);
  
  /**
    Add postprocessing scalar field to be saved.
    @param name       Real name (e.g. velocity_magnitude)
    @param expression Expression (e.g..sqrt(velocity0^2+velocity1^2+velocity2^2))
  */
  void addPPscalar(string name, string expression);
  
  /**
    Save solution in VTK format.
    @param file File name.
  */
  void saveVTKview(const char *file = 0);
  
  /**
    Load solution from file.
    @param file File name.
  */
  bool loadSolution(const char *file = 0);
  
  /**
    Save solution to file.
    @param file File name.
  */
  void saveSolution(const char *file = 0);
  
  /**
    Log message.
    @param window Window id.
    @param format Output format (printf-like syntax)
  */
  void log(int window, const char *format, ...);
  
  
};






#endif
