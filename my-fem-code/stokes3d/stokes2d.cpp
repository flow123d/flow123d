#include <math.h>
#include "fem.h"
#include "streamfunction2d.h"


class C2DStokesProblem : public C2DFlowProblem
{
 private:
  int  conv;
  CVar diverg;
 public:
  bool useVolumeIntegral() { return false; }
  bool useFaceIntegral() { return true; }
  bool useEdgeIntegral() { return false; }
  void afterReadConfig();
  void afterPrepare();
  void afterEvaluateFaceComponents(int e, int qp);
  CVar faceIntegral(int comp, int elem);
};


void C2DStokesProblem::afterReadConfig()
{
  registerQuadratures(&Q3d11p, &Q2d6p, &Q1d4p);
  addVector("velocity", UORDER3DQUAD);
  addScalar("pressure", UORDER3DLIN);
  
  C2DFlowProblem::afterReadConfig();
}

void C2DStokesProblem::afterPrepare()
{
  // parameter convection indicates whether convective term is to be included into the equations
  if (params.find("convection") != params.end() && params["convection"].compare("yes") == 0)
  {
    conv = 1;
  }
  else
  {
    conv = 0;
  }
}

void C2DStokesProblem::afterEvaluateFaceComponents(int elem, int qp)
{
  diverg = U[0]->grad(0) + U[1]->grad(1) + U[2]->grad(2);
}

CVar C2DStokesProblem::faceIntegral(int comp, int elem)
{
  switch (comp)
  {
    case 0:
    case 1:
    case 2:
        return (U[comp]->grad|gphi)  // grad v : grad phi
	      - U[3]->val*gphi(comp) // p div phi
	      +(U[comp]->grad(0)*U[0]->val.getVal()
	       +U[comp]->grad(1)*U[1]->val.getVal()
	       +U[comp]->grad(2)*U[2]->val.getVal())*phi*conv; // convective term
      break;
    case 3:
      return phi*diverg; // phi * div v
      break;
  }
  return 0;
}




int main(int argc, char **argv)
{
  if (argc != 2)
  {
    cout << "Usage: " << argv[0] << " problem.conf\n\n";
    exit(1);
  }
  
  C2DStokesProblem P;
  
  P.readConfig(argv[1]);
  
  P.prepare();
  P.nlSolve(1.0);
  P.calculateStreamFunction();
  
  P.saveSolution();
  P.saveVTKview();
  
  return 0;
}


