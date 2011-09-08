#include "fem.h"

class C3DStokesProblem : public CStateProblem
{
 private:
  CVar div;
 public:
  bool useVolumeIntegral() { return true; };
  bool useFaceIntegral() { return false; };
  bool useEdgeIntegral() { return false; };
  void afterEvaluateVolumeComponents(int e, int qp);
  CVar volumeIntegral(int comp, int elem, int qp);
};

void C3DStokesProblem::afterEvaluateVolumeComponents(int e, int qp)
{
  // calculate divergence of velocity
  div = U[0]->grad(0) + U[1]->grad(1) + U[2]->grad(2);
}

CVar C3DStokesProblem::volumeIntegral(int comp, int elem, int qp)
{
  switch (comp)
  {
    case 0:
    case 1:
    case 2:
        return (U[comp]->grad|gphi) - U[3]->val*gphi(comp);
	break;
    case 3:
      return phi*div;
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

  C3DStokesProblem P;
  
  P.readConfig(argv[1]);
  P.registerQuadratures(&Q3d11p, &Q2d4p, &Q1d4p);
  P.addVector("velocity", UORDER3DQUAD);
  P.addScalar("pressure", UORDER3DLIN);
  
  P.prepare();
  P.assemble();
  P.solve(1.0);
  
  P.saveSolution();
  P.saveVTKview();
  
  return 0;
}


