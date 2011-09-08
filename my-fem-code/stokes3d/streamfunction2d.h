/*
  class for copmutation of stream function of 2D velocity field
*/
#include "fem.h"

class C2DFlowProblem : public CStateProblem
{
 public:
  virtual void afterReadConfig();
  void calculateStreamFunction();
};

class C2DStreamFunctionProblem : public CStateProblem
{
 private:
  CVec *vq;
  CVar vx, vy;
 public:
  void afterEvaluateFaceComponents(int elem, int qp);
  CVec &getQ() { return q; };
  void setVelocityDOFs(CVec *_vq) { vq = _vq; };
  void setMesh(C3DMesh *m) { mesh = m; };
  bool useVolumeIntegral() { return false; };
  bool useFaceIntegral() { return true; };
  bool useEdgeIntegral() { return false; };
  CVar faceIntegral(int comp, int elem);
};


C2DStreamFunctionProblem SFP;

double getStreamF(double node)
{
  return SFP.getQ()((int)node).getVal();
}


void C2DFlowProblem::afterReadConfig()
{
  Parser *p = new Parser;
  p->SetExpr("StreamF(node)");
  p->DefineFun("StreamF", getStreamF, false);
  PPscalars["streamfunction"] = p;
  
  SFP.setMesh(mesh);
  SFP.registerQuadratures(Q3d, Q2d, Q1d);
  SFP.setVelocityDOFs(&q);
  SFP.addScalar("streamfunction", UORDER3DQUAD);
}

void C2DFlowProblem::calculateStreamFunction()
{
  SFP.prepare();
  SFP.assemble();
  SFP.solve(1.0);
}




void C2DStreamFunctionProblem::afterEvaluateFaceComponents(int elem, int qp)
{
  CVec g;
  U[0]->ip2d->eval(*vq, elem, 0, vx, g, EVAL_FG);
  U[0]->ip2d->eval(*vq, elem, mesh->getNnodes(), vy, g, EVAL_FG);
}

CVar C2DStreamFunctionProblem::faceIntegral(int comp, int elem)
{
  return (U[0]->grad|gphi) - (vx*gphi(1) - vy*gphi(0));
}


