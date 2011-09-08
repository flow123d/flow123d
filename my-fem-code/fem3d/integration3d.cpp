/* Integration points (quadratures)
   Jukka Toivanen, Jan Stebel
*/

#include <math.h>
#include "ad.h"
#include "integration3d.h"

void CIntPoint::init(C3DMesh *mesh, const CQuadrature *q, int ndof, int dm)
{
  imesh = mesh;
  Q     = q;
  indof = ndof;
  dim   = dm;
  
  J_a. resize(dim, 3); 
  IJ_a.resize(3, dim);
  x_a. resize(ndof, 3);
  
  J. resize(dim, 3); 
  IJ.resize(3, dim);
  x. resize(ndof, 3);
  
  SF. resize(ndof, 1);
  GSF.resize(dim, ndof);
  HSF.clear();
  
  dintvec = new int[indof];
  dvec.resize(indof);
  for (int i=0; i<indof; i++) dintvec[i]=i; 
}

void CIntPoint::eval(CVec &dofvec, int elem, int offset, CVar &val, CVec &grad, eval_t E)
{
  int dof[indof];
  for (int i=0; i<indof; i++)
  {
    dof[i] = (mapDof(elem, i)==-1)?-1:mapDof(elem, i)+offset;
  }
  evaluate(dofvec, dof, val, grad, E);
}

void CIntPoint::evalBasis(int basis, CVar &val, CVec &grad, eval_t E)
{
  dvec = 0.0;
  dvec(basis)= 1.0;
  evaluate(dvec, dintvec, val, grad, E);
}

CVar CIntPoint::getX()
{
  CVar _x = 0;
  for (int i=0; i<indof; i++) _x += x(i,0)*SF(i,0);
  return _x;
}

CVar CIntPoint::getY()
{
  CVar _y = 0;
  for (int i=0; i<indof; i++) _y += x(i,1)*SF(i,0);
  return _y;
}

CVar CIntPoint::getZ()
{
  CVar _z = 0;
  for (int i=0; i<indof; i++) _z += x(i,2)*SF(i,0);
  return _z;
}

void CIntPoint::calcJ(eval_t E)
{
  switch (E)
  {
    case EVAL_F:
    case EVAL_FG:
      for (int i=0; i<dim; i++)
      {
        for (int j=0; j<3; j++)
        {
          J(i,j) = 0;
          for (int k=0; k<indof; k++) J(i,j) += x(k,j)*GSF(i,k);
        }
      }
      break;
    case EVAL_FA:
      for (int i=0; i<dim; i++)
      {
        for (int j=0; j<3; j++)
        {
          J_a(i,j) = 0;
          for (int k=0; k<indof; k++) J_a(i,j) += x_a(k,j)*GSF(i,k);
        }
      }
      break;
  }
}

void CIntPoint::evaluate(CVec &dofvec, int *dof, CVar &val, CVec &grad, eval_t E)
{
  int i, j, k;
  val = 0;
  for (i=0; i<indof; i++) val += dofvec(dof[i])*SF(i,0);
  switch (E)
  {
    case EVAL_F: break;
    case EVAL_FG:
      grad.resize(3);
      for (i=0; i<3; i++)
      {
        grad(i) = 0;
        for (j=0; j<indof; j++)
          for (k=0; k<dim; k++) grad(i) += dofvec(dof[j])*(IJ(i,k)*GSF(k,j));
      }
      break;
    case EVAL_FA:
      grad.resize(3);
      for (i=0; i<3; i++)
      {
        grad(i) = 0;
        for (j=0; j<indof; j++)
          for (k=0; k<dim; k++) grad(i) += dofvec(dof[j])*(IJ_a(i,k)*GSF(k,j));
      }
      break;
  }
}



void C3DIntPoint::setX(int elem, eval_t E)
{
  switch (E)
  {
    case EVAL_F:
    case EVAL_FG:
      for (int i=0; i<indof; i++)
        for (int j=0; j<3; j++)
          x(i,j)=imesh->getCoord(imesh->getVolumeNode(elem,i),j).getVal();
      break;
    case EVAL_FA:
      for (int i=0; i<indof; i++)
        for (int j=0; j<3; j++)
          x_a(i,j)=imesh->getCoord(imesh->getVolumeNode(elem,i),j);
      break;
  }
}

void C3DIntPoint::calcE(int elem, eval_t E)
{
  setX(elem, E);
  if (isLinear())
  {
    switch (E)
    {
      case EVAL_F:
      case EVAL_FG:
        setGSF();
        calcJ(E);
        mat_inverse_3d<double>(J, IJ, dJ, IdJ);
        break;
      case EVAL_FA:
        setGSF();
        calcJ(E);
        mat_inverse_3d<CVar>(J_a, IJ_a, dJ_a, IdJ_a);
    }
  }
}

void C3DIntPoint::calcQP(int qp, eval_t E)
{
  ksi  = Q->x[qp];
  eta  = Q->y[qp];
  zeta = Q->z[qp];
  
  setSF();
  if (!isLinear())
  {
    switch (E)
    {
      case EVAL_F:
      case EVAL_FG:
        setGSF();
        calcJ(E);
        mat_inverse_3d<double>(J, IJ, dJ, IdJ);
        dx = dJ*Q->w[qp];
        break;
      case EVAL_FA:
        setGSF();
        calcJ(E);
        mat_inverse_3d<CVar>(J_a, IJ_a, dJ_a, IdJ_a);
        dx_a = dJ_a*Q->w[qp];
    }
  }
  else if (E == EVAL_FA)
  {
    dx_a = dJ_a*Q->w[qp];
  }
  else
  {
    dx = dJ*Q->w[qp];
  }
}


void C2DIntPoint::setX(int elem, eval_t E)
{
  switch (E)
  {
    case EVAL_F:
    case EVAL_FG:
      for (int i=0; i<indof; i++)
        for (int j=0; j<3; j++)
          x(i,j)=imesh->getCoord(imesh->getFaceNode(elem,i),j).getVal();
      break;
    case EVAL_FA:
      for (int i=0; i<indof; i++)
        for (int j=0; j<3; j++)
          x_a(i,j)=imesh->getCoord(imesh->getFaceNode(elem,i),j);
      break;
  }
}

/*void C2DIntPoint::calcJ(eval_t E)
{
  switch (E)
  {
    case EVAL_F:
    case EVAL_FG:
      for (int i=0; i<2; i++)
      {
        for (int j=0; j<3; j++)
        {
          J(i,j) = 0;
          for (int k=0; k<indof; k++) J(i,j) += x(k,j)*GSF(i,k);
        }
      }
      break;
    case EVAL_FA:
      for (int i=0; i<2; i++)
      {
        for (int j=0; j<3; j++)
        {
          J_a(i,j) = 0;
          for (int k=0; k<indof; k++) J_a(i,j) += x_a(k,j)*GSF(i,k);
        }
      }
      break;
  }
}*/

void C2DIntPoint::calcE(int elem, eval_t E)
{
  setX(elem, E);
  if (isLinear())
  {
    switch (E)
    {
      case EVAL_F:
      case EVAL_FG:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_2x3<double>(J, IJ, dJ, IdJ);
        break;
      case EVAL_FA:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_2x3<CVar>(J_a, IJ_a, dJ_a, IdJ_a);
	break;
    }
  }
}

void C2DIntPoint::calcQP(int qp, eval_t E)
{
  ksi  = Q->x[qp];
  eta  = Q->y[qp];
  
  setSF();
  if (!isLinear())
  {
    switch (E)
    {
      case EVAL_F:
      case EVAL_FG:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_2x3<double>(J, IJ, dJ, IdJ);
        dx = dJ*Q->w[qp];
        break;
      case EVAL_FA:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_2x3<CVar>(J_a, IJ_a, dJ_a, IdJ_a);
        dx_a = dJ_a*Q->w[qp];
	break;
    }
  }
  else if (E == EVAL_FA)
  {
    dx_a = dJ_a*Q->w[qp];
  }
  else
  {
    dx = dJ*Q->w[qp];
  }
}

/*void C2DIntPoint::evaluate(CVec &dofvec, int *dof, CVar &val, CVec &grad, eval_t E)
{
  int i, j, k;
  val = 0;
  for (i=0; i<indof; i++) val += dofvec(dof[i])*SF(i,0);
  switch (E)
  {
    case EVAL_F: break;
    case EVAL_FG:
      grad.resize(3);
      for (i=0; i<3; i++)
      {
        grad(i) = 0;
        for (j=0; j<indof; j++)
          for (k=0; k<2; k++) grad(i) += dofvec(dof[j])*(IJ(i,k)*GSF(k,j));
      }
      break;
    case EVAL_FA:
      grad.resize(3);
      for (i=0; i<3; i++)
      {
        grad(i) = 0;
        for (j=0; j<indof; j++)
          for (k=0; k<2; k++) grad(i) += dofvec(dof[j])*(IJ_a(i,k)*GSF(k,j));
      }
      break;
  }
}*/

void C2DIntPoint::evalHess(CVec &dofvec, int elem, int offset, CMat<CVar> &h, eval_t E)
{
  int i, j, k, l, m;
  int dof[indof];
  
  for (i=0; i<indof; i++)
  {
    dof[i] = (mapDof(elem, i)==-1)?-1:mapDof(elem, i)+offset;
  }
  
  h.resize(3,3);
  
  switch (E)
  {
    case EVAL_F: break;
    case EVAL_FG:
      for (i=0; i<3; i++)
      {
        for (j=0; j<3; j++)
        {
          h(i,j) = 0;
          for (map<int,map<int,map<int,double> > >::iterator ik=HSF.begin(); ik!=HSF.end(); ik++)
          {
            k = ik->first;
            for (map<int, map<int,double> >::iterator il=ik->second.begin(); il!=ik->second.end(); il++)
    	    {
    	      l = il->first;
    	      for (map<int,double>::iterator im=il->second.begin(); im!=il->second.end(); im++)
    	      {
    	        m = im->first;
	        h(i,j) = dofvec(dof[m])*IJ(i,k)*im->second*IJ(l,j);
	      }
	    }
          }
        }
      }
      break;
    case EVAL_FA:
      for (i=0; i<3; i++)
      {
        for (j=0; j<3; j++)
        {
          h(i,j) = 0;
          for (map<int,map<int,map<int,double> > >::iterator ik=HSF.begin(); ik!=HSF.end(); ik++)
          {
            k = ik->first;
            for (map<int, map<int,double> >::iterator il=ik->second.begin(); il!=ik->second.end(); il++)
    	    {
    	      l = il->first;
    	      for (map<int,double>::iterator im=il->second.begin(); im!=il->second.end(); im++)
    	      {
    	        m = im->first;
	        h(i,j) = IJ_a(i,k)*IJ_a(l,j)*dofvec(dof[m])*im->second;
	      }
	    }
          }
        }
      }
      break;
  }
}

void C1DIntPoint::setX(int elem, eval_t E)
{
  switch (E)
  {
    case EVAL_F:
    case EVAL_FG:
      for (int i=0; i<indof; i++)
        for (int j=0; j<3; j++)
          x(i,j)=imesh->getCoord(imesh->getEdgeNode(elem,i),j).getVal();
      break;
    case EVAL_FA:
      for (int i=0; i<indof; i++)
        for (int j=0; j<3; j++)
          x_a(i,j)=imesh->getCoord(imesh->getEdgeNode(elem,i),j);
      break;
  }
}

/*void C1DIntPoint::calcJ(eval_t E)
{
  switch (E)
  {
    case EVAL_F:
    case EVAL_FG:
      for (int i=0; i<1; i++)
      {
        for (int j=0; j<3; j++)
        {
          J(i,j) = 0;
          for (int k=0; k<indof; k++) J(i,j) += x(k,j)*GSF(i,k);
        }
      }
      break;
    case EVAL_FA:
      for (int i=0; i<1; i++)
      {
        for (int j=0; j<3; j++)
        {
          J_a(i,j) = 0;
          for (int k=0; k<indof; k++) J_a(i,j) += x_a(k,j)*GSF(i,k);
        }
      }
      break;
  }
}*/

void C1DIntPoint::calcE(int elem, eval_t E)
{
  setX(elem, E);
  if (isLinear())
  {
    switch (E)
    {
      case EVAL_F:
      case EVAL_FG:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_1x3<double>(J, IJ, dJ, IdJ);
        break;
      case EVAL_FA:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_1x3<CVar>(J_a, IJ_a, dJ_a, IdJ_a);
    }
  }
}

void C1DIntPoint::calcQP(int qp, eval_t E)
{
  ksi  = Q->x[qp];
  
  setSF();
  if (!isLinear())
  {
    switch (E)
    {
      case EVAL_F:
      case EVAL_FG:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_1x3<double>(J, IJ, dJ, IdJ);
        dx = dJ*Q->w[qp];
        break;
      case EVAL_FA:
        setGSF();
        calcJ(E);
        mat_pseudoinverse_1x3<CVar>(J_a, IJ_a, dJ_a, IdJ_a);
        dx_a = dJ_a*Q->w[qp];
	break;
    }
  }
  else if (E == EVAL_FA)
  {
    dx_a = dJ_a*Q->w[qp];
  }
  else
  {
    dx = dJ*Q->w[qp];
  }
}

/*void C1DIntPoint::evaluate(CVec &dofvec, int *dof, CVar &val, CVec &grad, eval_t E)
{
  int i, j, k;
  val = 0;
  for (i=0; i<indof; i++) val += dofvec(dof[i])*SF(i,0);
  switch (E)
  {
    case EVAL_F: break;
    case EVAL_FG:
      grad.resize(3);
      for (i=0; i<3; i++)
      {
        grad(i) = 0;
        for (j=0; j<indof; j++)
          for (k=0; k<1; k++) grad(i) += dofvec(dof[j])*(IJ(i,k)*GSF(k,j));
      }
      break;
    case EVAL_FA:
      grad.resize(3);
      for (i=0; i<3; i++)
      {
        grad(i) = 0;
        for (j=0; j<indof; j++)
          for (k=0; k<1; k++) grad(i) += dofvec(dof[j])*(IJ_a(i,k)*GSF(k,j));
      }
      break;
  }
}*/




void C3DLinPoint::setGSF()
{
  GSF(0,0) = -1;
  GSF(0,1) =  1;
  GSF(0,2) =  0;
  GSF(0,3) =  0;
  GSF(1,0) = -1;
  GSF(1,1) =  0;
  GSF(1,2) =  1;
  GSF(1,3) =  0;
  GSF(2,0) = -1;
  GSF(2,1) =  0;
  GSF(2,2) =  0;
  GSF(2,3) =  1;
}

void C3DLinPoint::setSF()
{
  SF(0,0) = 1-ksi-eta-zeta;
  SF(1,0) = ksi;
  SF(2,0) = eta;
  SF(3,0) = zeta;
}



void C3DLinBubblePoint::setSF()
{
  SF(0,0) = 1.0-ksi-eta-zeta;
  SF(1,0) = ksi;
  SF(2,0) = eta;
  SF(3,0) = zeta;
  SF(4,0) = ksi*eta*zeta;
}

void C3DLinBubblePoint::setGSF()
{
  GSF(0,0) = -1;
  GSF(0,1) = 1;
  GSF(0,2) = 0;
  GSF(0,3) = 0;
  GSF(0,4) = eta*zeta;
  GSF(1,0) = -1;
  GSF(1,1) = 0;
  GSF(1,2) = 1;
  GSF(1,3) = 0;
  GSF(1,4) = ksi*zeta;
  GSF(2,0) = -1;
  GSF(2,1) = 0;
  GSF(2,2) = 0;
  GSF(2,3) = 1;
  GSF(2,4) = -ksi*eta;
}

int C3DLinBubblePoint::mapDof(int elem, int n)
{
  return (n==4)?(imesh->getNnodes()+elem):(imesh->getVolumeNode(elem, n));
}



void C3DQuadPoint::setGSF()
{
  GSF(0,0) = -3 + 4*(ksi+eta+zeta);
  GSF(1,0) = -3 + 4*(ksi+eta+zeta);
  GSF(2,0) = -3 + 4*(ksi+eta+zeta);
  GSF(0,1) = -1 + 4*ksi;
  GSF(1,1) = 0;
  GSF(2,1) = 0;
  GSF(0,2) = 0;
  GSF(1,2) = -1 + 4*eta;
  GSF(2,2) = 0;
  GSF(0,3) = 0;
  GSF(1,3) = 0;
  GSF(2,3) = -1 + 4*zeta;
  GSF(0,4) = 4*(1-2*ksi-eta-zeta);
  GSF(1,4) = -4*ksi;
  GSF(2,4) = -4*ksi;
  GSF(0,5) = 4*eta;
  GSF(1,5) = 4*ksi;
  GSF(2,5) = 0;
  GSF(0,6) = -4*eta;
  GSF(1,6) = 4*(1-ksi-2*eta-zeta);
  GSF(2,6) = -4*eta;
  GSF(0,7) = -4*zeta;
  GSF(1,7) = -4*zeta;
  GSF(2,7) = 4*(1-ksi-eta-2*zeta);
  GSF(0,8) = 0;
  GSF(1,8) = 4*zeta;
  GSF(2,8) = 4*eta;
  GSF(0,9) = 4*zeta;
  GSF(1,9) = 0;
  GSF(2,9) = 4*ksi;
}

void C3DQuadPoint::setHSF()
{
  HSF[0][0][0] = 4;
  HSF[0][1][0] = 4;
  HSF[0][2][0] = 4;
  HSF[1][1][0] = 4;
  HSF[1][2][0] = 4;
  HSF[2][2][0] = 4;
  HSF[0][0][1] = 4;
  HSF[1][1][2] = 4;
  HSF[2][2][3] = 4;
  HSF[0][0][4] = -8;
  HSF[0][1][4] = -4;
  HSF[0][2][4] = -4;
  HSF[0][1][5] = 4;
  HSF[0][1][6] = -4;
  HSF[1][1][6] = -8;
  HSF[1][2][6] = -4;
  HSF[0][2][7] = -4;
  HSF[1][2][7] = -4;
  HSF[2][2][7] = -8;
  HSF[1][2][8] = 4;
  HSF[0][2][9] = 4;
}

void C3DQuadPoint::setSF()
{
  SF(0,0) = 1-3*ksi-3*eta-3*zeta+2*ksi*ksi+4*ksi*eta+4*ksi*zeta+2*eta*eta+4*eta*zeta+2*zeta*zeta;
  SF(1,0) = -ksi+2*ksi*ksi;
  SF(2,0) = -eta+2*eta*eta;
  SF(3,0) = -zeta+2*zeta*zeta;
  SF(4,0) = 4*ksi*(1-ksi-eta-zeta);
  SF(5,0) = 4*ksi*eta;
  SF(6,0) = 4*eta*(1-ksi-eta-zeta);
  SF(7,0) = 4*zeta*(1-ksi-eta-zeta);
  SF(8,0) = 4*eta*zeta;
  SF(9,0) = 4*ksi*zeta;
}

void C3DQuadPoint::calcQP(double l1, double l2, double l3)
{
  ksi  = l1;
  eta  = l2;
  zeta = l3;
  
  setGSF();
  setSF();
  calcJ(EVAL_FG);

  mat_inverse_3d<double>(J, IJ, dJ, IdJ);
  
//  dx   = dJ*Q->w[qp];
}



void C2DLinPoint::setSF()
{
  SF(0,0) = 1.0-ksi-eta;
  SF(1,0) = ksi;
  SF(2,0) = eta;
}

void C2DLinPoint::setGSF()
{
  GSF(0,0) = -1;
  GSF(0,1) = 1;
  GSF(0,2) = 0;
  GSF(1,0) = -1;
  GSF(1,1) = 0;
  GSF(1,2) = 1;
}






void C2DLinBubblePoint::setSF()
{
  SF(0,0) = 1.0-ksi-eta;
  SF(1,0) = ksi;
  SF(2,0) = eta;
  SF(3,0) = ksi*eta*(1-ksi-eta);
}

void C2DLinBubblePoint::setGSF()
{
  GSF(0,0) = -1;
  GSF(0,1) = 1;
  GSF(0,2) = 0;
  GSF(0,3) = eta*(1-2*ksi-eta);
  GSF(1,0) = -1;
  GSF(1,1) = 0;
  GSF(1,2) = 1;
  GSF(1,3) = ksi*(1-ksi-2*eta);
}

int C2DLinBubblePoint::mapDof(int elem, int n)
{
  return (n==3)?(imesh->getNlnodes()+elem):(imesh->getFaceNode(elem, n));
}

void C2DLinBubblePoint::evalHess(CVec &dof, int elem, int offset, CMat<CVar> &hess)
{
  int d;
  CMat<CVar> h(2,2);
  hess.resize(2,2);
  
  d = mapDof(elem, 3) + offset;
  
  h(0,0) = dof(d)*(-2.0*eta);
  h(1,0) = dof(d)*(1.0-2.0*ksi-2.0*eta);
  h(0,1) = h(1,0);
  h(1,1) = dof(d)*(-2.0*ksi);
  
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      hess(i,j) = 0;
      for (int k=0; k<2; k++) {
        for (int l=0; l<2; l++) {
          hess(i,j) = hess(i,j) + h(k,l)*IJ(i,k)*IJ(j,l);
	}
      }
    }
  }
}

void C2DLinBubblePoint::evalHessBasis(int basis, CMat<CVar> &hess)
{
  CMat<CVar> h(2,2);
  hess.resize(2,2);

  if (basis == 3) {
    h(0,0) = (-2.0*eta);
    h(1,0) = (1.0-2.0*ksi-2.0*eta);
    h(0,1) = h(1,0);
    h(1,1) = (-2.0*ksi);
  
    for (int i=0; i<2; i++) {
      for (int j=0; j<2; j++) {
        hess(i,j) = 0;
        for (int k=0; k<2; k++) {
          for (int l=0; l<2; l++) {
            hess(i,j) = hess(i,j) + h(k,l)*IJ(k,i)*IJ(l,j);
	  }
        }
      }
    }
  } else {
    hess(0,0) = 0.0; hess(0,1) = 0.0; hess(1,0) = 0.0; hess(1,1) = 0.0;
  }
}



void C2DNCLinPoint::setSF()
{
  SF(0,0) = 1-2*eta;
  SF(1,0) = -1+2*ksi+2*eta;
  SF(2,0) = 1-2*ksi;
}

void C2DNCLinPoint::setGSF()
{
  GSF(0,0) = 0;
  GSF(0,1) = 2;
  GSF(0,2) = -2;
  GSF(1,0) = -2;
  GSF(1,1) = 2;
  GSF(1,2) = 0;
}





void C2DQuadPoint::setSF()
{
  SF(0,0) = (1.0-ksi-eta)*(1.0-2.0*ksi-2.0*eta);
  SF(1,0) = ksi*(2.0*ksi-1.0);
  SF(2,0) = eta*(2.0*eta-1.0);
  SF(3,0) = 4.0*ksi*(1.0-ksi-eta);
  SF(4,0) = 4.0*ksi*eta;
  SF(5,0) = 4.0*eta*(1.0-ksi-eta);
}

void C2DQuadPoint::setGSF()
{
  GSF(0,0) = -3+4*ksi+4*eta;
  GSF(0,1) = 4*ksi-1;
  GSF(0,2) = 0;
  GSF(0,3) = 4*(1-2*ksi-eta);
  GSF(0,4) = 4*eta;
  GSF(0,5) = -4*eta;
  GSF(1,0) = -3+4*ksi+4*eta;
  GSF(1,1) = 0;
  GSF(1,2) = 4*eta-1;
  GSF(1,3) = -4*ksi;
  GSF(1,4) = 4*ksi;
  GSF(1,5) = 4*(1-ksi-2*eta);
}

void C2DQuadPoint::setHSF()
{
  HSF[0][0][0] = 4;
  HSF[0][1][0] = 4;
  HSF[1][1][0] = 4;
  HSF[0][0][1] = 4;
  HSF[1][1][2] = 4;
  HSF[0][0][3] = -8;
  HSF[0][1][3] = -4;
  HSF[0][1][4] = 4;
  HSF[0][1][5] = -4;
  HSF[1][1][5] = -8;
}

void C2DQuadPoint::evalHess(CVec &dof, int elem, int offset, CMat<CVar> &hess) {
  int d[6];
  CMat<CVar> h(2,2);
  hess.resize(2,2);
  
  for (int i=0; i<6; i++) {
    d[i] = imesh->getFaceNode(elem,i) + offset;
  }
  
  h(0,0) = (CVar)4*(dof(d[0])+dof(d[1])-(CVar)2*dof(d[5]));
  h(1,0) = (CVar)4*((CVar)2*dof(d[0])-dof(d[5])+dof(d[3])-dof(d[4]));
  h(0,1) = h(1,0);
  h(1,1) = (CVar)4*(dof(d[0])+dof(d[2])-(CVar)2*dof(d[4]));
  
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      hess(i,j) = 0;
      for (int k=0; k<2; k++) {
        for (int l=0; l<2; l++) {
          hess(i,j) = hess(i,j) + h(k,l)*IJ(i,k)*IJ(j,l);
	}
      }
    }
  }
}

void C2DQuadPoint::evalHessBasis(int basis, CMat<CVar> &hess) {
  CMat<CVar> h(2,2);
  hess.resize(2,2);
  dvec=0.0;
  dvec(basis) = 1.0;
  
  h(0,0) = (CVar)4*(dvec(0)+dvec(1)-(CVar)2*dvec(5));
  h(1,0) = (CVar)4*((CVar)2*dvec(0)-dvec(5)+dvec(3)-dvec(4));
  h(0,1) = h(1,0);
  h(1,1) = (CVar)4*(dvec(0)+dvec(2)-(CVar)2*dvec(4));
  
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      hess(i,j) = 0;
      for (int k=0; k<2; k++) {
        for (int l=0; l<2; l++) {
          hess(i,j) = hess(i,j) + h(k,l)*IJ(k,i)*IJ(l,j);
	}
      }
    }
  }
}




void C2DQuadBubblePoint::setSF()
{
  SF(0,0) = (1.0-ksi-eta)*(1.0-2.0*ksi-2.0*eta);
  SF(1,0) = ksi*(2.0*ksi-1.0);
  SF(2,0) = eta*(2.0*eta-1.0);
  SF(3,0) = 4.0*ksi*(1.0-ksi-eta);
  SF(4,0) = 4.0*ksi*eta;
  SF(5,0) = 4.0*eta*(1.0-ksi-eta);
  SF(6,0) = ksi*eta*(1-ksi-eta);
}

void C2DQuadBubblePoint::setGSF()
{
  GSF(0,0) = -3+4*ksi+4*eta;
  GSF(0,1) = 4*ksi-1;
  GSF(0,2) = 0;
  GSF(0,3) = 4*(1-2*ksi-eta);
  GSF(0,4) = 4*eta;
  GSF(0,5) = -4*eta;
  GSF(0,6) = eta*(1-2*ksi-eta);
  GSF(1,0) = -3+4*ksi+4*eta;
  GSF(1,1) = 0;
  GSF(1,2) = 4*eta-1;
  GSF(1,3) = -4*ksi;
  GSF(1,4) = 4*ksi;
  GSF(1,5) = 4*(1-ksi-2*eta);
  GSF(1,6) = ksi*(1-ksi-2*eta);
}

int C2DQuadBubblePoint::mapDof(int elem, int n)
{
  return (n==6)?(imesh->getNnodes()+elem):(imesh->getFaceNode(elem, n));
}




void C1DLinPoint::setSF()
{
  SF(0,0) = 0.5*(1.0-ksi);
  SF(1,0) = 0.5*(1.0+ksi);
}

void C1DLinPoint::setGSF()
{
  GSF(0,0) = -0.5;
  GSF(0,1) = 0.5;
}

int C1DLinPoint::mapDof(int be, int n)
{
  return imesh->getEdgeNode(be, n);
}





void C1DQuadPoint::setSF()
{
  SF(0,0) = 0.5*ksi*(ksi-1);
  SF(1,0) = 0.5*ksi*(ksi+1);
  SF(2,0) = -ksi*ksi+1;
}

void C1DQuadPoint::setGSF()
{
  GSF(0,0) = ksi-0.5;
  GSF(0,1) = ksi+0.5;
  GSF(0,2) = -2*ksi;
}

int C1DQuadPoint::mapDof(int be, int n)
{
  return imesh->getEdgeNode(be, n);
}




void C1DNCLinPoint::setSF()
{
  SF(0,0) = 1;
  SF(1,0) = ksi;
  SF(2,0) = -ksi;
}

void C1DNCLinPoint::setGSF()
{
  GSF(0,0) = 0;
  GSF(0,1) = 1;
  GSF(0,2) = -1;
}

int C1DNCLinPoint::mapDof(int be, int n)
{
  return imesh->getEdgeNode(be, n);
}










