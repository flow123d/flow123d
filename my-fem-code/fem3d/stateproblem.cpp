#include <stdlib.h>
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "stateproblem.h"
#include "solution.h"



void start_clock(clock_t &t)
{
  t = clock();
}

double stop_clock(clock_t &t)
{
  return ((double)(clock()-t))/CLOCKS_PER_SEC;
}

CComponent::CComponent()
{
}

CComponent::~CComponent()
{
  for (map<string,Parser*>::iterator pit=p.begin(); pit!=p.end(); pit++)
  {
    delete pit->second;
  }
  p.clear();
    
  delete[] node_value;
  delete[] node_grad;
}

CUnknown::CUnknown()
{
}

CUnknown::~CUnknown()
{
  for (map<string,Parser*>::iterator pit=normalbc.begin(); pit!=normalbc.end(); pit++) delete pit->second;
  normalbc.clear();
  for (map<string,Parser*>::iterator pit=tangentbc.begin(); pit!=tangentbc.end(); pit++) delete pit->second;
  tangentbc.clear();
}


CStateProblem::CStateProblem() :
  mesh(0),
  Q3d(0),
  Q2d(0),
  wm(0),
  tau(0),
  tmax(0),
  loaded_solution(0),
  calculated_solution(0),
  logfile(0)
{
  IP3D.clear();
  IP2D.clear();
  IP1D.clear();
  U.clear();
  Uinfo.clear();
  PPscalars.clear();
  
  wm = CWManager::getInstance();
}


CStateProblem::~CStateProblem()
{
  t_IP::iterator it;
  for (it=IP3D.begin(); it != IP3D.end(); it++) delete (*it).second;
  IP3D.clear();
  for (it=IP2D.begin(); it != IP2D.end(); it++) delete (*it).second;
  IP2D.clear();
  for (it=IP1D.begin(); it != IP1D.end(); it++) delete (*it).second;
  IP1D.clear();
  
  U.clear();
  
  for (int i=0; i<(int)Uinfo.size(); i++) delete Uinfo[i];
  Uinfo.clear();
  
  if (mesh == 0)
  {
    delete mesh;
    mesh = 0;
  }
  
  for (int i=0; i<q.size(); i++) drq[i].clear();
  delete[] drq;
  
  delete[] r;
  
  params.clear();
  
  for (map<string,Parser *>::iterator it = PPscalars.begin(); it != PPscalars.end(); it++) delete it->second;
  PPscalars.clear();
  
  if (CWManager::getInstanceNotCreate() == wm && wm != 0) delete wm;
  
  if (logfile != 0) fclose(logfile);
}

void CStateProblem::readConfig(char *file)
{
  ifstream f(file);
  string l, key, data, data1, data2;
  
  // check if file exists
  if ((bool)f == false) return;
  
  log(W_PROG, "Reading config file %s.\n", file);
  
  while (!f.eof())
  {
    getline(f, l);
    key = l.substr(0, l.find_first_of('='));
    data = l.substr(l.find_first_of('=')+1);
    
    if (key.compare("meshfile") == 0)
    {
      data1 = data.substr(0, data.find_first_of(','));
      data2 = data.substr(data.find_first_of(',')+1);
      if (data2.compare("linear") == 0)
      {
        loadMesh(UORDER3DLIN, data1.c_str());
      }
      else if (data2.compare("quadratic") == 0)
      {
        loadMesh(UORDER3DQUAD, data1.c_str());
      }
      else
      {
        log(W_PROG, "Error reading config file %s: Unsupported mesh type (%s).\n", file, data2.c_str());
	throw 1;
	exit(1);
      }
    }
    else if (key.compare("logfile") == 0)
    { // open log file
      logfile = fopen(data.c_str(), "w");
    }
    else if (key.compare("timestep") == 0)
    { // open log file
      tau = atof(data.c_str());
    }
    else if (key.compare("stoptime") == 0)
    { // open log file
      tmax = atof(data.c_str());
    }
    else if (key.compare(0, 11, "save_scalar") == 0)
    { // save post process field
      string name = key.substr(12, key.size()-13);
      addPPscalar(name, data);
    }
    else if (key.compare("") != 0 && key.compare(0, 1, "#") != 0)
    { // unrecognized keys are stored
      params[key] = data;
    }
  }
  afterReadConfig();
}

void CStateProblem::loadMesh(int meshtype, const char *file)
{
  switch (meshtype)
  {
    case UORDER3DLIN:
      mesh = new C3DLinMesh;
      break;
    case UORDER3DQUAD:
      mesh = new C3DQuadMesh;
      break;
    default:
      log(W_PROG, "Error: Unknown mesh type.\n");
      throw 1;
      exit(1);
  }
  mesh->readGMSH(file);
  log(W_PROG, "Loaded mesh file %s (%d volumes, %d faces, %d edges, %d nodes).\n", file, mesh->getNvolumes(), mesh->getNfaces(), mesh->getNedges(), mesh->getNnodes());
}



void CStateProblem::registerQuadratures(const CQuadrature *_q3d, const CQuadrature *_q2d, const CQuadrature *_q1d)
{
  Q3d = _q3d;
  Q2d = _q2d;
  Q1d = _q1d;
}

void CStateProblem::addComponent(int utype, bool _set_mv, int _mv_segment, double _mv)
{
  CComponent *u = new CComponent;
  t_IP::iterator i;
  Parser *p;
  char no[5];
  
  switch (utype)
  {
    case UORDER3DLIN:
      if (IP3D.find(UORDER3DLIN) == IP3D.end())
      {
        if (dynamic_cast<C3DLinMesh *>(mesh) != 0)
        {
	  IP3D[UORDER3DLIN] = new C3DLinPoint(mesh, Q3d);
	  IP2D[UORDER3DLIN] = new C2DLinPoint(mesh, Q2d);
	  IP1D[UORDER3DLIN] = new C1DLinPoint(mesh, Q1d);
	}
	else if (dynamic_cast<C3DQuadMesh *>(mesh) != 0)
	{
	  IP3D[UORDER3DLIN] = new C3DLinPointOnQuadMesh(mesh, Q3d);
	  IP2D[UORDER3DLIN] = new C2DLinPointOnQuadMesh(mesh, Q2d);
	  IP1D[UORDER3DLIN] = new C1DLinPointOnQuadMesh(mesh, Q1d);
	}
	else
	{
	  log(W_PROG, "Error: Component type is incompatible with mesh.\n");
	  throw 1;
	  exit(1);
	}
      }
      u->ip3d = IP3D[UORDER3DLIN];
      u->ip2d = IP2D[UORDER3DLIN];
      u->ip1d = IP1D[UORDER3DLIN];
      break;
    case UORDER3DQUAD:
      if (IP3D.find(UORDER3DQUAD) == IP3D.end())
      {
        IP3D[UORDER3DQUAD] = new C3DQuadPoint(mesh, Q3d);
	IP2D[UORDER3DQUAD] = new C2DQuadPoint(mesh, Q2d);
	IP1D[UORDER3DQUAD] = new C1DQuadPoint(mesh, Q1d);
      }
      u->ip3d = IP3D[UORDER3DQUAD];
      u->ip2d = IP2D[UORDER3DQUAD];
      u->ip1d = IP1D[UORDER3DQUAD];
      break;
    default:
      log(W_PROG, "Error: Component type %d not recognized.\n", utype);
      throw 1;
      exit(1);
  }
  
  u->set_mean_value = _set_mv;
  if (_set_mv)
  {
    u->mv_segment = _mv_segment;
    u->mean_value = _mv;
    u->volume     = 0;
    u->mean = new double[u->ip3d->getMeshNdof()];
    for (int j=0; j<u->ip3d->getMeshNdof(); j++) u->mean[j] = 0;
  }
  
  string bctag;
  for (map<int,string>::iterator bcname=mesh->beginPhysicalNames(); bcname!=mesh->endPhysicalNames(); bcname++)
  {
    sprintf(no, "%d", (int)U.size());
    bctag = "bc[" + bcname->second + "," + no + "]";
    if (params.find(bctag) != params.end())
    {
      p = new Parser;
      p->DefineVar("x", &x_d);
      p->DefineVar("y", &y_d);
      p->DefineVar("z", &z_d);
      p->SetExpr(params[bctag]);
      u->p[bcname->second] = p;
      params.erase(bctag);
    }
  }
  
  u->node_value = 0;
  u->node_grad  = 0;
  
  U.push_back(u);
}

void CStateProblem::addScalar(string name, int utype, bool _set_mv, int _mv_segment, double _mv)
{
  CUnknown *ui = new CUnknown;
  ui->name  = name;
  ui->type  = UT_SCALAR;
  ui->index = U.size();
  
  Uinfo.push_back(ui);
  
  addComponent(utype, _set_mv, _mv_segment, _mv);
}

void CStateProblem::addVector(string name, int utype, bool _set_mv, int _mv_segment, double _mv)
{
  CUnknown *ui = new CUnknown;
  ui->name  = name;
  ui->type  = UT_VECTOR;
  ui->index = U.size();
  
  Uinfo.push_back(ui);
  
  for (int i=0; i<3; i++) addComponent(utype, _set_mv, _mv_segment, _mv);
  
  string bctag;
  Parser *p;
  for (map<int,string>::iterator bcname=mesh->beginPhysicalNames(); bcname!=mesh->endPhysicalNames(); bcname++)
  {
    bctag = "normalbc[" + name + "," + bcname->second + "]";
    if (params.find(bctag) != params.end())
    {
      p = new Parser;
      p->DefineVar("x", &x_d);
      p->DefineVar("y", &y_d);
      p->DefineVar("z", &z_d);
      p->SetExpr(params[bctag]);
      ui->normalbc[bcname->second] = p;
      params.erase(bctag);
    }
    bctag = "tangentbc[" + name + "," + bcname->second + "]";
    if (params.find(bctag) != params.end())
    {
      p = new Parser;
      p->DefineVar("x", &x_d);
      p->DefineVar("y", &y_d);
      p->DefineVar("z", &z_d);
      p->SetExpr(params[bctag]);
      ui->tangentbc[bcname->second] = p;
      params.erase(bctag);
    }
  }
}

void CStateProblem::addTensor(string name, int utype, bool _set_mv, int _mv_segment, double _mv)
{
  CUnknown *ui = new CUnknown;
  ui->name  = name;
  ui->type  = UT_TENSOR;
  ui->index = U.size();
  
  Uinfo.push_back(ui);
  
  for (int i=0; i<9; i++) addComponent(utype, _set_mv, _mv_segment, _mv);
}

int CStateProblem::getProblemNdof()
{
  int n = 0;
  for (int i=0; i<(int)U.size(); i++)
  {
    n += U[i]->ip3d->getMeshNdof();
  }
  return n;
}

void CStateProblem::prepare()
{
  q.resize(getProblemNdof());

  if (loadSolution())
  {
    log(W_PROG, "Loaded solution from file %s.\n", params["loadfile"].c_str());
  }
  else
  {
    q = 0;
    log(W_PROG, "Starting from initial guess.\n");
  }
  q.setIndependent();
  
  drq = new t_drq[q.size()];
  r   = new double[q.size()];
  
  log(W_PROG, "Allocated data structures (problem size %d).\n", q.size());
  
  afterPrepare();
}

void CStateProblem::addResidual(CVar &r_loc, int ind)
{
  r[ind] += r_loc.getVal();
  for (int j=0; j<r_loc.getNder(); j++)
  {
    if (r_loc.getDerVar(j) >= q(0).getInd() && r_loc.getDerVar(j) < q(0).getInd()+q.size())
    {
      drq[ind][r_loc.getDerVar(j)-q(0).getInd()] += r_loc.getDer(j);
    }
  }
}

void CStateProblem::setResidual(CVar &r_loc, int ind)
{
  r[ind] = r_loc.getVal();
  drq[ind].clear();
  for (int j=0; j<r_loc.getNder(); j++)
  {
    if (r_loc.getDerVar(j) >= q(0).getInd() && r_loc.getDerVar(j) < q(0).getInd()+q.size())
    {
      drq[ind][r_loc.getDerVar(j)-q(0).getInd()] = r_loc.getDer(j);
    }
  }
}

void calculateTangentsFromNormal(CVec &n, CVec &t1, CVec &t2)
{
  CVar norm;
  t1.resize(3);
  t2.resize(3);
  if (fabs(n(2).getVal()) > 1e-3)
  {
    norm = sqrt(n(0)*n(0) + n(2)*n(2))/n(2);
    t1(0) = (CVar)1/norm;
    t1(1) = 0;
    t1(2) = -n(0)/(n(2)*norm);
  }
  else if (fabs(n(1).getVal()) > 1e-3)
  {
    norm = sqrt(n(0)*n(0) + n(1)*n(1))/n(1);
    t1(0) = (CVar)1/norm;
    t1(1) = -n(0)/(n(1)*norm);
    t1(2) = 0;
  }
  else
  {
    norm = sqrt(n(0)*n(0) + n(2)*n(2))/n(0);
    t1(0) = -n(2)/(n(0)*norm);
    t1(1) = 0;
    t1(2) = (CVar)1/norm;
  }
  t2(0) = n(1)*t1(2) - n(2)*t1(1);
  t2(1) = n(2)*t1(0) - n(0)*t1(2);
  t2(2) = n(0)*t1(1) - n(1)*t1(0);
}

void CStateProblem::assembleVolumeIntegral()
{
  int c, e, i, ind, k, offset, pct=0;
  CIntPoint *ip;
  CVar r_loc;
  t_IP::iterator it;

  log(W_SOLVE, "-Assembling volume integrals: ");

    // assemble volume integrals
    for (e=0; e<mesh->getNvolumes(); e++)
    {
      // calculate integration point data
      for (it=IP3D.begin(); it!=IP3D.end(); it++) it->second->calcE(e, EVAL_FG);
      for (k=0; k<Q3d->n; k++)
      {
        offset = 0;
        // calculate integration point data
        for (it=IP3D.begin(); it!=IP3D.end(); it++) it->second->calcQP(k, EVAL_FG);
        // evaluate unknowns at integration points
        for (c=0; c<(int)U.size(); c++)
        {
  	  U[c]->ip3d->eval(q, e, offset, U[c]->val, U[c]->grad, EVAL_FG);
	  if (tau>0)
	  {
	    U[c]->ip3d->eval(oq, e, offset, U[c]->o_val, U[c]->o_grad, EVAL_FG);
	    U[c]->ip3d->eval(ooq, e, offset, U[c]->oo_val, U[c]->oo_grad, EVAL_FG);
	  }
	  offset += U[c]->ip3d->getMeshNdof();
        }
	x = U[0]->ip3d->getX();
  	y = U[0]->ip3d->getY();
	z = U[0]->ip3d->getZ();
	afterEvaluateVolumeComponents(e, k);

        // assemble element integrals
        for (it=IP3D.begin(); it!=IP3D.end(); it++)
        {
          ip = it->second;
          for (i=0; i<ip->getNdof(); i++)
          {
            ip->evalBasis(i, phi, gphi, EVAL_FG);
            ind = ip->mapDof(e, i);
	    offset = 0;
	    for (c=0; c<(int)U.size(); c++)
	    {
	      if (U[c]->ip3d == ip && q(offset+ind).getDerOfVar(q(offset+ind).getInd()) == 1)
	      {
	        r_loc = volumeIntegral(c, e, k)*ip->dx;
		addResidual(r_loc, ind+offset);
	        if (U[c]->set_mean_value && U[c]->mv_segment == mesh->getVolumePhysicalNo(e))
		{
		  U[c]->mean[ind] += phi.getVal()*ip->dx;
		  U[c]->volume    += ip->dx/ip->getNdof();
		}
	      }
	      offset += U[c]->ip3d->getMeshNdof();
	    }
	  }
        }
      }
      if (pct < (e*100/mesh->getNvolumes()))
      {
        pct = (e*100/mesh->getNvolumes());
        logPct(W_SOLVE, pct);
      }
    } // for e
    logPct(W_SOLVE, 100);
}

void CStateProblem::assembleFaceIntegral()
{
  log(W_SOLVE, "Assembling face integrals: ");
  int c, e, i, ind, k, offset, pct=0;
  CIntPoint *ip;
  CVar r_loc;
  t_IP::iterator it;
    // assemble face integrals
    for (e=0; e<mesh->getNfaces(); e++)
    {
      // calculate integration point data
      for (it=IP2D.begin(); it!=IP2D.end(); it++) it->second->calcE(e, EVAL_FG);
      for (k=0; k<Q2d->n; k++)
      {
        offset = 0;
        // calculate integration point data
        for (it=IP2D.begin(); it!=IP2D.end(); it++) it->second->calcQP(k, EVAL_FG);
        // evaluate unknowns at integration points
        for (c=0; c<(int)U.size(); c++)
        {
      	  U[c]->ip2d->eval(q, e, offset, U[c]->val, U[c]->grad, EVAL_FG);
	  if (tau>0)
	  {
	    U[c]->ip2d->eval(oq, e, offset, U[c]->o_val, U[c]->o_grad, EVAL_FG);
	    U[c]->ip2d->eval(ooq, e, offset, U[c]->oo_val, U[c]->oo_grad, EVAL_FG);
	  }
  	  offset += U[c]->ip2d->getMeshNdof();
        }
	x = U[0]->ip2d->getX();
  	y = U[0]->ip2d->getY();
	z = U[0]->ip2d->getZ();
	afterEvaluateFaceComponents(e, k);
	
        // assemble element integrals
        for (it=IP2D.begin(); it!=IP2D.end(); it++)
        {
          ip   = it->second;
          
          for (i=0; i<ip->getNdof(); i++)
          {
            ip->evalBasis(i, phi, gphi, EVAL_FG);
            ind = ip->mapDof(e, i);
	    offset = 0;
  	    for (c=0; c<(int)U.size(); c++)
    	    {
  	      if (U[c]->ip2d == ip && q(offset+ind).getDerOfVar(q(offset+ind).getInd()) == 1)
	      {
	        r_loc = faceIntegral(c, e)*ip->dx;
		addResidual(r_loc, ind+offset);
		if (U[c]->set_mean_value && U[c]->mv_segment == mesh->getFacePhysicalNo(e))
		{
		  U[c]->mean[ind] += phi.getVal()*ip->dx;
		  U[c]->volume    += ip->dx/ip->getNdof();
		}
	      }
	      offset += U[c]->ip2d->getMeshNdof();
  	    }
	  }
        }
      }
      if (pct < (e*100/mesh->getNfaces()))
      {
        pct = (e*100/mesh->getNfaces());
        logPct(W_SOLVE, pct);
      }
    }
    logPct(W_SOLVE, 100);
}

void CStateProblem::assembleEdgeIntegral()
{
  int c, e, i, ind, k, offset, pct=0;
  CIntPoint *ip;
  CVar r_loc;
  t_IP::iterator it;
  
  // assemble edge integrals
  log(W_SOLVE, "Assembling edge integrals: ");
  for (e=0; e<mesh->getNedges(); e++)
  {
    // calculate integration point data
    for (it=IP1D.begin(); it!=IP1D.end(); it++) it->second->calcE(e, EVAL_FG);
    for (k=0; k<Q1d->n; k++)
    {
      offset = 0;
      // calculate integration point data
      for (it=IP1D.begin(); it!=IP1D.end(); it++) it->second->calcQP(k, EVAL_FG);
      // evaluate unknowns at integration points
      for (c=0; c<(int)U.size(); c++)
      {
    	U[c]->ip1d->eval(q, e, offset, U[c]->val, U[c]->grad, EVAL_FG);
	if (tau>0)
	{
	  U[c]->ip1d->eval(oq, e, offset, U[c]->o_val, U[c]->o_grad, EVAL_FG);
	  U[c]->ip1d->eval(ooq, e, offset, U[c]->oo_val, U[c]->oo_grad, EVAL_FG);
	}
	offset += U[c]->ip1d->getMeshNdof();
      }
      x = U[0]->ip1d->getX();
      y = U[0]->ip1d->getY();
      z = U[0]->ip1d->getZ();
      afterEvaluateEdgeComponents(e, k);
	
      // assemble element integrals
      for (it=IP1D.begin(); it!=IP1D.end(); it++)
      {
        ip   = it->second;
        for (i=0; i<ip->getNdof(); i++)
        {
          ip->evalBasis(i, phi, gphi, EVAL_FG);
          ind = ip->mapDof(e, i);
	  offset = 0;
	  for (c=0; c<(int)U.size(); c++)
  	  {
	    if (U[c]->ip1d == ip && q(offset+ind).getDerOfVar(q(offset+ind).getInd()) == 1)
	    {
	      r_loc = edgeIntegral(c, e)*ip->dx;
	      addResidual(r_loc, ind+offset);
	      if (U[c]->set_mean_value && U[c]->mv_segment == mesh->getEdgePhysicalNo(e))
	      {
	        U[c]->mean[ind] += phi.getVal()*ip->dx;
		U[c]->volume    += ip->dx/ip->getNdof();
              }
	    }
	    offset += U[c]->ip1d->getMeshNdof();
	  }
	}
      }
    }
    if (pct < (e*100/mesh->getNedges()))
    {
      pct = (e*100/mesh->getNedges());
      logPct(W_SOLVE, pct);
    }
  }
  logPct(W_SOLVE, 100);
}

void CStateProblem::assembleFaceNormalTangentBC()
{
  int offset, c, e, i, j, ind;
  CIntPoint *ip;
  for (map<int,string>::iterator b=mesh->beginPhysicalNames(); b != mesh->endPhysicalNames(); b++)
  {
    offset = 0;
    for (c=0; c<(int)Uinfo.size(); c++)
    {
      ip = U[Uinfo[c]->index]->ip2d;
      int dim=ip->getMeshNdof();
      if (Uinfo[c]->normalbc.find(b->second) != Uinfo[c]->normalbc.end() || Uinfo[c]->tangentbc.find(b->second) != Uinfo[c]->tangentbc.end())
      {
        map<int,int> index;
        map<int,CVec> n;
	CVec t1, t2, ln, lt1, lt2;
        index.clear();
        n.clear();
        for (e=0; e<mesh->getNfaces(); e++)
	{
	  if (mesh->getPhysicalName(mesh->getFacePhysicalNo(e)) == b->second)
	  {
	    mesh->calculateFaceNormal(e, ln, lt1, lt2);
	    for (i=0; i<ip->getNdof(); i++)
	    {
	      ind = ip->mapDof(e,i);
	      if (index.find(ind) == index.end())
	      {
	        index[ind] = mesh->getFaceNode(e, i);
		n    [ind] = ln;
	      }
	      else
	      {
		for (j=0; j<3; j++) n[ind](j) += ln(j);
	      }
	    }
	  }
	}
	CVar r1, r2, r3;
	double norm;
	t_drq d1, d2;
	int i1, i2, i3;
        for (map<int,int>::iterator it=index.begin(); it != index.end(); it++)
        {
          ind = it->first;
	  i1 = offset+ind; i2 = offset+ind+dim; i3 = offset+ind+2*dim;
	  x_d = mesh->getCoord(it->second, 0).getVal();
	  y_d = mesh->getCoord(it->second, 1).getVal();
	  z_d = mesh->getCoord(it->second, 2).getVal();
	  norm = n[ind].norm2();
	  n[ind](0) = n[ind](0)/norm; n[ind](1) = n[ind](1)/norm; n[ind](2) = n[ind](2)/norm;
	  calculateTangentsFromNormal(n[ind], t1, t2);
	  if (Uinfo[c]->normalbc.find(b->second) != Uinfo[c]->normalbc.end() && Uinfo[c]->tangentbc.find(b->second) != Uinfo[c]->tangentbc.end())
	  {
	    r1 = t1(0)*q(i1) + t1(1)*q(i2) + t1(2)*q(i3) - Uinfo[c]->tangentbc[b->second]->Eval();
	    r2 = t2(0)*q(i1) + t2(1)*q(i2) + t2(2)*q(i3) - Uinfo[c]->tangentbc[b->second]->Eval();
	    r3 = q(i1)*n[ind](0) + q(i2)*n[ind](1) + q(i3)*n[ind](2)-Uinfo[c]->normalbc[b->second]->Eval();
	    setResidual(r1, i1);
	    setResidual(r2, i2);
	    setResidual(r3, i3);
	  }
	  else if (Uinfo[c]->normalbc.find(b->second) != Uinfo[c]->normalbc.end())
	  {
	    d1.clear(); d2.clear();
	    for (t_drq::iterator it=drq[i1].begin();it!=drq[i1].end();it++)
	    {
	      d1[it->first] = t1(0).getVal()*it->second;
	      d2[it->first] = t2(0).getVal()*it->second;
	    }
	    for (t_drq::iterator it=drq[i2].begin();it!=drq[i2].end();it++)
	    {
	      d1[it->first] += t1(1).getVal()*it->second;
	      d2[it->first] += t2(1).getVal()*it->second;
	    }
	    for (t_drq::iterator it=drq[i3].begin();it!=drq[i3].end();it++)
	    {
	      d1[it->first] += t1(2).getVal()*it->second;
	      d2[it->first] += t2(2).getVal()*it->second;
	    }
	    drq[i1].clear();
	    drq[i1] = d1;
	    drq[i2].clear();
	    drq[i2] = d2;
	    r1 = t1(0)*r[i1] + t1(1)*r[i2] + t1(2)*r[i3];
	    r2 = t2(0)*r[i1] + t2(1)*r[i2] + t2(2)*r[i3];
	    r3 = q(i1)*n[ind](0) + q(i2)*n[ind](1) + q(i3)*n[ind](2)
	         - Uinfo[c]->normalbc[b->second]->Eval();
	    r[i1] = r1.getVal();
	    r[i2] = r2.getVal();
	    setResidual(r3, i3);
	  }
	  else if (Uinfo[c]->tangentbc.find(b->second) != Uinfo[c]->tangentbc.end())
	  {
	    d1.clear();
	    for (t_drq::iterator it=drq[offset+ind].begin();it!=drq[offset+ind].end();it++) d1[it->first] = n[ind](0).getVal()*it->second;
	    for (t_drq::iterator it=drq[offset+ind+dim].begin();it!=drq[offset+ind+dim].end();it++) d1[it->first] += n[ind](1).getVal()*it->second;
	    for (t_drq::iterator it=drq[offset+ind+2*dim].begin();it!=drq[offset+ind+2*dim].end();it++) d1[it->first] += n[ind](2).getVal()*it->second;
	    drq[ind+offset+2*dim].clear();
	    drq[ind+offset+2*dim] = d1;
	    r1 = t1(0)*q(offset+ind) + t1(1)*q(offset+ind+dim) + t1(2)*q(offset+ind+2*dim) - Uinfo[c]->tangentbc[b->second]->Eval();
	    r2 = t2(0)*q(offset+ind) + t2(1)*q(offset+ind+dim) + t2(2)*q(offset+ind+2*dim) - Uinfo[c]->tangentbc[b->second]->Eval();
	    r3 = n[ind](0)*r[offset+ind] + n[ind](1)*r[offset+dim+ind] + n[ind](2)*r[offset+2*dim+ind];
	    r[ind+offset+2*dim] = r3.getVal();
	    setResidual(r1, ind+offset);
	    setResidual(r2, ind+offset+dim);
	  }
        }
      }
      offset += Uinfo[c]->type*dim;
    }
  }
}

void CStateProblem::assembleEdgeNormalTangentBC()
{
  int offset, c, e, i, j, ind;
  CIntPoint *ip;
  for (map<int,string>::iterator b=mesh->beginPhysicalNames(); b != mesh->endPhysicalNames(); b++)
  {
    offset = 0;
    for (c=0; c<(int)Uinfo.size(); c++)
    {
      ip = U[Uinfo[c]->index]->ip1d;
      int dim=ip->getMeshNdof();
      if (Uinfo[c]->normalbc.find(b->second) != Uinfo[c]->normalbc.end() || Uinfo[c]->tangentbc.find(b->second) != Uinfo[c]->tangentbc.end())
      {
        map<int,int> index;
        map<int,CVec> t;
	CVec n1, n2, lt;
        index.clear();
        t.clear();
        for (e=0; e<mesh->getNedges(); e++)
	{
	  if (mesh->getPhysicalName(mesh->getEdgePhysicalNo(e)) == b->second)
	  {
	    mesh->calculateEdgeTangent(e, lt);
	    for (i=0; i<ip->getNdof(); i++)
	    {
	      ind = ip->mapDof(e,i);
	      if (index.find(ind) == index.end())
	      {
	        index[ind] = mesh->getEdgeNode(e, i);
		t    [ind] = lt;
	      }
	      else
	      {
		for (j=0; j<3; j++) t[ind](j) += lt(j);
	      }
	    }
	  }
	}
	CVar r1, r2, r3;
	double norm;
	t_drq d1, d2;
        for (map<int,int>::iterator it=index.begin(); it != index.end(); it++)
        {
          ind = it->first;
	  x_d = mesh->getCoord(it->second, 0).getVal();
	  y_d = mesh->getCoord(it->second, 1).getVal();
	  z_d = mesh->getCoord(it->second, 2).getVal();
	  norm = t[ind].norm2();
	  t[ind](0) = t[ind](0)/norm; t[ind](1) = t[ind](1)/norm; t[ind](2) = t[ind](2)/norm;
	  calculateTangentsFromNormal(t[ind], n1, n2);
	  if (Uinfo[c]->normalbc.find(b->second) != Uinfo[c]->normalbc.end() && Uinfo[c]->tangentbc.find(b->second) != Uinfo[c]->tangentbc.end())
	  {
	    r1 = n1(0)*q(offset+ind) + n1(1)*q(offset+ind+dim) + n1(2)*q(offset+ind+2*dim)-Uinfo[c]->normalbc[b->second]->Eval();
	    r2 = n2(0)*q(offset+ind) + n2(1)*q(offset+ind+dim) + n2(2)*q(offset+ind+2*dim)-Uinfo[c]->normalbc[b->second]->Eval();
	    r3 = q(offset+ind)*t[ind](0) + q(offset+dim+ind)*t[ind](1) + q(offset+2*dim+ind)*t[ind](2)-Uinfo[c]->tangentbc[b->second]->Eval();
	    setResidual(r1, ind+offset);
	    setResidual(r2, ind+offset+dim);
	    setResidual(r3, ind+offset+2*dim);
	  }
	  else if (Uinfo[c]->tangentbc.find(b->second) != Uinfo[c]->tangentbc.end())
	  {
	    d1.clear(); d2.clear();
	    for (t_drq::iterator it=drq[offset+ind].begin();it!=drq[offset+ind].end();it++)
	    {
	      d1[it->first] = n1(0).getVal()*it->second;
	      d2[it->first] = n2(0).getVal()*it->second;
	    }
	    for (t_drq::iterator it=drq[offset+ind+dim].begin();it!=drq[offset+ind+dim].end();it++)
	    {
	      d1[it->first] += n1(1).getVal()*it->second;
	      d2[it->first] += n2(1).getVal()*it->second;
	    }
	    for (t_drq::iterator it=drq[offset+ind+2*dim].begin();it!=drq[offset+ind+2*dim].end();it++)
	    {
	      d1[it->first] += n1(2).getVal()*it->second;
	      d2[it->first] += n2(2).getVal()*it->second;
	    }
	    drq[ind+offset].clear();
	    drq[ind+offset] = d1;
	    drq[ind+offset+dim].clear();
	    drq[ind+offset+dim] = d2;
	    r1 = n1(0)*r[offset+ind] + n1(1)*r[offset+ind+dim] + n1(2)*r[offset+ind+2*dim];
	    r2 = n2(0)*r[offset+ind] + n2(1)*r[offset+ind+dim] + n2(2)*r[offset+ind+2*dim];
	    r3 = q(offset+ind)*t[ind](0) + q(offset+dim+ind)*t[ind](1) + q(offset+2*dim+ind)*t[ind](2)
	         - Uinfo[c]->tangentbc[b->second]->Eval();
	    r[ind+offset] = r1.getVal();
	    r[ind+offset+dim] = r2.getVal();
	    setResidual(r3, ind+offset+2*dim);
	  }
	  else if (Uinfo[c]->normalbc.find(b->second) != Uinfo[c]->normalbc.end())
	  {
	    d1.clear();
	    for (t_drq::iterator it=drq[offset+ind].begin();it!=drq[offset+ind].end();it++) d1[it->first] = t[ind](0).getVal()*it->second;
	    for (t_drq::iterator it=drq[offset+ind+dim].begin();it!=drq[offset+ind+dim].end();it++) d1[it->first] += t[ind](1).getVal()*it->second;
	    for (t_drq::iterator it=drq[offset+ind+2*dim].begin();it!=drq[offset+ind+2*dim].end();it++) d1[it->first] += t[ind](2).getVal()*it->second;
	    drq[ind+offset+2*dim].clear();
	    drq[ind+offset+2*dim] = d1;
	    r1 = n1(0)*q(offset+ind) + n1(1)*q(offset+ind+dim) + n1(2)*q(offset+ind+2*dim) - Uinfo[c]->normalbc[b->second]->Eval();
	    r2 = n2(0)*q(offset+ind) + n2(1)*q(offset+ind+dim) + n2(2)*q(offset+ind+2*dim) - Uinfo[c]->normalbc[b->second]->Eval();
	    r3 = t[ind](0)*r[offset+ind] + t[ind](1)*r[offset+dim+ind] + t[ind](2)*r[offset+2*dim+ind];
	    r[ind+offset+2*dim] = r3.getVal();
	    setResidual(r1, ind+offset);
	    setResidual(r2, ind+offset+dim);
	  }
        }
      }
      offset += Uinfo[c]->type*dim;
    }
  }
}

void CStateProblem::setDirichletBC()
{
  int c, e, i, ind, offset;
  log(W_SOLVE, "-Setting boundary conditions.\n");
  for (e=0; e<mesh->getNfaces(); e++)
  {
    for (i=0; i<mesh->getNfaceDofs(); i++)
    {
      ind = mesh->getFaceNode(e, i);
      x = mesh->getCoord(ind, 0);
      y = mesh->getCoord(ind, 1);
      z = mesh->getCoord(ind, 2);
      x_d = x.getVal();
      y_d = y.getVal();
      z_d = z.getVal();
      
      offset = 0;
      for (c=0; c<(int)U.size(); c++)
      {
        if (U[c]->p.find(mesh->getPhysicalName(mesh->getFacePhysicalNo(e))) != U[c]->p.end())
        {
          q(offset+ind) = U[c]->p[mesh->getPhysicalName(mesh->getFacePhysicalNo(e))]->Eval();
	}
	offset += U[c]->ip3d->getMeshNdof();
      }
    }
  }
  for (e=0; e<mesh->getNedges(); e++)
  {
    for (i=0; i<mesh->getNedgeDofs(); i++)
    {
      ind = mesh->getEdgeNode(e, i);
      x = mesh->getCoord(ind, 0);
      y = mesh->getCoord(ind, 1);
      z = mesh->getCoord(ind, 2);
      x_d = x.getVal();
      y_d = y.getVal();
      z_d = z.getVal();
      
      offset = 0;
      for (c=0; c<(int)U.size(); c++)
      {
        if (U[c]->p.find(mesh->getPhysicalName(mesh->getEdgePhysicalNo(e))) != U[c]->p.end())
        {
          q(offset+ind) = U[c]->p[mesh->getPhysicalName(mesh->getEdgePhysicalNo(e))]->Eval();
	}
	offset += U[c]->ip3d->getMeshNdof();
      }
    }
  }
}

void CStateProblem::setResidualsForDirichletBC()
{
  int c, e, i, ind, offset;
  for (e=0; e<mesh->getNfaces(); e++)
  {
    for (i=0; i<mesh->getNfaceDofs(); i++)
    {
      ind = mesh->getFaceNode(e, i);
      offset = 0;
      for (c=0; c<(int)U.size(); c++)
      {
        if (U[c]->p.find(mesh->getPhysicalName(mesh->getFacePhysicalNo(e))) != U[c]->p.end())
        {
	  r[offset+ind] = 0;
	  drq[offset+ind].clear();
	  drq[offset+ind][offset+ind] = 1;
	}
	offset += U[c]->ip3d->getMeshNdof();
      }
    }
  }
  for (e=0; e<mesh->getNedges(); e++)
  {
    for (i=0; i<mesh->getNedgeDofs(); i++)
    {
      ind = mesh->getEdgeNode(e, i);
      offset = 0;
      for (c=0; c<(int)U.size(); c++)
      {
        if (U[c]->p.find(mesh->getPhysicalName(mesh->getEdgePhysicalNo(e))) != U[c]->p.end())
        {
	  r[offset+ind] = 0;
	  drq[offset+ind].clear();
	  drq[offset+ind][offset+ind] = 1;
	}
	offset += U[c]->ip3d->getMeshNdof();
      }
    }
  }
}

void CStateProblem::setMeanValues()
{
  int c, i, e, offset = 0, ind0 = -1, ind;
  for (c=0; c<(int)U.size(); c++)
  {
    if (U[c]->set_mean_value)
    {
      for (e=0; e<mesh->getNvolumes(); e++)
      {
        if (U[c]->mv_segment == mesh->getVolumePhysicalNo(e))
	{
	  for (i=0; i<U[c]->ip3d->getNdof(); i++)
	  {
	    ind = U[c]->ip3d->mapDof(e, i);
	    if (ind0 == -1)
	    {
	      ind0 = ind;
              r[offset+ind0] = -U[c]->mean_value*U[c]->volume;
              drq[offset+ind0].clear();
	    }
	    if (U[c]->mean[ind] != 0)
	    {
              r[offset+ind0] += q(offset+ind).getVal()*U[c]->mean[ind];
              drq[offset+ind0][offset+ind] = U[c]->mean[ind];
	      U[c]->mean[ind] = 0;
	    }
	  }
	}
      }
      for (e=0; e<mesh->getNfaces(); e++)
      {
        if (U[c]->mv_segment == mesh->getFacePhysicalNo(e))
	{
	  for (i=0; i<U[c]->ip2d->getNdof(); i++)
	  {
	    ind = U[c]->ip2d->mapDof(e, i);
	    if (ind0 == -1)
	    {
	      ind0 = ind;
              r[offset+ind0] = -U[c]->mean_value*U[c]->volume;
              drq[offset+ind0].clear();
	    }
	    if (U[c]->mean[ind] != 0)
	    {
              r[offset+ind0] += q(offset+ind).getVal()*U[c]->mean[ind];
              drq[offset+ind0][offset+ind] = U[c]->mean[ind];
	      U[c]->mean[ind] = 0;
	    }
	  }
	}
      }
      for (e=0; e<mesh->getNedges(); e++)
      {
        if (U[c]->mv_segment == mesh->getEdgePhysicalNo(e))
	{
	  for (i=0; i<U[c]->ip1d->getNdof(); i++)
	  {
	    ind = U[c]->ip1d->mapDof(e, i);
	    if (ind0 == -1)
	    {
	      ind0 = ind;
              r[offset+ind0] = -U[c]->mean_value*U[c]->volume;
              drq[offset+ind0].clear();
	    }
	    if (U[c]->mean[ind] != 0)
	    {
              r[offset+ind0] += q(offset+ind).getVal()*U[c]->mean[ind];
              drq[offset+ind0][offset+ind] = U[c]->mean[ind];
	      U[c]->mean[ind] = 0;
	    }
	  }
	}
      }
    }
    offset += U[c]->ip3d->getMeshNdof();
  }
}

void CStateProblem::assemble()
{
  clock_t t;
  
  log(W_SOLVE, "Starting assembly:\n");
  start_clock(t);
  
  for (int i=0; i<q.size(); i++)
  {
    r[i] = 0;
    drq[i].clear();
  }
  for (int c=0; c<(int)U.size(); c++)
  {
    if (!U[c]->set_mean_value || U[c]->mean == 0) continue;
    for (int i=0; i<U[c]->ip3d->getMeshNdof(); i++)
    {
      U[c]->mean[i] = 0;
    }
    U[c]->volume = 0;
  }
  
  // set boundary conditions
  setDirichletBC();

  // assemble integrals
  if (useVolumeIntegral()) assembleVolumeIntegral();
  if (useFaceIntegral()) assembleFaceIntegral();
  if (useEdgeIntegral()) assembleEdgeIntegral();
  
  // set normal and tangent boundary conditions
  assembleFaceNormalTangentBC();
  assembleEdgeNormalTangentBC();

  // set boundary conditions
  setResidualsForDirichletBC();

  // set mean values
  setMeanValues();
  
  log(W_SOLVE, "Assembly completed in %f s.\n", stop_clock(t));
}

int CStateProblem::solve(double relax)
{
  CSparseMat sp;
  t_drq::iterator it;
  int i, j;
  int nnz = 0;
  int *ia, *ja;
  double *v;
  clock_t t;
  
  log(W_SOLVE, "Allocating sparse matrix.\n");
  
  for (i=0; i<q.size(); i++)
  {
    nnz += drq[i].size();
  }
  sp.resize(q.size(), q.size(), nnz);
  ia = sp.getIa();
  ja = sp.getJa();
  v  = sp.getV();
  for (i=0; i<q.size(); i++)
  {
    ia[i+1] = ia[i] + drq[i].size();
    for (j=0, it=drq[i].begin(); it !=drq[i].end(); j++, it++)
    {
        ja[ia[i]+j] = it->first;
        v[ia[i]+j]  = it->second;
    }
    drq[i].clear();
  }
  
  int n = sp.getNrows();
  
  log(W_SOLVE, "Starting solution.\n");
  start_clock(t);
  
  if (UMFPACKSolve(sp, r, false) != 0)
  {
    log(W_SOLVE, "Solution failed.\n");
    return 1;
  }

  for (i=0; i<n; i++)
    q(i) = q(i)-r[i]*relax;

  log(W_SOLVE, "Solution completed in %f s.\n\n", stop_clock(t));
    
  return 0;
}

void CStateProblem::nlSolve(double relax)
{
  double rnorm = 0;
  int it = 0;
  clock_t t;
  double t_asm, t_sol;

  log(W_PROG, "\nStarting nonlinear loop.\n");
  log(W_PROG, "+----+---------+--------+--------+\n");
  log(W_PROG, "|iter|residuum |time_asm|time_sol|\n");
  log(W_PROG, "+----+---------+--------+--------+\n");
  
  
  start_clock(t);
  assemble();
  t_asm = stop_clock(t);
  
  for (int i=0; i<q.size(); rnorm += r[i]*r[i], i++);
  rnorm = sqrt(rnorm);
  
  while (rnorm > 1e-6)
  {
    start_clock(t);
    solve(relax);
    t_sol = stop_clock(t);
    
    log(W_PROG, "|%4d|%-9.3e|%8.2f|%8.2f|\n", it, rnorm, t_asm, t_sol);
    
    it++;
    log(W_SOLVE, "Nonlinear iteration no. %d:\n", it);
    
    start_clock(t);
    assemble();
    t_asm = stop_clock(t);
    
    rnorm = 0;
    for (int i=0; i<q.size(); rnorm += r[i]*r[i], i++);
    rnorm = sqrt(rnorm);
  }
  log(W_PROG, "|%4d|%-9.3e|%8.2f|  ----  |\n", it, rnorm, t_asm);
  log(W_PROG, "+----+---------+--------+--------+\n\n");
}

void CStateProblem::tSolve(double relax)
{
  int it = 0;
  clock_t t;
  double t_asm, t_sol;
  char fn[255], pvd_fn[255];

  log(W_PROG, "\nStarting time loop.\n");
  log(W_PROG, "+--------+--------+--------+\n");
  log(W_PROG, "|time_sim|time_asm|time_sol|\n");
  log(W_PROG, "+--------+--------+--------+\n");
  
  oq.resize(q.size());
  ooq.resize(q.size());
  sprintf(fn, "%s.%d.vtu", params["savebase"].c_str(), it);
  saveVTKview(fn);
  
  int pos = params["savebase"].find_last_of('/');
  if (pos > (int)params["savebase"].length()) pos = 0;
  sprintf(pvd_fn, "%s.pvd", params["savebase"].substr(pos+1).c_str());
  ofstream f(pvd_fn);
  f << "<?xml version=\"1.0\"?>\n"
    << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n"
    << "<Collection>\n"
    << "<DataSet timestep=\"0\" group=\"\" part=\"0\" file=\"" << params["savebase"].c_str() << ".0.vtu\"/>\n";
  
  while (tau*it < tmax)
  {
    it++;
    log(W_SOLVE, "Solving time %.2f:\n", tau*it);

    for (int i=0; i<q.size(); ooq(i) = oq(i).getVal(), oq(i) = q(i).getVal(), i++);
    calculated_solution = false;
    
    start_clock(t);
    assemble();
    t_asm = stop_clock(t);
    
    start_clock(t);
    solve(relax);
    t_sol = stop_clock(t);
    
    log(W_PROG, "|%8.2f|%8.2f|%8.2f|\n", tau*it, t_asm, t_sol);
    
    sprintf(fn, "%s.%d.vtu", params["savebase"].c_str(), it);
    saveVTKview(fn);
    sprintf(fn, "%s.%d.sol", params["savebase"].c_str(), it);
    saveSolution(fn);
    f << "<DataSet timestep=\"" << it << "\" group=\"\" part=\"0\" file=\"" << params["savebase"].c_str() << "." << it << ".vtu\"/>\n";
    f.flush();
  }
  
  log(W_PROG, "+--------+--------+--------+\n\n");
  
  f << "</Collection>\n"
    << "</VTKFile>";
  f.close();
  log(W_SOLVE, "Animation saved to %s.", fn);
}

void CStateProblem::saveMatrix(const char *file)
{
  int n = q.size();
  ofstream f(file);

  for (int i=0; i<n; i++)
  {
    for (t_drq::iterator it=drq[i].begin(); it!=drq[i].end(); it++)
    {
      if (it->second != 0)
      {
        f << i << " " << it->first  << " " << it->second << endl;
      }
    }
  }
  
  f.close();
}

void CStateProblem::addPPscalar(string name, string expression)
{
  Parser *p = new Parser;
  p->SetExpr(expression);
  
  PPscalars[name] = p;
}

void CStateProblem::calculateSolution()
{
  int e, k, c, offset, ind;
  t_IP::iterator it;
  CIntPoint *ip;
  CVar val;
  CVec grad;
  
  if (calculated_solution) return;
  
  for (c=0; c<(int)U.size(); c++)
  {
    if (U[c]->node_value == 0) U[c]->node_value = new double[  mesh->getNnodes()];
    if (U[c]->node_grad  == 0) U[c]->node_grad  = new double[3*mesh->getNnodes()];
  }
  
  for (e=0; e<mesh->getNvolumes(); e++)
  {
    for (t_IP::iterator it=IP3D.begin(); it != IP3D.end(); it++)
    {
      ip = it->second;
      ip->calcE(e, EVAL_FG);
      for (k=0; k<mesh->getNvolumeDofs(); k++)
      {
        ip->calcQP(k+ip->getQ()->n, EVAL_FG);
	ind = mesh->getVolumeNode(e, k);
	offset = 0;
	for (c=0; c<(int)U.size(); c++)
	{
	  if (U[c]->ip3d == ip)
	  {
	    ip->eval(q, e, offset, val, grad, EVAL_FG);
	    U[c]->node_value[ind] = val.getVal();
	    U[c]->node_grad[                    ind] = grad(0).getVal();
	    U[c]->node_grad[mesh->getNnodes()  +ind] = grad(1).getVal();
	    U[c]->node_grad[mesh->getNnodes()*2+ind] = grad(2).getVal();
	  }
	  offset += U[c]->ip3d->getMeshNdof();
	}
      }
    }
  }
  
  for (e=0; e<mesh->getNfaces(); e++)
  {
    for (t_IP::iterator it=IP2D.begin(); it != IP2D.end(); it++)
    {
      ip = it->second;
      ip->calcE(e, EVAL_FG);
      for (k=0; k<mesh->getNfaceDofs(); k++)
      {
        ip->calcQP(k+ip->getQ()->n, EVAL_FG);
	ind = mesh->getFaceNode(e, k);
	offset = 0;
	for (c=0; c<(int)U.size(); c++)
	{
	  if (U[c]->ip2d == ip)
	  {
	    ip->eval(q, e, offset, val, grad, EVAL_FG);
	    U[c]->node_value[ind] = val.getVal();
	    U[c]->node_grad[                    ind] = grad(0).getVal();
	    U[c]->node_grad[mesh->getNnodes()  +ind] = grad(1).getVal();
	    U[c]->node_grad[mesh->getNnodes()*2+ind] = grad(2).getVal();
	  }
	  offset += U[c]->ip2d->getMeshNdof();
	}
      }
    }
  }
  
  calculated_solution = true;
}


void CStateProblem::saveVTKview(const char *file)
{
  char fn[255];
  if (file == 0)
  {
    if (params.find("savebase") != params.end())
    {
      sprintf(fn, "%s.vtu", params["savebase"].c_str());
      file = fn;
    }
    else
    {
      if (tau>0)
      {
        log(W_SOLVE, "VTK XML view will not be saved.\n");
      }
      else
      {
        log(W_PROG, "VTK XML view will not be saved.\n");
      }
      return;
    }
  }
  
  calculateSolution();
  
  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
  ug->Allocate();
  
  vtkPoints *points = vtkPoints::New();
  vtkDoubleArray *da;
  int pts[10];
  int i, j, ind=0;
  
  // insert points
  for (i=0; i<mesh->getNnodes(); i++)
  {
    points->InsertNextPoint(mesh->getCoord(i, 0).getVal(), mesh->getCoord(i, 1).getVal(), mesh->getCoord(i, 2).getVal());
  }
  ug->SetPoints(points);
  points->Delete();
  // insert volumes
  for (i=0; i<mesh->getNvolumes(); i++)
  {
    for (j=0; j<mesh->getNvolumeDofs(); j++)
    {
      pts[j] = mesh->getVtkVolumeNode(i, j);
    }
    ug->InsertNextCell(mesh->getVtk3Dtype(), mesh->getNvolumeDofs(), pts);
  }
  // insert faces
  for (i=0; i<mesh->getNfaces(); i++)
  {
    for (j=0; j<mesh->getNfaceDofs(); j++)
    {
      pts[j] = mesh->getFaceNode(i, j);
    }
    ug->InsertNextCell(mesh->getVtk2Dtype(), mesh->getNfaceDofs(), pts);
  }
  // insert edges
  for (i=0; i<mesh->getNedges(); i++)
  {
    for (j=0; j<mesh->getNedgeDofs(); j++)
    {
      pts[j] = mesh->getEdgeNode(i, j);
    }
    ug->InsertNextCell(mesh->getVtk1Dtype(), mesh->getNedgeDofs(), pts);
  }
  // insert unknowns data
  for (vector<CUnknown *>::iterator ui = Uinfo.begin(); ui != Uinfo.end(); ui++)
  {
    da = vtkDoubleArray::New();
    da->SetName((*ui)->name.c_str());
    switch ((*ui)->type)
    {
      case UT_SCALAR:
        da->SetNumberOfComponents(1);
        for (i=0; i<mesh->getNnodes(); i++)
        {
          da->InsertTuple1(i, U[(*ui)->index]->node_value[i]);
        }
	ind = ug->GetPointData()->AddArray(da);
	ug->GetPointData()->SetActiveAttribute(ind, vtkDataSetAttributes::SCALARS);
        break;
	  
      case UT_VECTOR:
        da->SetNumberOfComponents(3);
        for (i=0; i<mesh->getNnodes(); i++)
        {
          da->InsertTuple3(i, U[(*ui)->index]->node_value[i], U[(*ui)->index+1]->node_value[i], U[(*ui)->index+2]->node_value[i]);
        }
	ind = ug->GetPointData()->AddArray(da);
	ug->GetPointData()->SetActiveAttribute(ind, vtkDataSetAttributes::VECTORS);
        break;
	
      case UT_TENSOR:
        da->SetNumberOfComponents(9);
        for (i=0; i<mesh->getNnodes(); i++)
        {
          da->InsertTuple9(i, U[(*ui)->index  ]->node_value[i], U[(*ui)->index+1]->node_value[i], U[(*ui)->index+2]->node_value[i],
                              U[(*ui)->index+3]->node_value[i], U[(*ui)->index+4]->node_value[i], U[(*ui)->index+5]->node_value[i],
                              U[(*ui)->index+6]->node_value[i], U[(*ui)->index+7]->node_value[i], U[(*ui)->index+8]->node_value[i]);
        }
	ind = ug->GetPointData()->AddArray(da);
	ug->GetPointData()->SetActiveAttribute(ind, vtkDataSetAttributes::TENSORS);
        break;
    } // end of unknowns switch
    da->Delete();
  } // end of unknowns loop
  // save additional scalar fields
  for (map<string,Parser *>::iterator it = PPscalars.begin(); it != PPscalars.end(); it++)
  {
    da = vtkDoubleArray::New();
    da->SetName(it->first.c_str());
    da->SetNumberOfComponents(1);

    Parser *p = it->second;
    char v[255];
    double node;
    
    for (i=0; i<mesh->getNnodes(); i++)
    {
      x_d = mesh->getCoord(i, 0).getVal();
      y_d = mesh->getCoord(i, 1).getVal();
      z_d = mesh->getCoord(i, 2).getVal();
      node = i;
      p->ClearVar();
      p->DefineVar("x", &x_d);
      p->DefineVar("y", &y_d);
      p->DefineVar("z", &z_d);
      sprintf(v, "node");
      p->DefineVar(v, &node);
      for (vector<CUnknown *>::iterator ui = Uinfo.begin(); ui != Uinfo.end(); ui++)
      {
        switch ((*ui)->type)
        {
          case UT_SCALAR:
  	    p->DefineVar((*ui)->name, &U[(*ui)->index]->node_value[i]);
	    for (j=0; j<3; j++)
             {
	      sprintf(v, "g%s%d", (*ui)->name.c_str(), j);
	      p->DefineVar(v, &U[(*ui)->index]->node_grad[mesh->getNnodes()*j+i]);
	    }
	    break;
	  case UT_VECTOR:
	    for (int c=0; c<3; c++)
	    {
  	      sprintf(v, "%s%d", (*ui)->name.c_str(), c);
	      p->DefineVar(v, &U[(*ui)->index+c]->node_value[i]);
	      for (int j=0; j<3; j++)
	      {
	        sprintf(v, "g%s%d%d", (*ui)->name.c_str(), c, j);
	        p->DefineVar(v, &U[(*ui)->index+c]->node_grad[mesh->getNnodes()*j+i]);
	      }
	    }
	    break;
	  case UT_TENSOR:
	    for (int tr=0; tr<3; tr++)
	    {
	      for (int tc=0; tc<3; tc++)
	      {
	        sprintf(v, "%s%d%d", (*ui)->name.c_str(), tr, tc);
	        p->DefineVar(v, &U[(*ui)->index+tr*3+tc]->node_value[i]);
	        for (int j=0; j<3; j++)
	        {
	          sprintf(v, "g%s%d%d%d", (*ui)->name.c_str(), tr, tc, j);
	          p->DefineVar(v, &U[(*ui)->index+tr*3+tc]->node_grad[mesh->getNnodes()*j+i]);
	        }
	      }
	    }
	    break;
        }
      }
      
      da->InsertTuple1(i, p->Eval());
    }
    ind = ug->GetPointData()->AddArray(da);
    ug->GetPointData()->SetActiveAttribute(ind, vtkDataSetAttributes::SCALARS);
    da->Delete();
  }
  // save data to VTK file
  vtkXMLUnstructuredGridWriter *ugw = vtkXMLUnstructuredGridWriter::New();
  ugw->SetInput(ug);
  ugw->SetFileName(file);
  // do not base64() encode the appended data
//  ugw->SetDataModeToAscii();
  ugw->SetEncodeAppendedData(false);
  ugw->Write();
  ugw->Delete();
  ug->Delete();
  
  if (tau>0)
  {
    log(W_SOLVE, "Results saved to VTK file %s.\n", file);
  }
  else
  {
    log(W_PROG, "Results saved to VTK file %s.\n", file);
  }
}

bool CStateProblem::loadSolution(const char *fname)
{
  if (fname == 0)
  {
    if (params.find("loadfile") != params.end())
    {
      fname = params["loadfile"].c_str();
    }
    else
    {
      loaded_solution = false;
      return false;
    }
  }
  
  ifstream f(fname);
  if ((bool)f == false)
  {
    loaded_solution = false;
    return false;
  }
  
  int size;
  double v;
  
  f >> size;
  if (size != q.size())
  {
    log(W_PROG, "Solution file %s is incompatible with problem size.\n", fname);
    loaded_solution = false;
    return false;
  }
  
  for (int i=0; i<size; f >> v, q(i) = v, i++);
  
  f.close();
  
  loaded_solution = true;
  return true;
}

void CStateProblem::saveSolution(const char *fname)
{
  FILE *f;
  char fn[255];
  
  if (fname == 0)
  {
    if (params.find("savebase") == params.end())
    {
      if (tau>0)
      {
        log(W_SOLVE, "Solution will not be saved.\n");
      }
      else
      {
        log(W_PROG, "Solution will not be saved.\n");
      }
      return;
    }
    else
    {
      sprintf(fn, "%s.sol", params["savebase"].c_str());
      fname = fn;
    }
  }
  
  f = fopen(fname, "w");
  
  fprintf(f, "%d\n", q.size());
  for (int i=0; i<q.size(); i++)
  {
    fprintf(f, "%.15f ", q(i).getVal());
  }
  
  fclose(f);
  
  if (tau>0)
  {
    log(W_SOLVE, "Solution saved to %s.\n", fname);
  }
  else
  {
    log(W_PROG, "Solution saved to %s.\n", fname);
  }
}

void CStateProblem::log(int window, const char *format, ...)
{
  if (wm == 0) return;
  
  WINDOW *w = 0;
  va_list l;
  
  switch (window)
  {
    case W_PROG:
      w = wm->getProgW();
      break;
    case W_SOLVE:
      w = wm->getSolveW();
      break;
  }
  
  va_start(l, format);
  vwprintw(w, format, l);
  va_end(l);
  
  // print the same text to log file
  if (logfile != 0 && window == W_PROG)
  {
    va_start(l, format);
    vfprintf(logfile, format, l);
    va_end(l);
    fflush(logfile);
  }
  
  wrefresh(w);
}

void CStateProblem::logPct(int window, int pct)
{
  if (wm == 0) return;
  
  WINDOW *w = 0;
  int posx, posy;
  
  switch (window)
  {
    case W_PROG:
      w = wm->getProgW();
      break;
    case W_SOLVE:
      w = wm->getSolveW();
      break;
  }
  
  getyx(w, posy, posx);
  wprintw(w, "%d%%", pct);
  
  if (pct >= 100)
  {
    wprintw(w, "\n");
    // print the same text to log file
    if (logfile != 0 && window == W_PROG)
    {
      fprintf(logfile, "\n");
    }
  }
  else
  {
    wmove(w, posy, posx);
  }
  wrefresh(w);
}


