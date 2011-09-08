/* Finite element subroutines: mesh data structures
   Jukka Toivanen, Jan Stebel
*/

#include <math.h>
#include <map>
#include "mesh3d.h" 


void C3DMesh::destroy()
{
  coords.clear();
  physicalNames.clear();
}

C3DMesh::~C3DMesh()
{
  destroy();
}

void C3DMesh::deform(CVar& func(double, double, double, int))
{
  int i;
  double x, y, z;
  
  // loop through all nodes
  for (i=0; i<nnodes; i++)
  {
    x = coords[3*i  ].getVal();
    y = coords[3*i+1].getVal();
    z = coords[3*i+2].getVal();
    
    coords[3*i  ] = func(x, y, z, 0);
    coords[3*i+1] = func(x, y, z, 1);
    coords[3*i+2] = func(x, y, z, 2);
  }
}

CVar& C3DMesh::getCoord(int node, int coord)
{
  return coords[3*node+coord];
}

void C3DMesh::setCoord(int node, int coord, CVar x)
{
  coords[3*node+coord] = x;
}

void C3DMesh::calculateFaceNormal(int fid, CVec &n, CVec &a, CVec &b)
{
  int i;
  double anorm, bnorm;
  a.resize(3);
  b.resize(3);
  n.resize(3);
  for (i=0; i<3; i++)
  {
    a(i) = getCoord(getFaceNode(fid, 1), i) - getCoord(getFaceNode(fid, 0), i);
    b(i) = getCoord(getFaceNode(fid, 2), i) - getCoord(getFaceNode(fid, 0), i);
  }
  anorm = a.norm2();
  bnorm = b.norm2();
  for (i=0; i<3; i++)
  {
    a(i) = a(i)/anorm;
    b(i) = b(i)/bnorm;
  }
  n(0) = a(1)*b(2) - a(2)*b(1);
  n(1) = a(2)*b(0) - a(0)*b(2);
  n(2) = a(0)*b(1) - a(1)*b(0);
}

void C3DMesh::calculateEdgeTangent(int eid, CVec &t)
{
  int i;
  double norm = 0;
  t.resize(3);
  for (i=0; i<3; i++)
  {
    t(i) = getCoord(getEdgeNode(eid, 1), i) - getCoord(getEdgeNode(eid, 0), i);
    norm += t(i).getVal()*t(i).getVal();
  }
  norm = sqrt(norm);
  for (i=0; i<3; i++) t(i) = t(i)/norm;
}

/*CVar C3DMesh::getElemSize(int ind) {
  CVar x1, x2, x3, y1, y2, y3;
  int pt;
  
  pt = getElem(ind, 0); x1 = getCoord(pt, 0); y1 = getCoord(pt, 1);
  pt = getElem(ind, 1); x2 = getCoord(pt, 0); y2 = getCoord(pt, 1);
  pt = getElem(ind, 2); x3 = getCoord(pt, 0); y3 = getCoord(pt, 1);
  
  return (CVar)0.5*abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
}

CVar C3DMesh::getElemDiam(int ind) {
  CVar x1, x2, x3, y1, y2, y3;
  int pt;
  
  pt = getElem(ind, 0); x1 = getCoord(pt, 0); y1 = getCoord(pt, 1);
  pt = getElem(ind, 1); x2 = getCoord(pt, 0); y2 = getCoord(pt, 1);
  pt = getElem(ind, 2); x3 = getCoord(pt, 0); y3 = getCoord(pt, 1);
  
  CVar d12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  CVar d23 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
  CVar d13 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
  
  if (d12.getVal() >= d23.getVal() && d12.getVal() >= d23.getVal()) return d12;
  if (d23.getVal() >= d13.getVal()) return d23;
  return d13;
}*/


C3DLinMesh::~C3DLinMesh()
{
  volumes.clear();
  faces.clear();
  edges.clear();
}

void C3DLinMesh::readGMSH(string fname)
{
  int i, j, ip, nelems, n, type, ntags, tag;
  int nphysical, ngeom, npart;
  double x, y, z;
  string s, s1;
  ifstream f(fname.c_str());
  map<int, int> mapNodes;
  CTetrahedron tet;
  CTriangle tri;
  CLine lin;
  
  try
  {
    getline(f, s);
    if (s.compare("$MeshFormat") != 0) throw 1;
    getline(f, s);
    sscanf(s.c_str(), "2 0 %d", &i);
    if (i != sizeof(double)) throw 2;
    getline(f, s);
    if (s.compare("$EndMeshFormat") != 0) throw 1;
  
    getline(f, s);
    if (s.compare("$PhysicalNames") == 0)
    { // read through the physical names block
      getline(f, s);
      sscanf(s.c_str(), "%d", &n);
      for (i=0; i<n; i++)
      {
        getline(f, s);
	sscanf(s.c_str(), "%d", &ip);
	j  = s.find("\"");
	s1 = s.substr(j+1, s.length()-j-2);
	physicalNames[ip] = s1;
	physNameIds[s1] = ip;
      }
      getline(f, s); // end tag "$EndPhysicalNames"
      getline(f, s);
    }
    if (s.compare("$Nodes") == 0)
    { // read nodes data
      getline(f, s);
      sscanf(s.c_str(), "%d", &nnodes);
      coords.resize(3*nnodes);

      for (i=0; i<nnodes; i++)
      {
        getline(f, s);
        sscanf(s.c_str(), "%d %lf %lf %lf", &j, &x, &y, &z);
        mapNodes[j] = i;
        setCoord(i, 0, x);
        setCoord(i, 1, y);
        setCoord(i, 2, z);
      }
  
      getline(f, s);
      sscanf(s.c_str(), "$EndNodes");
    }
    else
    {
      throw 3;
    }
  
    getline(f, s);
    if (s.compare("$Elements") != 0) throw 3;

    getline(f, s);
    sscanf(s.c_str(), "%d", &nelems);
  
    nvolumes = 0;
    nfaces   = 0;
    nedges   = 0;

    for (int i=0; i<nelems; i++)
    {
      f >> n >> type >> ntags;
      if (ntags >= 3) // read tags
      {
        f >> nphysical;
        f >> ngeom;
        f >> npart;
        for (j=3; j<ntags; f >> tag, j++);
      }
      else
      {
        npart     = 0;
        ngeom     = 0;
        nphysical = 0;
        for (j=0; j<ntags; f >> tag, j++);
      }
      switch (type)
      {
        case GMSH_TETRAHEDRON:
  	  tet.nphysical = nphysical;
	  tet.ngeom     = ngeom;
	  tet.npart     = npart;
	  for (j=0; j<4; j++) // read nodes
	  {
	    f >> n;
	    tet.nodes[j] = mapNodes[n];
	  }
	  volumes.push_back(tet);
	  nvolumes++;
          break;
        case GMSH_TRIANGLE:
          tri.nphysical = nphysical;
	  tri.ngeom     = ngeom;
	  tri.npart     = npart;
	  for (j=0; j<3; j++) // read nodes
	  {
	    f >> n;
	    tri.nodes[j] = mapNodes[n];
	  }
	  faces.push_back(tri);
	  nfaces++;
          break;
        case GMSH_LINE:
          lin.nphysical = nphysical;
	  lin.ngeom     = ngeom;
          lin.npart     = npart;
	  for (j=0; j<2; j++) // read nodes
 	  {
	    f >> n;
	    lin.nodes[j] = mapNodes[n];
	  }
	  edges.push_back(lin);
	  nedges++;
          break;
        case GMSH_POINT:
          f >> n;
	  break;
      }
    }
  }
  catch (int e)
  {
    switch (e)
    {
      case 1: cout << "Unsupported mesh file format.\n"; break;
      case 2: cout << "Data size does not match this system.\n"; break;
      case 3: cout << "Error reading mesh file.\n"; break;
    }
    f.close();
    mapNodes.clear();
    exit(1);
  }
  f.close();
  mapNodes.clear();
}

void C3DLinMesh::splitBoundaryVolumes()
{
  bool isBoundaryNode[nnodes];
  bool boundaryVolume;
  int i, j;
  
  for (i=0; i<nnodes; i++)
  {
    isBoundaryNode[i] = false;
  }
  
  // mark boundary nodes
  for (i=0; i<nfaces; i++)
  {
    for (j=0; j<3; j++)
    {
      isBoundaryNode[faces[i].nodes[j]] = true;
    }
  }
  
  // find volumes whose all vertices lie at the boundary
  int nvolumes_new = nvolumes;
  for (i=0; i<nvolumes; i++)
  {
    boundaryVolume = true;
    for (j=0; j<4; j++)
    {
      if (!isBoundaryNode[volumes[i].nodes[j]])
      {
        boundaryVolume = false;
	break;
      }
    }
    // split the volume
    if (boundaryVolume)
    {
      int node0 = volumes[i].nodes[0];
      int node1 = volumes[i].nodes[1];
      int node2 = volumes[i].nodes[2];
      int node3 = volumes[i].nodes[3];
      int new_node = nnodes;
      // add midpoint
      CVar x = (coords[3*node0  ]+coords[3*node1  ]+coords[3*node2  ]+coords[3*node3  ])*0.25;
      CVar y = (coords[3*node0+1]+coords[3*node1+1]+coords[3*node2+1]+coords[3*node3+1])*0.25;
      CVar z = (coords[3*node0+2]+coords[3*node1+2]+coords[3*node2+2]+coords[3*node3+2])*0.25;
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
      nnodes++;
      
      CTetrahedron tet(volumes[i]);
      // 0-1-2-4
      volumes[i].nodes[3] = new_node;
      // 0-1-4-3
      tet.nodes[2] = new_node;
      tet.nodes[3] = node3;
      volumes.push_back(tet);
      // 0-2-3-4
      tet.nodes[1] = node2;
      tet.nodes[2] = node3;
      tet.nodes[3] = new_node;
      volumes.push_back(tet);
      // 1-2-4-3
      tet.nodes[0] = node1;
      tet.nodes[2] = new_node;
      tet.nodes[3] = node3;
      volumes.push_back(tet);
      
      nvolumes_new += 3;
    }
  }
  nvolumes = nvolumes_new;
}



C3DQuadMesh::~C3DQuadMesh()
{
  volumes.clear();
  faces.clear();
  edges.clear();
  lnmap.clear();
}

void C3DQuadMesh::readGMSH(string fname)
{
  int i, j, ip, nelems, n, type, ntags, tag;
  int nphysical, ngeom, npart;
  double x, y, z;
  string s, s1;
  ifstream f(fname.c_str());
  map<int, int> mapNodes;
  CQuadTetrahedron tet;
  CQuadTriangle tri;
  CQuadLine lin;
  
  try
  {
    getline(f, s);
    if (s.compare("$MeshFormat") != 0) throw 1;
    getline(f, s);
    sscanf(s.c_str(), "2 0 %d", &i);
    if (i != sizeof(double)) throw 2;
    getline(f, s);
    if (s.compare("$EndMeshFormat") != 0) throw 1;
  
    getline(f, s);
    if (s.compare("$PhysicalNames") == 0)
    { // read through the physical names block
      getline(f, s);
      sscanf(s.c_str(), "%d", &n);
      for (i=0; i<n; i++)
      {
        getline(f, s);
	sscanf(s.c_str(), "%d", &ip);
	j  = s.find("\"");
	s1 = s.substr(j+1, s.length()-j-2);
	physicalNames[ip] = s1;
	physNameIds[s1] = ip;
      }
      getline(f, s); // end tag "$EndPhysicalNames"
      getline(f, s);
    }
    if (s.compare("$Nodes") == 0)
    { // read nodes data
      getline(f, s);
      sscanf(s.c_str(), "%d", &nnodes);
      coords.resize(3*nnodes);

      for (i=0; i<nnodes; i++)
      {
        getline(f, s);
        sscanf(s.c_str(), "%d %lf %lf %lf", &j, &x, &y, &z);
        mapNodes[j] = i;
        setCoord(i, 0, x);
        setCoord(i, 1, y);
        setCoord(i, 2, z);
      }
  
      getline(f, s);
      sscanf(s.c_str(), "$EndNodes");
    }
    else
    {
      throw 3;
    }
  
    getline(f, s);
    if (s.compare("$Elements") != 0) throw 3;

    getline(f, s);
    sscanf(s.c_str(), "%d", &nelems);
  
    nvolumes = 0;
    nfaces   = 0;
    nedges   = 0;
    nlnodes  = 0;
    lnmap.resize(nnodes);
    for (int i=0; i<nnodes; lnmap[i] = -1, i++);

    for (int i=0; i<nelems; i++)
    {
      f >> n >> type >> ntags;
      if (ntags >= 3) // read tags
      {
        f >> nphysical;
        f >> ngeom;
        f >> npart;
        for (j=3; j<ntags; f >> tag, j++);
      }
      else
      {
        npart     = 0;
        ngeom     = 0;
        nphysical = 0;
        for (j=0; j<ntags; f >> tag, j++);
      }
      switch (type)
      {
        case GMSH_QUAD_TETRAHEDRON:
  	  tet.nphysical = nphysical;
	  tet.ngeom     = ngeom;
	  tet.npart     = npart;
	  for (j=0; j<10; j++) // read nodes
	  {
	    f >> n;
	    tet.nodes[j] = mapNodes[n];
	  }
	  for (j=0; j<4; j++) // mark linear nodes
	  {
	    if (lnmap[tet.nodes[j]]==-1)
	    {
	      lnmap[tet.nodes[j]] = nlnodes;
	      nlnodes++;
	    }
	  }
	  volumes.push_back(tet);
	  nvolumes++;
          break;
        case GMSH_QUAD_TRIANGLE:
          tri.nphysical = nphysical;
	  tri.ngeom     = ngeom;
	  tri.npart     = npart;
	  for (j=0; j<6; j++) // read nodes
	  {
	    f >> n;
	    tri.nodes[j] = mapNodes[n];
	  }
	  for (j=0; j<3; j++) // mark linear nodes
	  {
	    if (lnmap[tri.nodes[j]]==-1)
	    {
	      lnmap[tri.nodes[j]] = nlnodes;
	      nlnodes++;
	    }
	  }
	  faces.push_back(tri);
	  nfaces++;
          break;
        case GMSH_QUAD_LINE:
          lin.nphysical = nphysical;
	  lin.ngeom     = ngeom;
          lin.npart     = npart;
	  for (j=0; j<3; j++) // read nodes
 	  {
	    f >> n;
	    lin.nodes[j] = mapNodes[n];
	  }
	  for (j=0; j<2; j++) // mark linear nodes
	  {
	    if (lnmap[lin.nodes[j]]==-1)
	    {
	      lnmap[lin.nodes[j]] = nlnodes;
	      nlnodes++;
	    }
	  }
	  edges.push_back(lin);
	  nedges++;
          break;
        case GMSH_POINT:
          f >> n;
	  break;
      }
    }
  }
  catch (int e)
  {
    switch (e)
    {
      case 1: cout << "Unsupported mesh file format.\n"; break;
      case 2: cout << "Data size does not match this system.\n"; break;
      case 3: cout << "Error reading mesh file.\n"; break;
    }
    f.close();
    mapNodes.clear();
    exit(1);
  }
  f.close();
  mapNodes.clear();
}


void C3DQuadMesh::writeGMSH(string fname)
{
  ofstream f(fname.c_str());
  
  f << "$MeshFormat\n";
  f << "2 0 " << sizeof(double) << endl;
  f << "$EndMeshFormat\n";
  if (getNphysicalNames() > 0)
  {
    f << "$PhysicalNames\n";
    f << getNphysicalNames() << endl;
    for (int i=0; i<getNphysicalNames(); i++)
    {
      f << i+1 << " \"" << physicalNames[i] << "\"\n";
    }
    f << "$EndPhysicalNames\n";
  }
  f << "$Nodes\n";
  f << nnodes << endl;
  for (int i=0; i<nnodes; i++)
  {
    f << i+1 << " ";
    for (int j=0; j<3; j++)
    {
      f << getCoord(i, j).getVal();
      if (j<2) f << " ";
    }
    f << endl;
  }
  f << "$EndNodes\n";
  f << "$Elements\n";
  f << nfaces + nvolumes << endl;
  int e=0;
  for (int i=0; i<nfaces; i++)
  {
    e++;
    f << e << " 9 3 " << getFacePhysicalNo(i) << " "
      << getFaceGeomNo(i) << " "
      << getFacePartNo(i) << " ";
    for (int j=0; j<6; j++)
    {
      f << getFaceNode(i, j)+1;
      if (j<5) f << " ";
    }
    f << endl;
  }
  for (int i=0; i<nvolumes; i++)
  {
    e++;
    f << e << " 11 3 " << getVolumePhysicalNo(i) << " "
      << getVolumeGeomNo(i) << " "
      << getVolumePartNo(i) << " ";
    for (int j=0; j<10; j++)
    {
      f << getVolumeNode(i, j)+1;
      if (j<9) f << " ";
    }
    f << endl;
  }
  f << "$EndElements\n";
  f.close();
}


void C3DQuadMesh::splitBoundaryVolumes()
{
  bool isBoundaryNode[nnodes];
  bool boundaryVolume;
  int i, j;
  
  for (i=0; i<nnodes; i++)
  {
    isBoundaryNode[i] = false;
  }
  
  // mark boundary nodes
  for (i=0; i<nfaces; i++)
  {
    for (j=0; j<6; j++)
    {
      isBoundaryNode[faces[i].nodes[j]] = true;
    }
  }
  
  // find volumes whose all vertices lie at the boundary
  int nvolumes_new = nvolumes;
  for (i=0; i<nvolumes; i++)
  {
    boundaryVolume = true;
    for (j=0; j<4; j++)
    {
      if (!isBoundaryNode[volumes[i].nodes[j]])
      {
        boundaryVolume = false;
	break;
      }
    }
    // split the volume
    if (boundaryVolume)
    {
      cout << "Volume " << i << " (" << volumes[i].nodes[0] << "," << volumes[i].nodes[1] << "," << volumes[i].nodes[2] << "," << volumes[i].nodes[3] << ") needs splitting.\n";
      int node[10];
      for (j=0; j<10; j++)
      {
        node[j] = volumes[i].nodes[j];
      }
      int noden  = nnodes;
      int node0n = nnodes+1;
      int node1n = nnodes+2;
      int node2n = nnodes+3;
      int node3n = nnodes+4;
      // add midpoint
      CVar x = (coords[3*node[0]  ]+coords[3*node[1]  ]+coords[3*node[2]  ]+coords[3*node[3]  ])*0.25;
      CVar y = (coords[3*node[0]+1]+coords[3*node[1]+1]+coords[3*node[2]+1]+coords[3*node[3]+1])*0.25;
      CVar z = (coords[3*node[0]+2]+coords[3*node[1]+2]+coords[3*node[2]+2]+coords[3*node[3]+2])*0.25;
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
      nnodes++;
      // add point between 0 and midpoint
      x = (coords[3*node[0]  ]+coords[3*noden  ])*0.5;
      y = (coords[3*node[0]+1]+coords[3*noden+1])*0.5;
      z = (coords[3*node[0]+2]+coords[3*noden+2])*0.5;
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
      nnodes++;
      // add point between 1 and midpoint
      x = (coords[3*node[1]  ]+coords[3*noden  ])*0.5;
      y = (coords[3*node[1]+1]+coords[3*noden+1])*0.5;
      z = (coords[3*node[1]+2]+coords[3*noden+2])*0.5;
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
      nnodes++;
      // add point between 2 and midpoint
      x = (coords[3*node[2]  ]+coords[3*noden  ])*0.5;
      y = (coords[3*node[2]+1]+coords[3*noden+1])*0.5;
      z = (coords[3*node[2]+2]+coords[3*noden+2])*0.5;
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
      nnodes++;
      // add point between 3 and midpoint
      x = (coords[3*node[3]  ]+coords[3*noden  ])*0.5;
      y = (coords[3*node[3]+1]+coords[3*noden+1])*0.5;
      z = (coords[3*node[3]+2]+coords[3*noden+2])*0.5;
      coords.push_back(x);
      coords.push_back(y);
      coords.push_back(z);
      nnodes++;
      // update linear nodes
      lnmap.resize(nnodes);
      lnmap[noden] = nlnodes;
      lnmap[node0n] = -1;
      lnmap[node1n] = -1;
      lnmap[node2n] = -1;
      lnmap[node3n] = -1;
      nlnodes++;
      
      CQuadTetrahedron tet = volumes[i];
      // 0-1-2-N
      volumes[i].nodes[3] = noden;
      volumes[i].nodes[7] = node0n;
      volumes[i].nodes[8] = node2n;
      volumes[i].nodes[9] = node1n;
      // 0-1-N-3
      tet.nodes[0] = node[0];
      tet.nodes[1] = node[1];
      tet.nodes[2] = noden;
      tet.nodes[3] = node[3];
      tet.nodes[4] = node[4];
      tet.nodes[5] = node1n;
      tet.nodes[6] = node0n;
      tet.nodes[7] = node[7];
      tet.nodes[8] = node3n;
      tet.nodes[9] = node[9];
      volumes.push_back(tet);
      // 0-2-3-N
      tet.nodes[0] = node[0];
      tet.nodes[1] = node[2];
      tet.nodes[2] = node[3];
      tet.nodes[3] = noden;
      tet.nodes[4] = node[6];
      tet.nodes[5] = node[8];
      tet.nodes[6] = node[7];
      tet.nodes[7] = node0n;
      tet.nodes[8] = node3n;
      tet.nodes[9] = node2n;
      volumes.push_back(tet);
      // 1-2-N-3
      tet.nodes[0] = node[1];
      tet.nodes[1] = node[2];
      tet.nodes[2] = noden;
      tet.nodes[3] = node[3];
      tet.nodes[4] = node[5];
      tet.nodes[5] = node2n;
      tet.nodes[6] = node1n;
      tet.nodes[7] = node[9];
      tet.nodes[8] = node3n;
      tet.nodes[9] = node[8];
      volumes.push_back(tet);
      
      nvolumes_new += 3;
    }
  }
  nvolumes = nvolumes_new;
}






