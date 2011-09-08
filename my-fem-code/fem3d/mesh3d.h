/* Finite element subroutines: mesh data structures
   Jan Stebel
*/

#ifndef _MESH3D_H__
#define _MESH3D_H__

#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <vector>

#include "ad.h"

#define GMSH_LINE             1
#define GMSH_TRIANGLE         2
#define GMSH_TETRAHEDRON      4
#define GMSH_QUAD_LINE        8
#define GMSH_QUAD_TRIANGLE    9
#define GMSH_QUAD_TETRAHEDRON 11
#define GMSH_POINT            15

#define VTK_LINE               3
#define VTK_TRIANGLE           5
#define VTK_TETRA              10
#define VTK_QUAD_EDGE     21
#define VTK_QUAD_TRIANGLE 22
#define VTK_QUAD_TETRA    24

/** VTK ordering of nodes differs from GMSH */
const int VTK_QUAD_TETRA_NODE_ORDER[] = { 0, 1, 2, 3, 4, 5, 6, 7, 9, 8 };


using namespace std;

/** Structure for elementary geometric entity */
class CElement
{
 public:
  int nphysical;  /**< Number of physical entity */
  int ngeom;      /**< Number of elementary geometrical entity */
  int npart;      /**< Number of mesh partition */
};

/** Structure for tetrahedra */
class CTetrahedron : public CElement
{
 public:
  int nodes[4]; /**< Array of node numbers. */
};

/** Structure for quadratic tetrahedra */
class CQuadTetrahedron : public CElement
{
 public:
  int nodes[10]; /**< Array of node numbers. */
};

/** Structure for triangles */
class CTriangle : public CElement
{
 public:
  int nodes[3]; /**< Array of node numbers. */
};

/** Structure for quadratic triangles */
class CQuadTriangle : public CElement
{
 public:
  int nodes[6]; /**< Array of node numbers. */
};

/** Structure for lines */
class CLine : public CElement
{
 public:
  int nodes[2]; /**< Array of node numbers. */
};

/** Structure for quadratic lines */
class CQuadLine : public CElement
{
 public:
  int nodes[3]; /**< Array of node numbers. */
};

/** Abstract class for 3D mesh. */
class C3DMesh
{
 protected:
  vector<CVar>    coords;        /**< Node coordinates.                   */
  map<int,string> physicalNames; /**< Names of physical entities.         */
  map<string,int> physNameIds;   /**< Reverse array of physical entities. */

  int     nnodes;    /**< Number of nodes.    */
  int     nedges;    /**< Number of 1D edges. */
  int     nfaces;    /**< Number of 2D faces. */
  int     nvolumes;  /**< Number of 3D cells. */
  
 public:
  /** Implicit constructor. */
  C3DMesh() {};
  
  /** Destructor. */
  virtual ~C3DMesh();
  
  /** Delete allocated arrays and structures. */
  virtual void destroy();
  
  /**
    Virtual method for reading mesh from file.
    @param fname File name.
  */
  virtual void readGMSH (string fname) = 0;
  
  /**
    Virtual method for writing mesh to file.
    @param fname File name.
  */
  virtual void writeGMSH(string fname) = 0;
  
  /**
    Deform the mesh coordinates according to given function.
    @param func Function performing deformation of comp-th coordinate
                of point located ad (x,y,z).
  */
  void         deform(CVar& func(double x, double y, double z, int comp));
  
  /**
    Return coordinate of mesh node.
    @param node  Node number.
    @param coord Coordinate number (0=x, 1=y, 2=z).
  */
  CVar& getCoord(int node, int coord);
  
  /**
    Set mesh node coordinate.
    @param node  Node number.
    @param coord Coordinate number.
    @param x     New coordinate value.
  */
  void  setCoord(int node, int coord, CVar x);
  
  /** Return number of physical entities. */
  int    getNphysicalNames()    { return physicalNames.size(); }
  
  /**
    Return physical name.
    @param n Number of physical entity.
  */
  string getPhysicalName(int n) { return physicalNames[n];     }
  
  /**
    Return number of physical entity given by name.
    @param name Name of physical entity.
  */
  int    getPhysNameId(string name) { return physNameIds[name]; }
  
  /** Iterator to beginning of physical names. */
  map<int,string>::iterator beginPhysicalNames() { return physicalNames.begin(); }
  
  /** Iterator to end of physical names. */
  map<int,string>::iterator endPhysicalNames() { return physicalNames.end(); }
  
  /**
    Calculate normal and tangent vectors on face.
    @param fid Face number.
    @param n   Normal vector.
    @param t1  First tangent vector.
    @param t2  Second tangent vector.
  */
  void   calculateFaceNormal(int fid, CVec &n, CVec &t1, CVec &t2);
  
  /**
    Calculate tangent vector on edge.
    @param eid Edge number.
    @param t   Tangent vector.
  */
  void   calculateEdgeTangent(int eid, CVec &t);
  
  /** Return number of nodes. */
  virtual int getNnodes()     { return nnodes; }
  
  /** Return number of linear nodes. */
  virtual int getNlnodes()    { return nnodes; }
  
  /** Return number of edges. */
  virtual int getNedges()     { return nedges; }
  
  /** Return number of faces. */
  virtual int getNfaces()     { return nfaces; }
  
  /** Return number of volumes. */
  virtual int getNvolumes()   { return nvolumes; }

  /**
    Return physical entity of given volume.
    @param v Volume number.
  */
  virtual int getVolumePhysicalNo(int v) = 0;
  
  /**
    Return physical entity of given face.
    @param f Face number.
  */
  virtual int getFacePhysicalNo  (int f) = 0;
  
  /**
    Return physical entity of given edge.
    @param e Edge number.
  */
  virtual int getEdgePhysicalNo  (int e) = 0;
  
  /**
    Return geometric entity of given volume.
    @param v Volume number.
  */
  virtual int getVolumeGeomNo(int v) = 0;
  
  /**
    Return geometric entity of given face.
    @param f Face number.
  */
  virtual int getFaceGeomNo  (int f) = 0;
  
  /**
    Return geometric entity of given edge.
    @param e Edge number.
  */
  virtual int getEdgeGeomNo  (int e) = 0;
  
  /**
    Return partition number of given volume.
    @param v Volume number.
  */
  virtual int getVolumePartNo(int v) = 0;
  
  /**
    Return partition number of given face.
    @param f Face number.
  */
  virtual int getFacePartNo  (int f) = 0;
  
  /**
    Return partition number of given edge.
    @param e Edge number.
  */
  virtual int getEdgePartNo  (int e) = 0;
  
  /**
    Return global index of given volume node.
    @param volume Volume number.
    @param node   Local node number.
  */
  
  virtual int getVolumeNode(int volume, int node) = 0;
  /**
    Return global index of given face node.
    @param face Face number.
    @param node Local node number.
  */
  
  virtual int getFaceNode  (int face, int node)   = 0;
  
  /**
    Return global index of given edge node.
    @param edge Edge number.
    @param node Local node number.
  */
  virtual int getEdgeNode  (int edge, int node)   = 0;
  
  /** Return number of d.o.f. per volume. */
  virtual int getNvolumeDofs() = 0;
  
  /** Return number of d.o.f. per face. */
  virtual int getNfaceDofs()   = 0;
  
  /** Return number of d.o.f. per edge. */
  virtual int getNedgeDofs()   = 0;
  
  /** Return VTK type of volume element. */
  virtual int getVtk3Dtype() = 0;
  
  /** Return VTK type of face element. */
  virtual int getVtk2Dtype() = 0;
  
  /** Return VTK type of edge element. */
  virtual int getVtk1Dtype() = 0;
  
  /**
    Return global index of volume node in VTK ordering.
    @param volume Volume number.
    @param node   Local node number (in VTK ordering).
  */
  virtual int getVtkVolumeNode(int volume, int node) { return getVolumeNode(volume, node); }
  
//  virtual CVar getFaceSize(int face) = 0;
//  virtual CVar getElemSize(int ind);
//  virtual CVar getElemDiam(int ind);
  
};


/** 3D Linear mesh. */
class C3DLinMesh : public C3DMesh
{
 private:
  vector<CTetrahedron> volumes;
  vector<CTriangle>    faces;
  vector<CLine>        edges;
  
 public:
  ~C3DLinMesh();
 
  void readGMSH(string fname);
  void writeGMSH(string fname) {}
  void splitBoundaryVolumes();
  
  int getVolumePhysicalNo(int v) { return volumes[v].nphysical; }
  int getFacePhysicalNo  (int f) { return faces[f].nphysical;   }
  int getEdgePhysicalNo  (int e) { return edges[e].nphysical;   }
  
  int getVolumeGeomNo(int v) { return volumes[v].ngeom; }
  int getFaceGeomNo  (int f) { return faces[f].ngeom;   }
  int getEdgeGeomNo  (int e) { return edges[e].ngeom;   }
  
  int getVolumePartNo(int v) { return volumes[v].npart; }
  int getFacePartNo  (int f) { return faces[f].npart;   }
  int getEdgePartNo  (int e) { return edges[e].npart;   }
  
  int getVolumeNode(int volume, int node) { return volumes[volume].nodes[node]; }
  int getFaceNode  (int face,   int node) { return faces[face].nodes[node]; }
  int getEdgeNode  (int edge,   int node) { return edges[edge].nodes[node]; }
  
  int getNvolumeDofs() { return 4; }
  int getNfaceDofs()   { return 3; }
  int getNedgeDofs()   { return 2; }
  
  int getVtk3Dtype() { return VTK_TETRA; }
  int getVtk2Dtype() { return VTK_TRIANGLE; }
  int getVtk1Dtype() { return VTK_LINE; }
  
//  CVar getFaceSize(int face) { return 0; }
  
};


/** 3D Quadratic mesh. */
class C3DQuadMesh : public C3DMesh
{
 private:
  vector<CQuadTetrahedron> volumes;
  vector<CQuadTriangle>    faces;
  vector<CQuadLine>        edges;
  
  int nlnodes;
  vector<int> lnmap;
  
 public:
  ~C3DQuadMesh();
 
  void readGMSH (string fname);
  void writeGMSH(string fname);
  void splitBoundaryVolumes ();
  
  int getVolumePhysicalNo(int v) { return volumes[v].nphysical; }
  int getFacePhysicalNo  (int f) { return faces[f].nphysical;   }
  int getEdgePhysicalNo  (int e) { return edges[e].nphysical;   }
  
  int getVolumeGeomNo(int v) { return volumes[v].ngeom; }
  int getFaceGeomNo  (int f) { return faces[f].ngeom;   }
  int getEdgeGeomNo  (int e) { return edges[e].ngeom;   }
  
  int getVolumePartNo(int v) { return volumes[v].npart; }
  int getFacePartNo  (int f) { return faces[f].npart;   }
  int getEdgePartNo  (int e) { return edges[e].npart;   }
  void setVolumePartNo(int v, int no) { volumes[v].npart = no; }
  void setFacePartNo  (int f, int no) { faces[f].npart = no;   }
  
  int getVolumeNode(int volume, int node) { return volumes[volume].nodes[node]; }
  int getFaceNode  (int face,   int node) { return faces[face].nodes[node]; }
  int getEdgeNode  (int edge,   int node) { return edges[edge].nodes[node]; }
  
  int getNvolumeDofs() { return 10; }
  int getNfaceDofs()   { return 6; }
  int getNedgeDofs()   { return 3; }
  
  int getNlnodes() { return nlnodes; }
  int getNoOfLnode(int node) { return lnmap[node]; }  // translate node no. within all nodes into linear node no.
  
  int getVtk3Dtype() { return VTK_QUAD_TETRA; }
  int getVtk2Dtype() { return VTK_QUAD_TRIANGLE; }
  int getVtk1Dtype() { return VTK_QUAD_EDGE; }
  int getVtkVolumeNode(int volume, int node) { return getVolumeNode(volume, VTK_QUAD_TETRA_NODE_ORDER[node]); }
  
//  CVar getFaceSize(int face) { return 0; }
  
};





#endif

