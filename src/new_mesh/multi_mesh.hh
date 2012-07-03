#ifndef MULTI_MESH_HH_
#define MULTI_MESH_HH_
  
  
// forward declarations  
template <int dim>
class Mesh;

template <int dim>
class Matching;


/**
 * Mesh of multiple dimensions with:
 * - matching - Cell<dim> match face of an Cell<dim+1>
 * - intersection - arbitrary
 * 
 * Future: 
 * - have arbitrary number of meshes of every dimension
 */
class MultiMesh {
public:
  MultiMesh();
  
private:  
  Mesh<1>       mesh_1_;
  Mesh<2>       mesh_2_;
  Mesh<3>       mesh_3_;
  
  Matching<1>   match_1_;
  Matching<2>   match_2_;
  
};

#endif

