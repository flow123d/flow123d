
#ifndef MESH_HH_
#define MESH_HH_

#include <vector>

class Point;
  
/**
 * Mesh of given dimension. Always assumed in 3d ambient space.
 */
template <int dim>
class Mesh {
public:  
    /**
     * Default constructor just constructs an empty mesh.
     */ 
    Mesh();
private:
    std::vector< Point >  nodes;
  
};  


#include "new_mesh/mesh.impl.hh"

#endif