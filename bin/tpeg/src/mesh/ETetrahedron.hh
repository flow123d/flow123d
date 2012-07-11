#ifndef F123DTPEG_ETETRAHEDRON_HH
#define F123DTPEG_ETETRAHEDRON_HH

#include<vector>
#include"../mesh/Node.hh"
#include"../mesh/Element.hh"

using namespace std;

class ETetrahedron : public Element {
// =============================================================================
   private:
      double volume;
// =============================================================================
   public:
      ETetrahedron(int, short, int*, Node**);
      ~ETetrahedron();
// -----------------------------------------------------------------------------
      void update();
      double getSize();
};

#endif
