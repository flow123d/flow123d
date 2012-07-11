#ifndef F123DTPEG_ETRIANGLE_HH
#define F123DTPEG_ETRIANGLE_HH

#include<vector>
#include"../mesh/Node.hh"
#include"../mesh/Element.hh"

using namespace std;

class ETriangle : public Element {
// =============================================================================
   private:
      double area;
// =============================================================================
   public:
      ETriangle(int, short, int*, Node**);
      ~ETriangle();
// -----------------------------------------------------------------------------
      void update();
      double getSize();
};

#endif
