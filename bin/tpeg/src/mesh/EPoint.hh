#ifndef F123DTPEG_EPOINT_HH
#define F123DTPEG_EPOINT_HH

#include<vector>
#include"../mesh/Node.hh"
#include"../mesh/Element.hh"

using namespace std;

class EPoint : public Element {
// =============================================================================
   private:
      double size;
// =============================================================================
   public:
      EPoint(int, short, int*, Node**);
      ~EPoint();
// -----------------------------------------------------------------------------
      void update();
      double getSize();
};

#endif
