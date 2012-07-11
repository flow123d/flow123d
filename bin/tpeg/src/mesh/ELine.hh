#ifndef F123DTPEG_ELINE_HH
#define F123DTPEG_ELINE_HH

#include<vector>
#include"../mesh/Node.hh"
#include"../mesh/Element.hh"

using namespace std;

class ELine : public Element {
// =============================================================================
   private:
      double length;
// =============================================================================
   public:
      ELine(int, short, int*, Node**);
      ~ELine();
// -----------------------------------------------------------------------------
      void update();
      double getSize();
};

#endif
