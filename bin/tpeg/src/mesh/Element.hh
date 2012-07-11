#ifndef F123DTPEG_ELEMENT_HH
#define F123DTPEG_ELEMENT_HH

#include<vector>
#include"../mesh/Node.hh"

#define ELEMENT_TYPE_POINT       15
#define ELEMENT_TYPE_LINE         1
#define ELEMENT_TYPE_TRIANGLE     2
#define ELEMENT_TYPE_TETRAHEDRON  4

using namespace std;

class Element {
   static const short NUMBER_OF_NODES[];
// =============================================================================
   private:
      int    label;
      short  elementType;
      short  numTags;

      int*   tags;
      Node** nodes;

      bool   swIsFictive;

// =============================================================================
   public:
      Node* T;
// -----------------------------------------------------------------------------
      static short getNumNodes(short);
// -----------------------------------------------------------------------------
      Element(short, int, short, int*, Node**);
      ~Element();
// -----------------------------------------------------------------------------
      int getLabel();

      short getElementType();

      short getNumTags();
      int getTag(short);

      short getNumNodes();
      Node* getNode(short);

      void setFictive(bool);
      bool isFictive();

// -----------------------------------------------------------------------------
      virtual void update() = 0;
      virtual double getSize() = 0;
};

#endif
