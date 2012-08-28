#ifndef nodeH
#define nodeH

#include <vector>
#include "point.h"

#define FOR_NODES(i, mesh)      for((i)=mesh->firstNode;(i)!=NULL;(i)=(i)->getNext())
#define FOR_NODE_SIDES(i,j)     for((j)=(i)->getSideBegin();(j)!=(i)->getSideEnd();(j)++)
#define FOR_NODE_ELEMENTS(i,j)  for((j)=(i)->getElementBegin();(j)!=(i)->getElementEnd();(j)++)

class TSide;
class TElement;

class TNode : public TPoint {
private:
    static int numberInstance;
    int id;

    int label;

    TNode* prev;
    TNode* next;

    std::vector<TSide*> sides;
    std::vector<TElement*> elements;

    int generateId();

    void Initialize();
    void setLabel(int);

public:
    TNode();
    TNode(double, double, double);
    TNode(int, double, double, double);
    ~TNode();

    void ParseLine(char*);
    int getLabel();

    void setPrev(TNode*);
    void setNext(TNode*);
    TNode* getNext();

    std::vector<TSide*>::iterator getSideBegin();
    std::vector<TSide*>::iterator getSideEnd();
    void addSide(TSide*);

    std::vector<TElement*>::iterator getElementBegin();
    std::vector<TElement*>::iterator getElementEnd();
    void addElement(TElement*);

    static int getNumInstances() {
        return TNode::numberInstance;
    }
};

#endif
