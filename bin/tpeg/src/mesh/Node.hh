#ifndef F123DTPEG_NODE_HH
#define F123DTPEG_NODE_HH

using namespace std;

#define COOR_X 0
#define COOR_Y 1
#define COOR_Z 2

class Node {
private:
    int label;
    double* coor;
    double level;

public:
    Node(int, double, double, double);
    ~Node();

    int getLabel();

    void setCoor(int, double);
    double getCoor(int);

    void setLevel(double);
    double getLevel();

};

#endif
