
#ifndef F123DTPEG_BOUNDARYDATA_HH
#define F123DTPEG_BOUNDARYDATA_HH

using namespace std;

typedef enum BoundaryDataType {
    DIRICHLET_PIEZOMETRIC_HEAD, // Piezometricka vyska [m]
    DIRICHLET_HEAD_PRESSURE,    // Tlakova vyska [m]
    NEUMANN_PER_SQUARE,         // Tok na jednotku plochy [(m^3 s^-1)/m^2]
    NEUMANN_SUM,                // Celkovy tok hranici [m^3 s^-1]
    NEWTON_PIEZOMETRIC_HEAD,    // Newtonova Piezometricka vyska [?]
    NEWTON_HEAD_PRESSURE,       // Newtonova Tlakova vyska [?]
    NUM_OF_BOUNDARY_TYPES
};

class BoundaryData {
private:
   void init(int, int, double, int);

   int     id;
   int     type;
   double  value;
   double  sigma;
   int     tag;

public:
   BoundaryData(int, int, double, int);
   BoundaryData(int, int, double, double, int);
   ~BoundaryData();

   int     getId();
   int     getType();
   double  getValue();
   double  getSigma();
   int     getTag();
};

#endif
