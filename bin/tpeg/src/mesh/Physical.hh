#ifndef F123DTPEG_PHYSICAL_HH
#define F123DTPEG_PHYSICAL_HH

#define PHYSICAL_GROUP 0
#define ELEMENT_GROUP 1

using namespace std;

class Physical {
private:
      int       dim;
      int       id;
      char*     name;

   public:
      Physical(int, int, char*);
      ~Physical();

      int getDim();
      int getId();
      char* getName();

};

#endif
