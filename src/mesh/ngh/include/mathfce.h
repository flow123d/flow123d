#ifndef mainH
#define mainH

#include <limits>

const double epsilon = std::numeric_limits<double>::epsilon();


//#define epsilon 1e-6
namespace mathfce{
        bool IsZero(double);
        bool IsEqual(double, double);
}
        double Determinant3(double [3][3]);
#endif
