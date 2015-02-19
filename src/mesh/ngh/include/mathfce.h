#ifndef mainH
#define mainH

#include <limits>

const double epsilon = std::numeric_limits<double>::epsilon();


namespace mathfce{
        bool IsZero(double);
        bool IsEqual(double, double);
}
        double Determinant3(double [3][3]);
#endif
