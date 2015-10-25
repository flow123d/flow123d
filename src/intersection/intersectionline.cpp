/*
 * IntersectionLine.cpp
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#include "intersectionline.h"
#include <armadillo>
#include "system/global_defs.h"


namespace computeintersection{
 
double IntersectionLine::compute_length()
{
    ASSERT(i_points.size() > 1, "Not enough intersetion points to define a line.");
    double length = 0;
    for(unsigned int i=0; i < i_points.size()-1; i++)
    {
        length += std::abs(i_points[i].get_local_coords1()[0] - i_points[i+1].get_local_coords1()[0]);
    }
    return length;
}

}