/*
 * IntersectionLine.cpp
 *
 *  Created on: 8.4.2015
 *      Author: viktor
 */

#include "intersectionline.h"
#include "intersectionpoint.h"
#include <armadillo>
#include "system/system.hh"


namespace computeintersection{
 
double IntersectionLine::compute_length()
{
    ASSERT(i_points_.size() > 1, "Not enough intersetion points to define a line.");
    double length = 0;
    for(unsigned int i=0; i < i_points_.size()-1; i++)
    {
        length += std::abs(i_points_[i].local_coords1()[0] - i_points_[i+1].local_coords1()[0]);
    }
    return length;
}

// void IntersectionLine::trace_line()
// {
//     if(i_points_.size() > 1){
//         unsigned int j = (i_points_[0].side_idx2() + i_points_[0].orientation())%2;
// 
//         if(j == 0){
//             std::vector<IntersectionPoint<1,3>> new_points(2);
//             new_points[0] = i_points_[1];
//             new_points[1] = i_points_[0];
//             i_points_ = new_points;
//         }
//     }
// }


}