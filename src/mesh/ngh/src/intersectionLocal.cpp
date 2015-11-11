/*!
 *
﻿ * Copyright (C) 2015 Technical University of Liberec.  All rights reserved.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3 as published by the
 * Free Software Foundation. (http://www.gnu.org/licenses/gpl-3.0.en.html)
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * 
 * @file    intersectionLocal.cpp
 * @brief   
 */

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <limits>

#include "system/exc_common.hh"
#include "mesh/ngh/include/intersectionLocal.h"

int IntersectionLocal::numberInstance = 0;

int IntersectionLocal::generateId() {
    return IntersectionLocal::numberInstance++;
}

bool eps_double_equal(double a, double b) {
	return fabs(a-b) < numeric_limits<double>::epsilon();
}

bool IntersectionPoint::operator ==(const IntersectionPoint &IP) {
	if (! equal (IP.coord1.begin(), IP.coord1.end(), el1_coord().begin(), eps_double_equal)) return false;
	if (! equal (IP.coord2.begin(), IP.coord2.end(), el2_coord().begin(), eps_double_equal)) return false;
	return true;
}

IntersectionPoint *interpolate(const IntersectionPoint &A1, const IntersectionPoint &A2, double t) {
    if (! (A1.el1_coord().size() == 1 && A2.el1_coord().size() == 1 ) ) {
        THROW( ExcAssertMsg() << EI_Message( "Interpolation of IntersectionPoints with non line first element.") );
    }
    if (! (A1.el2_coord().size() == A2.el2_coord().size()) ) {
        THROW( ExcAssertMsg() << EI_Message( "Interpolation of IntersectionPoints with non matching second element type.") );
    }
    std::vector<double> el2_coord(A1.el2_coord().size());
	for(unsigned int  i = 0; i < el2_coord.size(); ++i) {
    	el2_coord[i] = ((t - A1.el1_coord()[0]) / (A2.el1_coord()[0] - A1.el1_coord()[0])) * (A2.el2_coord()[i] - A1.el2_coord()[i]) + A1.el2_coord()[i];
	}
    return new IntersectionPoint(std::vector<double>(1,t), el2_coord);
}

IntersectionLocal::IntersectionLocal(IntersectionType i_type)
: type(i_type)
{
	id = generateId();
}


void IntersectionLocal::add_local_coord(const std::vector<double> &coordin1, const std::vector<double> &coordin2) {
	i_points.push_back(new IntersectionPoint(coordin1, coordin2));
}

void IntersectionLocal::add_local_point(IntersectionPoint *InPoint) {
	i_points.push_back(InPoint);
}

void IntersectionLocal::print(FILE *out_file) {
	int size_0 = i_points.size();
	fprintf(out_file," %d", size_0); //pocet bodu pruniku N

	// cyklus pres body pruniku
	for (std::vector<IntersectionPoint *>::iterator i_point = i_points.begin();
		i_point != i_points.end();++i_point) {

			// pocet lokalnich souradnic 1. elementu a jejich vypis
	        int size_1 = (*i_point)->coord1.size();
	        fprintf(out_file, " %d", size_1);
	    	for(std::vector<double>::iterator l_coord = (*i_point)->coord1.begin();
	    		l_coord != (*i_point)->coord1.end();++l_coord) {
	    			float f_1 = *l_coord;
	    			fprintf(out_file, " %1.7e", f_1);
	    	}

	    	// pocet lokalnich souradnic 2. elementu a jejich vypis
	    	int size_2 = (*i_point)->coord2.size();
	    	fprintf(out_file, " %d", size_2);
		    for(std::vector<double>::iterator l_coord = (*i_point)->coord2.begin();
		    	l_coord != (*i_point)->coord2.end();++l_coord) {
		    		float f_2 = *l_coord;
		    		fprintf(out_file, " %1.7e", f_2);
		    }
	}
}

IntersectionLocal::~IntersectionLocal() {
	for (std::vector<IntersectionPoint *>::iterator i_point = i_points.begin();
		i_point != i_points.end();++i_point) {
		delete *i_point;
	}
}
