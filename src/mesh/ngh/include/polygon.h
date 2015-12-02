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
 * @file    polygon.h
 * @brief   
 */

#ifndef polygonH
#define polygonH

#include "vertex.h"
#include "triangle.h"
#include <vector>

#define FOR_POL_VERTECES(i,j)     \
	for((j)=(i)->verteces.begin();(j)!=(i)->verteces.end();(j)++)

class TPolygon {
private:
    static int numberInstance;
    int id;

    double area;
    TPoint center;
    TVector normal_vector;
    bool area_is_actual;
    bool center_is_actual;

    std::vector<TVertex*> verteces;

    int generateId();

    void ComputeCenter();
    void ComputeArea();

    /**
     * Find position of new vertex
     *
     * @param Vx New vertex
     * @return position of new vertex
     */
    int InsertPosition(const TVertex& Vx);

public:
    TPolygon();
    ~TPolygon();

    friend std::ostream & operator <<(std::ostream&, const TPolygon&);

    void Add(const TPoint&);
    double GetArea();
    TPoint GetCenter();

    /**
     * Get count of vertexes
     */
    int vertexes_count() {
    	return verteces.size();
    }

    static int getNumInstances() {
        return TPolygon::numberInstance;
    }
};

#endif

