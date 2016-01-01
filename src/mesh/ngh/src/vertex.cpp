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
 * @file    vertex.cpp
 * @brief   
 */

#include "mesh/ngh/include/vertex.h"
#include "mesh/ngh/include/myvector.h"

int TVertex::numberInstance = 0;

int TVertex::generateId() {
    return TVertex::numberInstance++;
}

TVertex::TVertex(const TPoint& PP) {
    id = generateId();

    X = new TPoint(PP);
}

TVertex::~TVertex() {
    delete X;
}

TPoint TVertex::GetPoint() const {
    TPoint tmp;
    tmp = *X;
    return tmp;
}
