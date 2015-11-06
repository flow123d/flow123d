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
 * @file    matrix.h
 * @brief   
 */

#ifndef matrixH
#define matrixH

#include <iostream>

typedef enum {
    no_solution,
    one_solution,
    inf_solutions,
    badconditioned,
    singular
} TNSolutions;

class TMatrix {
private:
    int nc;
    int nr;

    double* elm;

public:
    TMatrix(int);
    TMatrix(const TMatrix &);
    TMatrix(int, int);
    ~TMatrix();

    int NRows() const;
    int NCols() const;

    void Set(int, int, double);
    double Get(int, int) const;

    void SwapRows(int, int);

    friend std::ostream & operator <<(std::ostream&, const TMatrix&);
};

class TMVector {
private:
    int size;

    double* elm;

public:
    TMVector(int);
    TMVector(const TMVector &);
    ~TMVector();

    void Set(int, double);
    double Get(int);
    void SwapElements(int, int);

    friend std::ostream & operator <<(std::ostream&, const TMVector&);
};

TNSolutions Gauss(const TMatrix&, TMVector*, const TMVector&);

#endif
