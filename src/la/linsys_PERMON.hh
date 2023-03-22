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
 * @file    linsys_PERMON.hh
 * @brief   PERMON QP solvers and FETI
 * @author  Jakub Kruzik
 */

#ifndef LA_LINSYS_PERMON_HH_
#define LA_LINSYS_PERMON_HH_

#include <functional>          // for unary_function
#include <string>              // for string
#include <vector>              // for vector
#include "la/linsys_PETSC.hh"  // for LinSys
#include "permonqps.h"         // for QPS and whole PERMON/PETSc stack
#include "permonqpfeti.h"      // ^except for FETI implementation

class Distribution;
namespace Input {
	class Record;
	namespace Type {
		class Record;
	}
}
namespace la {
    class BddcmlWrapper;
}

class LinSys_PERMON : public LinSys_PETSC
{

public:
	typedef LinSys FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    LinSys_PERMON(const  Distribution * rows_ds, const std::string &params = "");

    LinSys_PERMON(const DOFHandlerMultiDim &dh, const std::string &params = "");

    /**
     * Copy constructor.
     */
    LinSys_PERMON( LinSys_PERMON &other );

    void set_inequality(Mat matrix_ineq, Vec ineq);

    LinSys_PETSC::SolveInfo solve() override;

    /**
     * Returns information on absolute solver accuracy
     */
    inline double get_absolute_accuracy() override {
       return a_tol_;
    };

    void view(string text="") override;

    void set_from_input(const Input::Record in_rec) override;

    double get_solution_precision() override;

    double compute_residual() override;

    ~LinSys_PERMON( );

private:
    /// Registrar of class to factory
    static const int registrar;

    // make a pointer to the data array out of a std::vector
    template<typename T>
    T *  makePetscPointer_( std::vector<T> & array )
    {
        if ( array.size() ) return &(array[0]);
        return PETSC_NULL;
    }

    // PetscScalar to double casting functor
    struct PetscScalar2Double_ : public std::unary_function< PetscScalar, double >
    {
        double operator()( PetscScalar arg ) 
        {
            return static_cast<double>( arg );
        }
    };

protected:

    Mat     matrix_ineq_;        //!< PETSc matrix of inequality constraint.
    Vec     ineq_;               //!< PETSc vector of inequality constraint.
    Vec     warm_solution_;

    QP      system;
    QPS     solver;

    PetscReal maxeig_;

    bool    warm_start_;
};

#endif /* LA_LINSYS_PERMON_HH_ */
