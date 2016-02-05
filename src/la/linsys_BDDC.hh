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
 * @file    linsys_BDDC.hh
 * @brief   Solver based on Multilevel BDDC - using corresponding class of OpenFTL package
 * @author  Jakub Sistek
 */

#ifndef LA_LINSYS_BDDC_HH_
#define LA_LINSYS_BDDC_HH_

// derived from base linsys
#include "mesh/mesh.h"
#include "la/linsys.hh"
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"

#include <vector>

namespace la {
    class BddcmlWrapper; 
};

class LinSys_BDDC : public LinSys
{

public:
	typedef LinSys FactoryBaseType;

    static const Input::Type::Record & get_input_type();

    LinSys_BDDC( const unsigned numDofsSub,
                 const Distribution * rows_ds,
                 const int matrixTypeInt = 0,
                 const int  numSubLoc = 1,
                 const bool swap_sign = false );

    void load_mesh( const int nDim, const int numNodes, const int numDofs,
                    const std::vector<int> & inet, 
                    const std::vector<int> & nnet, 
                    const std::vector<int> & nndf, 
                    const std::vector<int> & isegn, 
                    const std::vector<int> & isngn, 
                    const std::vector<int> & isvgvn,
                    const std::vector<double> & xyz,
                    const std::vector<double> & element_permeability,
                    const int meshDim );

    void mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals );

    void rhs_set_values( int nrow, int *rows, double *vals );

    void diagonal_weights_set_value( int global_index, double value );

    PetscErrorCode mat_zero_entries() override;

    PetscErrorCode rhs_zero_entries() override;

    void finish_assembly( );

    void apply_constrains( double scalar = 1. );

    int solve();

    void set_from_input(const Input::Record in_rec);

    double get_solution_precision();

    double compute_residual() override
    {
        // until we have correct function for residual we
        // return a practical infinity
        return numeric_limits<double>::max();
    }

    ~LinSys_BDDC( );

//private:

    //void gatherSolution_( );

private:
    /// Registrar of class to factory
    static const int registrar;

    //! parameters expected from input file:
    int  max_nondecr_it_;         //!< maximum number of iterations of linear solver with non-decreasing residual
    int  number_of_levels_;       //!< number of levels in the multilevel method
    bool use_adaptive_bddc_;      //!< should adaptive BDDC be used?
    int  bddcml_verbosity_level_; //!< level of verbosity of BDDCML library 
                                  //!< ( 0 - only fatal errors reported, 
                                  //!<   1 - mild output, 
                                  //!<   2 - detailed output )//!< should adaptive BDDC be used?
                                  
    const bool                        swap_sign_;      //!< swap sign of matrix and rhs entries, e.g. to make the matrix SPD

    std::vector<int>                  isngn_;          //!< indices of subdomain nodes in global numbering
    std::vector<double>               locSolution_;    //!< subdomain solution
    Vec                               locSolVec_;      //!< local solution PETSc vector - sequential
    VecScatter                        VSpetscToSubScatter_; //!< scatter from solution_ to locSolVec_

    typedef la::BddcmlWrapper         Bddcml_;
    Bddcml_ *                         bddcml_;         //!< BDDCML wrapper
};

#endif /* LA_LINSYS_BDDC_HH_ */
