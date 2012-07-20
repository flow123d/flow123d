/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id: la_linsys.hh 1299 2011-08-23 21:42:50Z jakub.sistek $
 * $Revision: 1299 $
 * $LastChangedBy: jakub.sistek $
 * $LastChangedDate: 2011-08-23 23:42:50 +0200 (Tue, 23 Aug 2011) $
 *
 * @file
 * @brief   Solver based on Multilevel BDDC - using corresponding class of OpenFTL package
 * @author  Jakub Sistek
 *
 *
 */

#ifndef LA_LINSYS_BDDC_HH_
#define LA_LINSYS_BDDC_HH_

// derived from base linsys
#include "mesh/mesh.h"
#include "la/linsys.hh"

#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace la{
    class SystemSolveBddc; 
};

class LinSys_BDDC : public LinSys
{

public:

    LinSys_BDDC( const unsigned lsize,
                 const unsigned numDofsSub,
                 const MPI_Comm comm = MPI_COMM_WORLD, 
                 const int matrixTypeInt = 0,
                 const int  numSubLoc = 1 );

    void load_mesh( Mesh *mesh,
                    Distribution *edge_ds,  
                    Distribution *el_ds,        
                    Distribution *side_ds,     
                    Distribution *rows_ds,    
                    int *el_4_loc,    
                    int *row_4_el,     
                    int *side_id_4_loc, 
                    int *side_row_4_id, 
                    int *edge_4_loc,   
                    int *row_4_edge );

    void mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals );

    void rhs_set_values( int nrow, int *rows, double *vals );

    void finish_assembly( );

    void apply_constrains( double scalar = 1. );

    int solve( );

    void get_whole_solution( std::vector<double> & globalSolution );

    void set_whole_solution( std::vector<double> & globalSolution );

    ~LinSys_BDDC( );

private:

    void loadFlowMesh_( );

    void gatherSolution_( );

private:

    std::vector<double>               globalSolution_; //!< global solution in numbering for linear system
    unsigned                          numDofs_;        //!< number of degrees of freedom on subdomain

    std::vector<int>                  isngn_;

    typedef la::SystemSolveBddc       Solver_;
    Solver_ *                         solver_;         //!< OpenFTL solver

};

#endif /* LA_LINSYS_BDDC_HH_ */
