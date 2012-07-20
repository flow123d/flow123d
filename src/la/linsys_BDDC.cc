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

#include <mpi.h>

// need OpenFTL BDDC wrapper
#include "la/SystemSolveBddc.hpp"

#include "la/linsys_BDDC.hh"


LinSys_BDDC::LinSys_BDDC( const unsigned lsize,
                          const unsigned numDofsSub,
                          const MPI_Comm comm, 
                          const int matrixTypeInt,
                          const int  numSubLoc )
        : LinSys( comm )
{
    // set type
    type = LinSys::BDDC;

    int lsizeInt = static_cast<int>( lsize );
    int numDofsInt;
    MPI_Allreduce ( &lsizeInt, &numDofsInt, 1, MPI_INT, MPI_SUM, comm_ );

    numDofs_ = static_cast<unsigned>( numDofsInt );

    la::SystemSolveBddc::MatrixType matrixType;
    switch ( matrixTypeInt ) {
        case 0:
            matrixType = la::SystemSolveBddc::GENERAL;
            break;
        case 1:
            matrixType = la::SystemSolveBddc::SPD;
            break;
        case 2:
            matrixType = la::SystemSolveBddc::SYMMETRICGENERAL;
            break;
        case 3:
            matrixType = la::SystemSolveBddc::SPD_VIA_SYMMETRICGENERAL;
            break;
        default:
            ASSERT( true, "Unknown matrix type %d", matrixTypeInt );
    }

    solver_ = new Solver_( numDofs_,
                           numDofsSub,
                           matrixType,
                           comm, 
                           numSubLoc );

    // set type
    type = LinSys::BDDC;
}

void LinSys_BDDC::load_mesh( Mesh *mesh,
                             Distribution *edge_ds,  
                             Distribution *el_ds,        
                             Distribution *side_ds,     
                             Distribution *rows_ds,    
                             int *el_4_loc,    
                             int *row_4_el,     
                             int *side_id_4_loc, 
                             int *side_row_4_id, 
                             int *edge_4_loc,   
                             int *row_4_edge )     
{
    mesh_          =  mesh;
    edge_ds_       =  edge_ds;
    el_ds_         =  el_ds;  
    side_ds_       =  side_ds;
    rows_ds_       =  rows_ds;
    el_4_loc_      =  el_4_loc;    
    row_4_el_      =  row_4_el;     
    side_id_4_loc_ =  side_id_4_loc;
    side_row_4_id_ =  side_row_4_id;
    edge_4_loc_    =  edge_4_loc;   
    row_4_edge_    =  row_4_edge;   

    // load mesh right away
    this -> loadFlowMesh_( );
}

void LinSys_BDDC::mat_set_values( int nrow, int *rows, int ncol, int *cols, double *vals )
{
    namespace ublas = boost::numeric::ublas;

    std::vector< unsigned >  myRows( nrow ); 
    std::vector< unsigned >  myCols( ncol ); 
    ublas::matrix< double >  mat( nrow, ncol ); 

    std::copy( &(rows[0]), &(rows[nrow]), myRows.begin() );
    std::copy( &(cols[0]), &(cols[ncol]), myCols.begin() );

    for ( unsigned i = 0; i < nrow; i++ ) {
        for ( unsigned j = 0; j < ncol; j++ ) {
            mat( i, j ) = vals[i*ncol + j];
        }
    }

    solver_ -> insertToMatrix( mat, myRows, myCols );
} 

void LinSys_BDDC::rhs_set_values( int nrow, int *rows, double *vals)
{
    namespace ublas = boost::numeric::ublas;

    std::vector< unsigned >  myRows( nrow ); 
    ublas::vector< double >  vec( nrow ); 

    std::copy( &(rows[0]), &(rows[nrow]), myRows.begin() );

    for ( unsigned i = 0; i < nrow; i++ ) {
        vec( i ) = vals[i];
    }

    solver_ -> insertToRhs( vec, myRows );
}

void LinSys_BDDC::finish_assembly( )
{
    solver_ -> finishMatAssembly( );
}

void LinSys_BDDC::apply_constrains( double scalar )
{
    solver_ -> applyConstraints( constraints_, 1., scalar );
}

int LinSys_BDDC::solve( )
{
    double              tol            = 1.e-7; //!< tolerance on relative residual ||res||/||rhs||
    int                 numLevels      = 2;     //!< number of levels
    std::vector<int> *  numSubAtLevels = NULL;  //!< number of subdomains at levels
    int                 verboseLevel   = 1;     //!< level of verbosity of BDDCML library 
                                                //!< ( 0 - only fatal errors reported, 
                                                //!<   1 - mild output, 
                                                //!<   2 - detailed output )
    int                 maxIt          = 1000;  //!< maximum number of iterations
    int                 ndecrMax       = 30;    //!< maximum number of iterations with non-decreasing residual 
                                                //!< ( used to stop diverging process )
    bool                use_adaptive   = false; //!< should adaptive BDDC be used?

    solver_ -> solveSystem( tol, numLevels, numSubAtLevels, verboseLevel, maxIt, ndecrMax, use_adaptive );

    DBGMSG("BDDCML converged reason: %d ( 0 means OK ) \n", solver_ -> giveConvergedReason() );
    DBGMSG("BDDCML converged in %d iterations. \n", solver_ -> giveNumIterations() );
    DBGMSG("BDDCML estimated condition number is %f \n", solver_ -> giveCondNumber() );

    return solver_ -> giveConvergedReason();
}

void LinSys_BDDC::get_whole_solution( std::vector<double> & globalSolution )
{
    this -> gatherSolution_( );
    globalSolution.resize( globalSolution_.size( ) );
    std::copy( globalSolution_.begin(), globalSolution_.end(), globalSolution.begin() );
}

void LinSys_BDDC::set_whole_solution( std::vector<double> & globalSolution )
{
    globalSolution_.resize( globalSolution.size( ) );
    std::copy( globalSolution.begin(), globalSolution.end(), globalSolution_.begin() );
}

LinSys_BDDC::~LinSys_BDDC()
{ 
    delete solver_; 
};

void LinSys_BDDC::loadFlowMesh_( )
{
    // initialize arrays
    std::map<int,arma::vec3> localDofMap;
    std::vector<int> inet;
    std::vector<int> nnet;
    std::vector<int> isegn;
    for ( int i_loc = 0; i_loc < el_ds_->lsize(); i_loc++ ) {
        // for each element, create local numbering of dofs as fluxes (sides), pressure (element centre), Lagrange multipliers (edges), compatible connections
        Element *el = mesh_->element(el_4_loc_[i_loc]);
        int e_idx = ELEMENT_FULL_ITER( mesh_, el ).index();

        isegn.push_back( e_idx );
        int nne = 0;

        FOR_ELEMENT_SIDES(el,si) {
            // insert local side dof
            int side_row = side_row_4_id_[el->side[si]->id];
            arma::vec3 coord = (el->side[si])->centre;

            localDofMap.insert( std::make_pair( side_row, coord ) );
            inet.push_back( side_row );
            nne++;
        }

        // insert local pressure dof
        int el_row  = row_4_el_[ el_4_loc_[i_loc] ];
        arma::vec3 coord = el->centre;
        localDofMap.insert( std::make_pair( el_row, coord ) );
        inet.push_back( el_row );
        nne++;

        FOR_ELEMENT_SIDES(el,si) {
            Edge *edg=el->side[si]->edge; 

            // insert local edge dof
            int edge_row = row_4_edge_[mesh_->edge.index(el->side[si]->edge)];
            arma::vec3 coord = (el->side[si])->centre;

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        // insert dofs related to compatible connections
        for ( int i_neigh = 0; i_neigh < el->n_neighs_vb; i_neigh++) {
            int edge_row = row_4_edge_[mesh_->edge.index(el->neigh_vb[i_neigh]->edge)];
            arma::vec3 coord = ( el->neigh_vb[i_neigh]->edge->side[0] )->centre;

            localDofMap.insert( std::make_pair( edge_row, coord ) );
            inet.push_back( edge_row );
            nne++;
        }

        nnet.push_back( nne );
    }
    //convert set of dofs to vectors
    int numNodeSub = localDofMap.size();
    isngn_.resize(  numNodeSub );
    std::vector<double> xyz( numNodeSub * 3 ) ;
    int ind = 0;
    std::map<int,arma::vec3>::iterator itB = localDofMap.begin();
    for ( ; itB != localDofMap.end(); ++itB ) {
        isngn_[ind] = itB -> first;

        arma::vec3 coord = itB -> second;
        for ( int j = 0; j < 3; j++ ) {
            xyz[ j*numNodeSub + ind ] = coord[j];
        }

        ind++;
    }
    localDofMap.clear();

    // nndf is trivially one
    std::vector<int> nndf( numNodeSub, 1 );

    // prepare auxiliary map for renumbering nodes 
    typedef std::map<int,int> Global2LocalMap_; //! type for storage of global to local map
    Global2LocalMap_ global2LocalNodeMap;
    for ( unsigned ind = 0; ind < isngn_.size(); ++ind ) {
        global2LocalNodeMap.insert( std::make_pair( static_cast<unsigned>( isngn_[ind] ), ind ) );
    }

    //std::cout << "INET: \n";
    //std::copy( inet.begin(), inet.end(), std::ostream_iterator<int>( std::cout, " " ) );
    //std::cout << std::endl;
    //std::cout << "ISNGN: \n";
    //std::copy( isngn_.begin(), isngn_.end(), std::ostream_iterator<int>( std::cout, " " ) );
    //std::cout << std::endl << std::flush;
    //std::cout << "ISEGN: \n";
    //std::copy( isegn.begin(), isegn.end(), std::ostream_iterator<int>( std::cout, " " ) );
    //std::cout << std::endl << std::flush;
    //MPI_Barrier( PETSC_COMM_WORLD );

    // renumber nodes in the inet array to locals
    int indInet = 0;
    for ( int iEle = 0; iEle < isegn.size(); iEle++ ) {
        int nne = nnet[ iEle ];
        for ( unsigned ien = 0; ien < nne; ien++ ) {

            int indGlob = inet[indInet];
            // map it to local node
            Global2LocalMap_::iterator pos = global2LocalNodeMap.find( indGlob );
            ASSERT( pos != global2LocalNodeMap.end(),
                    "Cannot remap node index %d to local indices. \n ", indGlob );
            int indLoc = static_cast<int> ( pos -> second );

            // store the node
            inet[ indInet++ ] = indLoc;
        }
    }

    int numNodes    = numDofs_;
    int numDofsInt  = numDofs_;
    int spaceDim = 3; // TODO: what is the proper value here?
    int meshDim  = 1; // TODO: what is the proper value here?

    solver_ -> loadRawMesh( spaceDim, numNodes, numDofsInt, inet, nnet, nndf, isegn, isngn_, isngn_, xyz, meshDim );
}

// construct global solution
void LinSys_BDDC::gatherSolution_( )
{
    std::vector<double> sol_disordered( numDofs_ );
    
    // download local solution
    std::vector<double> locSolution( isngn_.size() );
    solver_ -> giveSolution( isngn_, locSolution ); 
    int ierr;

    // merge solution on root
    int rank;
    MPI_Comm_rank( comm_, &rank );
    int nProc;
    MPI_Comm_size( comm_, &nProc );

    if ( rank == 0 ) {
        // merge my own data
        for ( int i = 0; i < isngn_.size(); i++ ) {
            int ind = isngn_[i];
            sol_disordered[ind] = locSolution[i];
        }
        for ( int iProc = 1; iProc < nProc; iProc++ ) {
            // receive length
            int length;
            MPI_Status status;
            ierr = MPI_Recv( &length, 1, MPI_INT, iProc, iProc, comm_, &status ); 

            // receive indices
            std::vector<int> inds( length );
            ierr = MPI_Recv( &(inds[0]), length, MPI_INT, iProc, iProc, comm_, &status ); 

            // receive local solution
            locSolution.resize( length );
            ierr = MPI_Recv( &(locSolution[0]), length, MPI_DOUBLE, iProc, iProc, comm_, &status ); 

            // merge data
            for ( int i = 0; i < length; i++ ) {
                int ind = inds[i];
                sol_disordered[ind] = locSolution[i];
            }
        }
    }
    else {
        // send my solution to root
        int length = isngn_.size();
        ierr = MPI_Send( &length,                1, MPI_INT,    0, rank, comm_ ); 
        ierr = MPI_Send( &(isngn_[0]),      length, MPI_INT,    0, rank, comm_ ); 
        ierr = MPI_Send( &(locSolution[0]), length, MPI_DOUBLE, 0, rank, comm_ ); 
    }
    // broadcast global solution from root
    ierr = MPI_Bcast( &(sol_disordered[0]), sol_disordered.size(), MPI_DOUBLE, 0, comm_ );

    //reorder solution
    std::vector<unsigned> indices;
    this->create_renumbering_( indices );
    globalSolution_.resize( numDofs_ );
    for ( int i = 0; i < numDofs_; i++ ) {
        globalSolution_[i] = sol_disordered[indices[i]];
    }

    //if ( myp == 0 ) {
    //    std::copy( solution_.begin(), solution_.end(), std::ostream_iterator<double>( std::cout, " " ) );
    //    std::cout << std::endl;
    //}
}


