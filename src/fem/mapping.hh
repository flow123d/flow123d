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
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Class Mapping calculates data related to the mapping
 *        of the reference cell to the actual cell, such as Jacobian
 *        and normal vectors.
 * @author Jan Stebel
 */

#ifndef MAPPING_HH_
#define MAPPING_HH_

#include <armadillo>
#include <vector>
#include "fem/dofhandler.hh"
#include "mesh/side_impl.hh"
#include "fem/update_flags.hh"
#include "mesh/ref_element.hh"



template<unsigned int dim> class Quadrature;
template<unsigned int dim, unsigned int spacedim> class FEValuesData;




/**
 * @brief Calculates determinant of a rectangular matrix.
 */
template<class T>
double determinant(const T &M);


template<> inline double determinant(const arma::mat::fixed<1,2> &M)
{
    return sqrt(M(0,0)*M(0,0)+M(0,1)*M(0,1));
}

template<> inline double determinant(const arma::mat::fixed<2,1> &M)
{
    return sqrt(M(0,0)*M(0,0)+M(1,0)*M(1,0));
}

template<> inline double determinant(const arma::mat::fixed<0,3> &M)
{
    return 0;
}

template<> inline double determinant(const arma::mat::fixed<3,0> &M)
{
    return 0;
}

template<> inline double determinant(const arma::mat::fixed<1,3> &M)
{
    return sqrt(M(0,0)*M(0,0)+M(0,1)*M(0,1)+M(0,2)*M(0,2));
}

template<> inline double determinant(const arma::mat::fixed<3,1> &M)
{
    return sqrt(M(0,0)*M(0,0)+M(1,0)*M(1,0)+M(2,0)*M(2,0));
}

template<> inline double determinant(const arma::mat::fixed<2,3> &M)
{
    return sqrt((M(0,0)*M(0,0)+M(0,1)*M(0,1)+M(0,2)*M(0,2))*(M(1,0)*M(1,0)+M(1,1)*M(1,1)+M(1,2)*M(1,2))
               -(M(0,0)*M(1,0)+M(0,1)*M(1,1)+M(0,2)*M(1,2))*(M(0,0)*M(1,0)+M(0,1)*M(1,1)+M(0,2)*M(1,2)));
}

template<> inline double determinant(const arma::mat::fixed<3,2> &M)
{
    return sqrt((M(0,0)*M(0,0)+M(1,0)*M(1,0)+M(2,0)*M(2,0))*(M(0,1)*M(0,1)+M(1,1)*M(1,1)+M(2,1)*M(2,1))
               -(M(0,0)*M(0,1)+M(1,0)*M(1,1)+M(2,0)*M(2,1))*(M(0,0)*M(0,1)+M(1,0)*M(1,1)+M(2,0)*M(2,1)));
}

template<unsigned int n> inline double determinant(const arma::mat::fixed<n,n> &M)
{
    return det(M);
}


/**
 * @brief Mapping data that can be precomputed on the actual cell.
 *
 * So far this involves only the (local) barycentric coordinates of quadrature points.
 */
class MappingInternalData
{
public:
    /**
     * @brief Auxiliary array of barycentric coordinates of quadrature points.
     */
    std::vector<arma::vec> bar_coords;
};



/**
 * @brief Abstract class for the mapping between reference and actual cell.
 *
 * Class Mapping calculates data related to the mapping of the
 * reference cell to the actual cell, such as Jacobian and normal
 * vectors.
 */
template<unsigned int dim, unsigned int spacedim>
class Mapping
{
public:

    /**
     * @brief Calculates the mapping data on the reference cell.
     *
     * @param q Quadrature rule.
     * @param flags Update flags.
     */
    virtual MappingInternalData *initialize(const Quadrature<dim> &q, UpdateFlags flags) = 0;

    /**
     * @brief Decides which additional quantities have to be computed
     * for each cell.
     *
     * @param flags Flags of required quantities.
     * @return Flags of all necessary quantities.
     */
    virtual UpdateFlags update_each(UpdateFlags flags) = 0;

    /**
     * @brief Calculates the mapping data and stores them in the provided
     * structures.
     *
     * @param cell The actual cell.
     * @param q Quadrature rule.
     * @param data Precomputed mapping data.
     * @param fv_data Data to be computed.
     */
    virtual void fill_fe_values(const typename DOFHandler<dim,spacedim>::CellIterator &cell,
                        const Quadrature<dim> &q,
                        MappingInternalData &data,
                        FEValuesData<dim,spacedim> &fv_data) = 0;

    /**
     * @brief Calculates the mapping data related to a given side, namely the
     * jacobian determinants and the normal vectors.
     *
     * @param cell The actual cell.
     * @param side Number of the side.
     * @param q Quadrature rule.
     * @param data Precomputed mapping data.
     * @param fv_data Data to be computed.
     */
    virtual void fill_fe_side_values(const typename DOFHandler<dim,spacedim>::CellIterator &cell,
                            unsigned int sid,
                            const Quadrature<dim> &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data) = 0;

    /**
     * @brief Creates a cell dim-dimensional quadrature from side (dim-1)-dimensional quadrature.
     *
     * @param cell The actual cell.
     * @param sid  The side id.
     * @param pid  Permutations index.
     * @param subq @p dim-1 dimensional quadrature for integration on the side.
     * @param q    The computed @p dim dimensional quadrature.
     */
    void transform_subquadrature(unsigned int sid,
    		unsigned int pid,
    		const Quadrature<dim-1> &subq,
    		Quadrature<dim> &q);

    /// Destructor.
    virtual ~Mapping() {};
};



template<unsigned int dim, unsigned int spacedim> inline
void Mapping<dim,spacedim>::transform_subquadrature(unsigned int sid,
		unsigned int pid,
        const Quadrature<dim - 1> & subq,
        Quadrature<dim> &q)
{
//    ASSERT(side.dim()==dim-1, "Side dimension mismatch.");
    q.resize(subq.size());

    double lambda;

    // vectors of barycentric coordinates of quadrature points
    arma::vec::fixed<dim+1> el_bar_coords;
    arma::vec::fixed<dim> side_bar_coords;

    for (unsigned int k=0; k<q.size(); k++)
    {
        // Calculate barycentric coordinates on the side of the k-th
        // quadrature point.
        el_bar_coords.zeros();
        lambda = 0;
        // Apply somewhere permutation of indices!
        for (int j=0; j<dim-1; j++)
        {
            side_bar_coords(j) = (subq.point(k))[j];
            lambda += (subq.point(k))[j];
        }
        side_bar_coords(dim-1) = 1.-lambda;

        // transform to element coordinates
        for (unsigned int i=0; i<dim; i++)
            el_bar_coords((RefElement<dim>::side_nodes[sid][RefElement<dim>::side_permutations[pid][i]]+dim)%(dim+1)) = side_bar_coords((i+dim-1)%dim);
        q.set_point(k, el_bar_coords.subvec(0,dim-1));
        q.set_weight(k, subq.weight(k));
    }
}












#endif /* MAPPING_HH_ */
