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
 * @file    mapping.hh
 * @brief   Class Mapping calculates data related to the mapping
 *          of the reference cell to the actual cell, such as Jacobian
 *          and normal vectors.
 * @author  Jan Stebel
 */

#ifndef MAPPING_HH_
#define MAPPING_HH_

#include <armadillo>
#include <vector>
#include "fem/dofhandler.hh"
#include "mesh/accessors.hh"
#include "fem/update_flags.hh"
#include "mesh/ref_element.hh"
#include "system/fmt/posix.h"           // for FMT_UNUSED



class Quadrature;
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

template<> inline double determinant(FMT_UNUSED const arma::mat::fixed<0,3> &M)
{
    return 0;
}

template<> inline double determinant(FMT_UNUSED const arma::mat::fixed<3,0> &M)
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

template<arma::uword n> inline double determinant(const arma::mat::fixed<n,n> &M)
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
    virtual MappingInternalData *initialize(const Quadrature &q, UpdateFlags flags) = 0;

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
    virtual void fill_fe_values(const ElementAccessor<3> &cell,
                        const Quadrature &q,
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
    virtual void fill_fe_side_values(const ElementAccessor<3> &cell,
                            unsigned int sid,
                            const Quadrature &q,
                            MappingInternalData &data,
                            FEValuesData<dim,spacedim> &fv_data) = 0;

    /// Destructor.
    virtual ~Mapping() {};
};

#endif /* MAPPING_HH_ */
