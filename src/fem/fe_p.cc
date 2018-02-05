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
 * @file    fe_p.cc
 * @brief   
 */

// !! implementation of specializations has to be i *.cc file to avoid multiple definition error during linking
#include "fe_p.hh"
#include "mesh/ref_element.hh"




PolynomialSpace::PolynomialSpace(unsigned int degree, unsigned int dim)
    : degree_(degree)
{
// computes powers of all monomials up to given @p degree
// the order is: 1,x,x^2, y, yx,y^2
//
// TODO: - check and possibly rewrite to be more clear (use sum_degree temporary
//       - change order of monomials: 1, x, y, xy, x^2 , y^2 (increasing order)
//       - allow Q polynomials: 1,x, y, xy
//       - can use tensor products

    this->space_dim_ = dim;
    this->n_components_ = 1;
	arma::uvec pows(dim);
	pows.zeros();

    unsigned int degree_sum=0;
    unsigned int i_dim;


    while (true) {
        powers.push_back(pows);

        // increment pows
        for(i_dim=0; i_dim < dim; i_dim++) {
            if (degree_sum < degree) {
                pows[i_dim]++;
                degree_sum++;
                break;
            } else {                    // if degree_sum == degree, we find first non empty power, free it, and raise the next one
                degree_sum-=pows[i_dim];
                pows[i_dim]=0;
            }
        }
        if (i_dim == dim) break; // just after pow == (0, 0, .., degree)
    }
}


const double PolynomialSpace::basis_value(unsigned int i,
                                          const arma::vec &point,
                                          unsigned int comp_index
                                         ) const
{
    ASSERT(comp_index == 0);
	OLD_ASSERT(i<=powers.size(), "Index of basis function is out of range.");
    ASSERT(point.size()==space_dim_);

    double v = 1;
    for (unsigned int j=0; j<this->space_dim_; j++)
        v *= pow(point[j], (int) powers[i][j]);

    return v;
}


const arma::vec PolynomialSpace::basis_grad(unsigned int i,
                                            const arma::vec &p,
                                            unsigned int comp_index
                                           ) const
{
    ASSERT(comp_index == 0);
	OLD_ASSERT(i<=powers.size(), "Index of basis function is out of range.");

    arma::vec grad(this->space_dim_);

    for (unsigned int j=0; j<this->space_dim_; j++)
    {
        grad[j] = powers[i][j];
        if (powers[i][j] == 0) continue;

        for (unsigned int k=0; k<this->space_dim_; k++)
        {
            grad[j] *= pow(p[k], (int) (k==j?powers[i][k]-1:powers[i][k]));
        }
    }
    return grad;
}











template<unsigned int dim, unsigned int spacedim>
void FE_P<dim,spacedim>::init_dofs()
{
    if (degree_ == 0)
    {
        this->number_of_dofs = 1;
        this->number_of_single_dofs[dim] = 1;
        // we define nodal dof:
        // coords = barycentric coordinates of the support point,
        // coefs  = 1 (dof value = function value at the point)
        arma::vec coords = arma::ones<arma::vec>(dim+1)/(dim+1);
        this->dofs_.push_back(Dof(0, coords, { 1 }, Value));
    }
    else
    {
        // Create vector of barycentric coordinates.
        // First we make vector uvbc which contains barycentric coordinates
        // multiplied by degree_, so that its values are unsigned ints.
        // Then by counting the nonzero barycentric coordinates we can decide
        // whether the dof lies on node, line, triangle or tetrahedron.
        std::vector<arma::uvec> uvbc;
        arma::uvec ubc = arma::zeros<arma::uvec>(dim+1);
        ubc[dim] = degree_;
        bool finish = false;
        do {
            uvbc.push_back(ubc);
            if (ubc[dim] > 0)
            {
                // by default increment the first coordinate
                ubc[0] += 1;
                ubc[dim] -= 1;
            }
            else
            {
                // if sum of coordinates is maximal (last coordinate is zero)
                // then find first nonzero coordinate,
                // set it to zero, and increment the following coordinate.
                unsigned int c = 0;
                while (ubc[c] == 0) c++;
                // if the first nonzero coordinate is the last but one, we reach the end
                if (c == dim-1) finish = true;
                else {
                    ubc[dim] = ubc[c]-1;
                    ubc[c] = 0;
                    ubc[c+1] += 1;
                }
            }
        } while (!finish);
        
        // define dofs
        for (auto ubc : uvbc)
        {
            // we define nodal dof:
            // coords = barycentric coordinates of the support point,
            // coefs  = 1 (dof value = function value at the point)
            
            // count nonzero coordinates in ubc: 1=>nodal dof, 2=>dof on line etc.
            arma::uvec nonzeros = find(ubc);
            // convert "unsigned barycentric coordinates" to real ones
            arma::vec coords = arma::conv_to<arma::vec>::from(ubc);
            coords /= degree_;
            this->dofs_.push_back(Dof(nonzeros.size()-1, coords, { 1 }, Value));
        }
        
        this->number_of_dofs = this->dofs_.size();
    }
}




template<unsigned int dim, unsigned int spacedim>
FE_P<dim,spacedim>::FE_P(unsigned int degree)
  : FiniteElement<dim,spacedim>(),
    degree_(degree)
{
    this->function_space_ = new PolynomialSpace(degree,dim);
    
    init_dofs();

    this->setup_components();
    
    this->compute_node_matrix();
}















template<unsigned int dim, unsigned int spacedim>
FE_P_disc<dim,spacedim>::FE_P_disc(unsigned int degree)
    : FE_P<dim,spacedim>(degree)
{
    for (unsigned int i = 0; i <= dim; i++)
    {
        this->number_of_single_dofs[i] = 0;
        this->number_of_pairs[i] = 0;
        this->number_of_triples[i] = 0;
        this->number_of_sextuples[i] = 0;
    }
    this->number_of_single_dofs[dim] = this->number_of_dofs;

    for (auto dof : this->dofs_)
        dof.dim = dim;

    this->setup_components();

    this->compute_node_matrix();
}





template class FE_P_disc<0, 3>;
template class FE_P_disc<1, 3>;
template class FE_P_disc<2, 3>;
template class FE_P_disc<3, 3>;
