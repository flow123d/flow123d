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
    : FunctionSpace(dim, 1),
      degree_(degree)
{
// computes powers of all monomials up to given @p degree
// the order is: 1,x,x^2, y, yx,y^2
//
// TODO: - check and possibly rewrite to be more clear (use sum_degree temporary
//       - change order of monomials: 1, x, y, xy, x^2 , y^2 (increasing order)
//       - allow Q polynomials: 1,x, y, xy
//       - can use tensor products

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


double PolynomialSpace::basis_value(unsigned int i,
                                    const arma::vec &point,
                                    unsigned int comp_index
                                    ) const
{
    ASSERT_EQ_DBG(comp_index, 0);
	ASSERT_LE_DBG(i, powers.size());
    ASSERT_EQ_DBG(point.size(), space_dim_);

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
    ASSERT_EQ_DBG(comp_index, 0);
	ASSERT_LE_DBG(i, powers.size());

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











template<unsigned int dim>
void FE_P<dim>::init_dofs()
{
    if (degree_ == 0 || dim == 0)
    {
        // we define nodal dof:
        // coords = barycentric coordinates of the support point,
        // coefs  = 1 (dof value = function value at the point)
        arma::vec coords = arma::ones<arma::vec>(dim+1)/(dim+1);
        this->dofs_.push_back(Dof(dim, 0, coords, { 1 }, Value));
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
        ubc[0] = degree_;
        bool finish = false;
        do {
            uvbc.push_back(ubc);
            if (ubc[0] > 0)
            {
                // by default increment the first coordinate
                ubc[1] += 1;
                ubc[0] -= 1;
            }
            else
            {
                // if sum of coordinates is maximal (0-th coordinate is zero)
                // then find first nonzero coordinate,
                // set it to zero, and increment the following coordinate.
                unsigned int c = 1;
                while (ubc[c] == 0) c++;
                // if the first nonzero coordinate is the last but one, we reach the end
                if (c == dim) finish = true;
                else {
                    ubc[0] = ubc[c]-1;
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
            
            // find index of n-face
            std::pair<unsigned int, unsigned int> zeros = RefElement<dim>::zeros_positions(coords);
            unsigned int n_face_idx = -1;
            switch (dim-zeros.first) {
                case 0:
                    n_face_idx = RefElement<dim>::template topology_idx<0>(zeros.second);
                    break;
                case 1:
                    n_face_idx = RefElement<dim>::template topology_idx<1>(zeros.second);
                    break;
                case 2:
                    n_face_idx = RefElement<dim>::template topology_idx<2>(zeros.second);
                    break;
                case 3:
                    n_face_idx = RefElement<dim>::template topology_idx<3>(zeros.second);
                    break;
            }
            this->dofs_.push_back(Dof(nonzeros.size()-1, n_face_idx, coords, { 1 }, Value));
        }
    }
}




template<unsigned int dim>
FE_P<dim>::FE_P(unsigned int degree)
  : FiniteElement<dim>(),
    degree_(degree)
{
    this->function_space_ = std::make_shared<PolynomialSpace>(degree,dim);
    
    init_dofs();

    this->setup_components();
    
    this->compute_node_matrix();
}















template<unsigned int dim>
FE_P_disc<dim>::FE_P_disc(unsigned int degree)
    : FE_P<dim>(degree)
{
    // all dofs are "inside" the cell (not shared with neighbours)
    for (unsigned int i=0; i<this->dofs_.size(); i++)
        this->dofs_[i].dim = dim;
    
    this->setup_components();

    this->compute_node_matrix();
}






template<unsigned int dim>
FE_CR<dim>::FE_CR()
: FiniteElement<dim>()
{
    this->function_space_ = std::make_shared<PolynomialSpace>(1,dim);
    
    if (dim == 0)
    {
        this->dofs_.push_back(Dof(0, 0, { 1 }, { 1 }, Value));
    }
    else
    {
        arma::vec::fixed<dim> sp; // support point
        for (unsigned int sid=0; sid<RefElement<dim>::n_sides; ++sid)
        {
            sp.fill(0);
            for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; ++i)
                sp += RefElement<dim>::node_coords(RefElement<dim>::interact(Interaction<0,dim-1>(sid))[i]);
            sp /= RefElement<dim>::n_nodes_per_side;
            // barycentric coordinates
            arma::vec::fixed<dim+1> bsp;
            bsp.subvec(1,dim) = sp;
            bsp[0] = 1. - arma::sum(sp);
            this->dofs_.push_back(Dof(dim-1, sid, bsp, { 1 }, Value));
        }
    }

    this->setup_components();
    this->compute_node_matrix();
}




template<unsigned int dim>
FE_CR_disc<dim>::FE_CR_disc()
: FiniteElement<dim>()
{
    this->function_space_ = std::make_shared<PolynomialSpace>(1,dim);

    if (dim == 0)
    {
        this->dofs_.push_back(Dof(0, 0, { 1 }, { 1 }, Value));
    }
    else
    {
        arma::vec::fixed<dim> sp; // support point
        for (unsigned int sid=0; sid<RefElement<dim>::n_sides; ++sid)
        {
            sp.fill(0);
            for (unsigned int i=0; i<RefElement<dim>::n_nodes_per_side; ++i)
                sp += RefElement<dim>::node_coords(RefElement<dim>::interact(Interaction<0,dim-1>(sid))[i]);
            sp /= RefElement<dim>::n_nodes_per_side;
            // barycentric coordinates
            arma::vec::fixed<dim+1> bsp;
            bsp.subvec(1,dim) = sp;
            bsp[0] = 1. - arma::sum(sp);
            this->dofs_.push_back(Dof(dim, 0, bsp, { 1 }, Value));
        }
    }

    this->setup_components();
    this->compute_node_matrix();
}





template class FE_P<0>;
template class FE_P<1>;
template class FE_P<2>;
template class FE_P<3>;


template class FE_P_disc<0>;
template class FE_P_disc<1>;
template class FE_P_disc<2>;
template class FE_P_disc<3>;


template class FE_CR<0>;
template class FE_CR<1>;
template class FE_CR<2>;
template class FE_CR<3>;


template class FE_CR_disc<0>;
template class FE_CR_disc<1>;
template class FE_CR_disc<2>;
template class FE_CR_disc<3>;

