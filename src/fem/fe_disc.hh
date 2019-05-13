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
 * @file    fe_disc.hh
 * @brief   Auxiliary templates for FE with discontinuous dofs.
 * @author  Jan Stebel
 */

#ifndef FE_DISC_HH_
#define FE_DISC_HH_

#include "fem/finite_element.hh"





/**
 * @brief Discontinuous finite element.
 *
 * The dofs are not shared across faces/edges/vertices.
 * 
 * @param FE   Base finite element class.
 * @param Args Optional argument types for constructor of FE.
 * 
 * Example use:
 * 
 *   // polynomial finite element in 3D
 *   typedef FE_disc<FE_P<3>, unsigned int> FE_P_3D_disc;
 *   FE_P_3D_disc fe_p(2);
 * 
 *   // lowest order Raviart-Thomas element in 2D
 *   typedef FE_disc<FE_RT0<2> > FE_RT0_2D_disc;
 *   FE_RT0_2D_disc fe_rt;
 * 
 */
template <class FE, typename ... Args>
class FE_disc : public FE
{
public:

    FE_disc(Args... args) : FE(args ...)
    {
        unsigned int dim = this->function_space_->space_dim();
        // all dofs are "inside" the cell (not shared with neighbours)
        for (unsigned int i=0; i<this->dofs_.size(); i++)
            this->dofs_[i].dim = dim;
        
        this->setup_components();

        this->compute_node_matrix();
    }

};













#endif /* FE_DISC_HH_ */
