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
 * @file    qxfem_factory.hh
 * @brief   XFEM adaptive quadrature factory.
 * @author  Pavel Exner
 */

#ifndef QXFEM_FACTORY_HH_
#define QXFEM_FACTORY_HH_

#include "system/global_defs.h"

#include "quadrature/qxfem.hh"

#include "fem/singularity.hh"

#include "mesh/point.hh"
#include "mesh/mesh_types.hh"

/**
 * @brief XFEM adaptive quadrature factory.
 * 
 * 1) refine simplex in real coordinates into auxsimplices
 *      - single vector of auxsimplices, dont delete old ones, just mark non-active
 *      - got through active on current level, set flag for refinement according to criteria
 *      - refine auxsimplices
 * 2) go through active auxsimplices:
 *      - compute measure; divide by measure of element -> weight of auxsimplex
 *      - select QGauss, map quarature points onto auxsimplex -> real quad points
 *      - quad points weight = QGauss wheight * weight of auxsimplex
 * 3) go through all quadrature points a map them to RefElement
 * 4) possibly gnuplot, separately active auxsimplices, then real quad points
 * 5) delete auxiliary objects and return pure quadrature
 */
template<int dim, int spacedim>
class QXFEMFactory {
public:
    typedef typename Space<spacedim>::Point Point;
    /**
     * @brief Create a formula of given order.
     *
     * The formula is exact for polynomials of degree @p order.
     */
    QXFEMFactory(unsigned int max_level = 10)
    : max_level_(max_level), level_offset_(0){}
    
    
    std::shared_ptr<QXFEM<dim,spacedim>> create_singular(const std::vector<Singularity0D<spacedim>> & sing,
                                                         ElementFullIter ele);
    
    /// Clears temporary refinement simplices, resets the factory to initial state.
    void clear();
    
    /** @brief Calls gnuplot to create image of refined element.
     * 
      * Also can save the gnuplot script to file.
      * @param output_dir is the directory for output_dir
      * @param real is true then the element is printed in real coordinates
      * @param show is true then the gnuplot utility is started and plots the refinement on the screen
      */ 
    void gnuplot_refinement(ElementFullIter ele,
                            const string& output_dir,
                            const QXFEM<dim,spacedim>& quad,
                            const std::vector<Singularity0D<spacedim>> & sing);
    
    
private:
    
    struct AuxSimplex{
        std::vector<Point> nodes;
        int sing_id;
        bool active;
        bool refine;
    };

    
    /** @brief Marks simplices which are to be refined.
     * Criterions:
     * 1] refine simplex which intersects with the edge of the singularity
     */
    unsigned int refine_edge(const std::vector<Singularity0D<spacedim>> & sing);
    
    /// Refines marked simplices.
    void refine_level(unsigned int n_simplices_to_refine);
    /// Refines single simplex.
    void refine_simplex(AuxSimplex &aux_element);
    
    void distribute_qpoints(std::vector< Point >& real_points, 
                            std::vector< double >& weights,
                            const std::vector<Singularity0D<spacedim>> & sing);
    
    void map_real_to_unit_points(const std::vector<Point>& real_points,
                                 std::vector< typename Space< dim >::Point >& unit_points,
                                 ElementFullIter ele);
    
    /** @brief Tests intersection of 2d simplex and singularity (2d circle) in spacedim-dimensional (2 or 3) space.
     * @return intersection type:
     * 
     */
    int simplex_sigularity_intersection(const Singularity0D<spacedim>& w, AuxSimplex &s);
    
    /// Level of current refinement.
    unsigned int level_;
    unsigned int max_level_;
    
    std::vector<AuxSimplex> simplices_;
    unsigned int level_offset_;
    
//     // Square refinement criteria constant on the cells without well inside.
//     static const double square_refinement_criteria_factor_;
};


#endif // QXFEM_FACTORY_HH_