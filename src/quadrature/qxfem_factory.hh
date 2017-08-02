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
class QXFEMFactory {
public:
    typedef typename Space<3>::Point Point;
    typedef typename std::shared_ptr<Singularity0D> Singularity0DPtr;
    typedef typename std::shared_ptr<Singularity1D> Singularity1DPtr;
    /**
     * @brief Create a formula of given order.
     *
     * The formula is exact for polynomials of degree @p order.
     */
    QXFEMFactory(unsigned int max_level = 10)
    : max_level_(max_level), level_offset_(0){}
    
    
    std::shared_ptr<QXFEM<2,3>> create_singular(const std::vector<Singularity0DPtr> & sing,
                                                         ElementFullIter ele);
    
    std::shared_ptr<QXFEM<3,3>> create_singular(const std::vector<Singularity1DPtr> & sing,
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
    template<int dim>
    void gnuplot_refinement(ElementFullIter ele,
                            const string& output_dir,
                            const QXFEM<dim,3>& quad,
                            const std::vector<Singularity0DPtr> & sing);
    
    template<int dim>
    void gnuplot_refinement(ElementFullIter ele,
                            const string& output_dir,
                            const QXFEM<dim,3>& quad);
    
    
protected:
    
    struct AuxSimplex{
        AuxSimplex() {
            active = true; 
            refine = false;
            sing_id = -1;
        }
        std::vector<Point> nodes;
        int sing_id;
        bool active;
        bool refine;
        template<int dim> double measure() const;
        template<int dim> double compute_max_h() const;
    };

    
    /** @brief Marks simplices which are to be refined.
     * Criterions:
     * 1] refine simplex which intersects with the edge of the singularity
     */
    unsigned int mark_refinement_2D(const std::vector<Singularity0DPtr> & sing);
    
    unsigned int mark_refinement_3D(const std::vector<Singularity1DPtr> & sing);
    
    /// Refines marked simplices.
    template<int dim>
    void refine_level(unsigned int n_simplices_to_refine);
    /// Refines single simplex.
    template<int dim>
    void refine_simplex(const AuxSimplex &aux_element);
    
    template<int dim, class S>
    void distribute_qpoints(std::vector< Point >& real_points, 
                            std::vector< double >& weights,
                            const std::vector<S> & sing,
                            double ele_measure);
    
    template<int dim, class S>
    bool point_in_singularity(const Point& p, S s);
    
    template<int dim>
    void map_real_to_unit_points(const std::vector<Point>& real_points,
                                 std::vector< typename Space< dim >::Point >& unit_points,
                                 ElementFullIter ele);
    
    /** @brief Tests intersection of 2d simplex and singularity (2d circle) in spacedim-dimensional (2 or 3) space.
     * @return intersection type:
     * 
     */
    int sigularity0D_distance(const Singularity0D& w,
                              const AuxSimplex &s,
                              double& distance);
    
    /** @brief Tests intersection of 3d simplex and singularity (cylinder) in 3d space.
     * Uses bounding boxes as an approximation.
     * @return intersection type:
     * -1:  tetrahedron of @p max_h at @p distance from singularity
     * 0:   tetrahedron inside singularity (all nodes inside)
     * 1:   one or more nodes of tetrahedron inside singularity
     * 2:   singularity is inside tetrahedron
     * 3:   distance is less then radius
     */
    int singularity1D_distance(const Singularity1D& w,
                              const AuxSimplex &s,
                              double& distance);
    
    /// Level of current refinement.
    unsigned int level_;
    unsigned int max_level_;
    double max_h_level_0_;
    
    std::vector<AuxSimplex> simplices_;
    unsigned int level_offset_;
    
    // Refinement criteria constant on elements without well inside. [2D]
    static const double distance_criteria_factor_2d_;
    // Refinement criteria constant on elements without well inside. [3D]
    static const double distance_criteria_factor_3d_;
};



#endif // QXFEM_FACTORY_HH_