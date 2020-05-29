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
 * @file    surface_depth.cc
 * @brief
 */

#include <vector>
#include <armadillo>
#include <limits>

#include "fields/surface_depth.hh"
#include "mesh/accessors.hh"
#include "mesh/mesh.h"
#include "mesh/bc_mesh.hh"
#include "mesh/ref_element.hh"
#include "system/sys_profiler.hh"
#include "mesh/accessors.hh"
#include "mesh/range_wrapper.hh"


SurfaceDepth::SurfaceDepth(const Mesh *mesh, std::string surface_region, std::string surface_direction)
: surface_region_(surface_region) {
    this->create_projection_matrix( arma::vec3(surface_direction) );
    this->construct_bih_tree( const_cast<Mesh*>(mesh), surface_region);
    this->projection_search_radius_ = mesh->global_snap_radius();
    this->searched_elements_.reserve(BIHTree::default_leaf_size_limit);
}


void SurfaceDepth::create_projection_matrix(arma::vec3 surface_vec)
{
    if (surface_vec(2)==0) {
        THROW( ExcSurfaceProjection() << EI_Message("Vertical plane of surface is forbidden!") );
    }
    surface_norm_vec_ = surface_vec / arma::norm(surface_vec, 2); // normalize to size == 1

    arma::vec3 ex("1 0 0"), ey("0 1 0");
    arma::vec3 ta = ex - arma::dot(ex, surface_norm_vec_) * surface_norm_vec_;
    ta /= arma::norm(ta, 2); // normalize
    arma::vec3 tb = ey - arma::dot(ey, surface_norm_vec_) * surface_norm_vec_;
    tb /= arma::norm(tb, 2); // normalize
    m_ = arma::mat(2,3);
    m_.row(0) = ta.t();
    m_.row(1) = tb.t();
}


void SurfaceDepth::construct_bih_tree(Mesh *mesh, std::string surface_region)
{
	std::vector<BoundingBox> boxes;
	inv_projection_.clear();
	b_vecs_.clear();

	RegionSet region_set = mesh->region_db().get_region_set(surface_region);
    if (region_set.size() == 0)
        THROW( RegionDB::ExcUnknownSet() << RegionDB::EI_Label(surface_region) );
    for (auto reg : region_set) {
    	if (reg.dim() != 2)
    		THROW( ExcSurfaceProjection()
    		    			<< EI_Message("Dimension of surface region " + surface_region + " must be 2!") );
    }

    // make element boxes
    unsigned int i=0;
    unsigned int i_node;
    arma::vec3 project_node("0 0 0");
    for( auto ele : mesh->get_bc_mesh()->elements_range() ) {
        if (ele.region().is_in_region_set(region_set)) {
        	ASSERT_EQ(ele->n_nodes(), 3);

        	arma::vec projection = m_ * (*ele.node(0));
        	project_node(0) = projection(0); project_node(1) = projection(1);
            BoundingBox bb(project_node);
            for(i_node=1; i_node<ele->n_nodes(); i_node++) {
                arma::vec project_coords = m_ * (*ele.node(i_node));
                project_node(0) = project_coords(0); project_node(1) = project_coords(1);
                bb.expand(project_node);
            }
            boxes.push_back(bb);

            arma::mat a_mat(3,3);
        	a_mat.col(0) = *ele.node(1) - *ele.node(0);
        	a_mat.col(1) = *ele.node(2) - *ele.node(0);
        	a_mat.col(2) = surface_norm_vec_;

            inv_projection_.push_back( a_mat.i() );
            b_vecs_.push_back( *ele.node(0) );
        }
        i++;
    }

    if ( boxes.size() == 0) {
    	THROW( ExcSurfaceProjection()
    			<< EI_Message("Region " + surface_region + " contains no boundary element! Probably bulk region was set.") );
    }

    bih_tree_.add_boxes( boxes );
    bih_tree_.construct();
}

double SurfaceDepth::compute_distance(arma::vec3 point)
{
	double distance = std::numeric_limits<double>::max();
	bool found_surface_projection = false;
	arma::vec3 project_point;
	arma::vec3 x;

	project_point.subvec(0,1) = m_ * point;
	project_point(2) = 0;

	searched_elements_.clear();
	bih_tree_.find_point(project_point, searched_elements_);
	for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++) {
		prepare_distance_solve( (*it), point, x );
		if ( (x(0)>=0) && (x(1)>=0) && (x(0)+x(1)<=1) ) { // point is in triangle
			double new_distance = -x(2);
			if ( (new_distance>=0) && (new_distance<distance) ) distance = new_distance;
			found_surface_projection = true;
		}
	}

	if (!found_surface_projection) {
		double snap_dist = std::numeric_limits<double>::max();
		arma::vec3 proj_to_surface_plane;

		for (std::vector<unsigned int>::iterator it = searched_elements_.begin(); it!=searched_elements_.end(); it++) {
			// check snap distance point - triangle
			prepare_distance_solve( (*it), point, x );

			proj_to_surface_plane = point - x(2)*surface_norm_vec_;
			arma::vec local_point = x.subvec(0,1);
			auto bary_point = RefElement<2>::local_to_bary(local_point);
			auto clip_point = RefElement<2>::clip(bary_point);
			auto proj_ref = RefElement<2>::bary_to_local(clip_point);
			auto proj_3d = inv_projection_[*it].submat(0,0,2,1) * proj_ref;
			double new_snap_dist = arma::norm(proj_3d - proj_to_surface_plane, 2);
			if (new_snap_dist<snap_dist) {
				snap_dist = new_snap_dist;
				distance = -x(2);
			}
		}

		if (snap_dist > projection_search_radius_) {
	        THROW(ExcTooLargeSnapDistance() << EI_Xcoord(proj_to_surface_plane(0)) << EI_Ycoord(proj_to_surface_plane(1))
	        		<< EI_Zcoord(proj_to_surface_plane(2)) << EI_SnapDistance(snap_dist) << EI_RegionName(surface_region_) );
		}
	}

	return distance;
}

void SurfaceDepth::prepare_distance_solve(unsigned int elem_idx, arma::vec3 &point, arma::vec3 &x)
{
	x = inv_projection_[elem_idx] * (point - b_vecs_[elem_idx]);
}
