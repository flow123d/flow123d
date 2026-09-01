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
 * @file    eval_points.cc
 * @brief
 * @author  David Flanderka
 */

#include "fem/eval_points.hh"
#include "fem/integral_acc.hh"
#include "fem/patch_fe_values.hh"
#include "quadrature/quadrature.hh"
#include "mesh/ref_element.hh"
#include <memory>


const unsigned int EvalPoints::undefined_dim = 10;

EvalPoints::EvalPoints()
: dim_eval_points_({DimEvalPoints(0), DimEvalPoints(1), DimEvalPoints(2), DimEvalPoints(3)}),
  max_size_(0)
{}

template <unsigned int dim>
std::shared_ptr<internal_integrals::Bulk> EvalPoints::add_bulk_internal(Quadrature *quad) {
    ASSERT_EQ(dim, quad->dim());

    auto tpl = IntegralTplHash::integral_tuple(dim, quad->size());
    auto map_it = bulk_integrals_.find( tpl );
    if (map_it == bulk_integrals_.end()) {
        dim_eval_points_[dim].add_local_points<dim>( quad->get_points() );
        uint i_subset = dim_eval_points_[dim].add_subset(points_domain::bulk_points, quad->size());

        bulk_integrals_[tpl] = std::make_shared<internal_integrals::Bulk>(quad, quad->dim(), shared_from_this(), i_subset);
        map_it = bulk_integrals_.find( tpl );
        this->set_max_size();
    }
    return map_it->second;
}

template <>
std::shared_ptr<internal_integrals::Bulk> EvalPoints::add_bulk_internal<0>(Quadrature *quad)
{
    ASSERT_EQ(0, quad->dim());

    auto tpl = IntegralTplHash::integral_tuple(0, quad->size());
    auto map_it = bulk_integrals_.find( tpl );
    if (map_it == bulk_integrals_.end()) {
        uint i_subset = dim_eval_points_[0].add_subset(points_domain::bulk_points, quad->size());

        bulk_integrals_[tpl] = std::make_shared<internal_integrals::Bulk>(quad, quad->dim(), shared_from_this(), i_subset);
        map_it = bulk_integrals_.find( tpl );
        this->set_max_size();
    }
    return map_it->second;
}

template <unsigned int dim>
std::shared_ptr<internal_integrals::Edge> EvalPoints::add_edge_internal(Quadrature *quad)
{
    ASSERT_EQ(dim, quad->dim()+1);

    auto tpl = IntegralTplHash::integral_tuple(dim, quad->size());
    auto map_it = edge_integrals_.find( tpl );
    if (map_it == edge_integrals_.end()) {
        for (unsigned int i=0; i<dim+1; ++i) {  // sides
            Quadrature high_dim_q = quad->make_from_side<dim>(i);
            dim_eval_points_[dim].add_local_points<dim>( high_dim_q.get_points() );
        }
        uint i_subset = dim_eval_points_[dim].add_subset(points_domain::side_points, quad->size(), quad->size() * dim);

        edge_integrals_[tpl] = std::make_shared<internal_integrals::Edge>(quad, quad->dim()+1, shared_from_this(), i_subset);
        map_it = edge_integrals_.find( tpl );
        this->set_max_size();
    }
    return map_it->second;
}

uint EvalPoints::get_max_bulk_quad_size(unsigned int dim) const {
    return get_max_integral_quad_size<internal_integrals::Bulk>(bulk_integrals_, dim);
}

uint EvalPoints::get_max_side_quad_size(unsigned int dim) const {
    return get_max_integral_quad_size<internal_integrals::Edge>(edge_integrals_, dim);
}

template<class Integral>
uint EvalPoints::get_max_integral_quad_size(IntegralPtrMap<Integral> integrals, unsigned int dim) const {
    uint max_qsize=0;
    for (auto integral_it : integrals)
        if (integral_it.second->dim() == dim)
            if (integral_it.second->quad()->size() > max_qsize)
                max_qsize = integral_it.second->quad()->size();
    return max_qsize;
}

std::vector<Quadrature *> EvalPoints::get_bulk_quad_vector(unsigned int dim) const {
    return get_quad_vector<internal_integrals::Bulk>(bulk_integrals_, dim);
}

std::vector<Quadrature *> EvalPoints::get_side_quad_vector(unsigned int dim) const {
    return get_quad_vector<internal_integrals::Edge>(edge_integrals_, dim);
}

template<class Integral>
std::vector<Quadrature *> EvalPoints::get_quad_vector(IntegralPtrMap<Integral> integrals, unsigned int dim) const {
    std::vector<Quadrature *> quad_vec;
    for (auto integral_it : integrals)
        if (integral_it.second->dim() == dim)
            quad_vec.push_back( integral_it.second->quad() );
    return quad_vec;
}

EvalPoints::DimEvalPoints::DimEvalPoints(unsigned int dim)
: local_points_(dim), n_subsets_(0), dim_(dim)
{
	subset_starts_[0] = 0;
	if (dim>0) local_points_.reinit(EvalPoints::max_subsets*EvalPoints::max_subset_points);
}

template <unsigned int dim>
void EvalPoints::DimEvalPoints::add_local_points(const Armor::Array<double> & quad_points) {
    ASSERT_GT(dim, 0).error("Dimension 0 not supported!\n");
    unsigned int local_points_old_size = local_points_.size();
    local_points_.resize(quad_points.size() + local_points_old_size);
    for (unsigned int i=0; i<quad_points.size(); ++i) {
        //DebugOut() << "add i: " << i << " : " <<  quad_points.vec<dim>(i);

        local_points_.set(i+local_points_old_size) = quad_points.vec<dim>(i);
	}
}


uint EvalPoints::DimEvalPoints::add_subset(points_domain point_domain, unsigned int quad_size, unsigned int repeated_points) {
    ASSERT_LT(n_subsets_, EvalPoints::max_subsets).error("Maximal number of subsets exceeded!\n");

    n_subsets_++;
    subset_starts_[n_subsets_] = this->size();
    for (uint i=0; i<quad_size; ++i)
        points_domains_.push_back(point_domain);
    for (uint i=0; i<repeated_points; ++i)
        points_domains_.push_back(points_domain::repeated_side_points);
    return n_subsets_ - 1;
}


template std::shared_ptr<internal_integrals::Bulk> EvalPoints::add_bulk_internal<0>(Quadrature *);
template std::shared_ptr<internal_integrals::Bulk> EvalPoints::add_bulk_internal<1>(Quadrature *);
template std::shared_ptr<internal_integrals::Bulk> EvalPoints::add_bulk_internal<2>(Quadrature *);
template std::shared_ptr<internal_integrals::Bulk> EvalPoints::add_bulk_internal<3>(Quadrature *);
template std::shared_ptr<internal_integrals::Edge> EvalPoints::add_edge_internal<1>(Quadrature *);
template std::shared_ptr<internal_integrals::Edge> EvalPoints::add_edge_internal<2>(Quadrature *);
template std::shared_ptr<internal_integrals::Edge> EvalPoints::add_edge_internal<3>(Quadrature *);
template unsigned int EvalPoints::get_max_integral_quad_size<internal_integrals::Bulk>(IntegralPtrMap<internal_integrals::Bulk>, unsigned int) const;
template unsigned int EvalPoints::get_max_integral_quad_size<internal_integrals::Edge>(IntegralPtrMap<internal_integrals::Edge>, unsigned int) const;

template void EvalPoints::DimEvalPoints::add_local_points<1>(const Armor::Array<double> &);
template void EvalPoints::DimEvalPoints::add_local_points<2>(const Armor::Array<double> &);
template void EvalPoints::DimEvalPoints::add_local_points<3>(const Armor::Array<double> &);
