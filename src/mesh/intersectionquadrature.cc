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
 * @file    intersection.cc
 * @brief   
 */

#include "mesh/intersection.hh"
#include "mesh/mesh.h"

#include <boost/tokenizer.hpp>
#include "boost/lexical_cast.hpp"
#include <armadillo>

// inicializovat objekt, cist zbytek tokenu z tok a naplnit map a shift pro master a slave
// viz dokumentace k Armadillu
Intersection::Intersection(const  ElementFullIter ele_master, const ElementFullIter ele_slave,
    			 const IntersectionLocal *isec)
:	dim(isec->n_points() - 1),
	master(ele_master), slave(ele_slave),
	master_map(master->dim(), dim), slave_map(slave->dim(), dim),
	master_shift(master->dim()), slave_shift(slave->dim())
{
	///otestuje se jestli dimenze masteru je mensi nez dimenze slave - chybova hlaska (vyjimka - throw)
	///pocet pointu=dim+1
	if (master->dim() > slave->dim()) {
		cout << "Exception: master->dim() > slave->dim()" << endl;
		//throw((char*) "master->dim > slave->dim");
	}

	intersection_point_to_vectors(isec->get_point(0),master_shift, slave_shift);

	arma::vec master_tmp(master_shift), slave_tmp(slave_shift);
	// cyklus pres body pruniku
	for (unsigned int i = 1; i < (dim + 1); ++i) {
		intersection_point_to_vectors(isec->get_point(i),master_tmp, slave_tmp);
		master_tmp -= master_shift;
		slave_tmp -= slave_shift;

		master_map.col(i-1) = master_tmp;
		slave_map.col(i-1) = slave_tmp;
	}
}



unsigned int Intersection::master_dim()
    {return master->dim();}



unsigned int Intersection::slave_dim()
    {return slave->dim();}



void Intersection::intersection_point_to_vectors(const IntersectionPoint *point, arma::vec &vec1, arma::vec &vec2)
{
	const vector<double> &coord_el1 = point->el1_coord();
	OLD_ASSERT_EQUAL(coord_el1.size() , vec1.n_elem);
	vec1=arma::vec(coord_el1);

	const vector<double> &coord_el2 = point->el2_coord();
	OLD_ASSERT_EQUAL(coord_el2.size() , vec2.n_elem);
	vec2=arma::vec(coord_el2);
}


arma::vec Intersection::map_to_master(const arma::vec &point) const
{
	//dim = dimenze intersec elementu
	OLD_ASSERT(( point.n_elem == dim ),"Map to slave: point.n_elem(%d) != dim(%d) \n", point.n_elem, dim);
    int result_dim = master->dim();
    arma::vec result(result_dim+1);
	result(0)=1.0;
	result.subvec(1, result_dim) = (master_map * point + master_shift);
	return result;
}

arma::vec Intersection::map_to_slave(const arma::vec &point) const
{
	OLD_ASSERT(( point.n_elem == dim ),"Map to slave: point.n_elem(%d) != dim(%d) \n", point.n_elem, dim);
	int result_dim = slave->dim();
	arma::vec result(result_dim+1);
	result(0)=1.0;
	result.subvec(1, result_dim) = (slave_map * point + slave_shift);
	return result;
}

double Intersection::intersection_true_size() const {

    static const double factorial[4] = {1.0, 1.0, 2.0, 6.0};
	return (master->measure() * det(master_map) / factorial[dim]);
}


