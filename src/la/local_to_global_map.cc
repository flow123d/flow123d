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
 * @file    local_to_global_map.cc
 * @brief   
 */

#include <la/local_to_global_map.hh>
#include <la/distribution.hh>

#include "system/system.hh"


LocalToGlobalMap::LocalToGlobalMap(const Distribution &any_distr)
: nonlocal_indices_(new std::set<unsigned int>()),
  distr_(boost::make_shared<Distribution>(any_distr)),
  global_indices_(0)
{}

LocalToGlobalMap::LocalToGlobalMap(boost::shared_ptr<Distribution> any_distr)
: nonlocal_indices_(new std::set<unsigned int>()),
  distr_(any_distr),
  global_indices_(0)
{}

void LocalToGlobalMap::insert(const unsigned int global_idx) {
	OLD_ASSERT( global_indices_.size() == 0, "Insertion into the map after finalize.");
    if (! distr_->is_local( global_idx ) ) nonlocal_indices_->insert( global_idx );
}

void LocalToGlobalMap::insert(const std::vector<unsigned int> &indices) {
	OLD_ASSERT( global_indices_.size() == 0, "Insertion into the map after finalize.");
    for(std::vector<unsigned int>::const_iterator it=indices.begin(); it != indices.end(); ++it)
        if (! distr_->is_local( *it ) ) nonlocal_indices_->insert( *it );
}

void LocalToGlobalMap::finalize() {
    global_indices_.resize( distr_->lsize() + nonlocal_indices_->size() );
    unsigned int i;

    for(i=0; i< distr_->lsize(); i++) global_indices_[i]=distr_->begin()+i;
    for(std::set<unsigned int>::iterator it = nonlocal_indices_->begin();
        it != nonlocal_indices_->end(); ++it)
        global_indices_[i++]=*it;
    delete nonlocal_indices_;
}



