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
 * @file    lazy_dependecy.hh
 * @brief   
 */

#ifndef LAZY_DEPENDECY_HH_
#define LAZY_DEPENDECY_HH_

#include <boost/lambda/lambda.hpp>
#include <list>
#include <algorithm>

using namespace boost;

/**
 * This class implements lazy dependency of objects that inherits from it.
 *
 * Consider you have equation E which depends on input vector fields A and B which are possibly results of another computations
 * and can change over time. When someone ask you for result of equation E you have to look if your current result is up to date according
 * to input fields A and B. This is purpose of this class.
 *
 * How to use it:
 * 1) make all classes in the dependency graph inherited form LazyDependency class. Multiple inheritance is OK in this case.
 * 2) equation E should call add_dependency for A and B
 * 3) every time you change A or B call update()
 * 4) when you are asked for result of equation check needs_update() and update the result only if it returns true
 *    Possibly call also update() for possible objects that depends on E.
 */
class LazyDependency {
public:
    /// Default constructor.
    LazyDependency()
    : change_set_(1)
    {}

    /// Increase the change set and set actual values of change sets of objects we depend on.
    void update() {
        change_set_ ++;
        std::for_each(dependencies_.begin(), dependencies_.end(),
                _1.second() = _1.first().change_set_
                );
    }

    /**
     * Adds new object into dependency list. This do not update this dependency!
     * @param object - any instance of LazyDependency
     */
    void add_dependency(LazyDependency &object) {
        dependencies_.push_back(std::pair<LazyDependency&, unsigned int>(object, 0));
    }

    /**
     * Returns true if there is at least one object in the dependency list that changed its change set
     * since the last call of update().
     */
    bool needs_update() const {
        bool no_chnage=true;
        std::for_each(dependencies_.begin(), dependencies_.end(),
                no_change = no_change && ( _1.first().change_set_ == _1.second() )
                );
        return ! no_change;
    }

private:
    unsigned int change_set_;
    std::list<std::pair<LazyDependency &, unsigned int> > dependencies_;
};

#endif /* LAZY_DEPENDECY_HH_ */
