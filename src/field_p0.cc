/*
 * field_p0.cc
 *
 *  Created on: Sep 12, 2011
 *      Author: jb
 */

#include "field_p0.hh"
#include "system/parsed_function.hh"
#include "mesh/mesh.h"

// specialization
template <>
void FieldP0<double>::setup_from_function(const std::string &expr)
{
    ParsedFunction func;
    func.set_expression(expr);
    func.set_time(0.0);

    data.resize(mesh->element.size());
    FOR_ELEMENTS(mesh,ele) {
        arma::vec3 center(ele->centre());
        data[ele.index()]=func.value( center ); // has to make explicit local copy
    }
}
