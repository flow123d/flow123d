/*
 * field_values.cc
 *
 *  Created on: Dec 11, 2012
 *      Author: jb
 */


#include "fields/field_values.hh"

namespace internal {
// Helper functions to get scalar type name
std::string type_name_(double) { return "Real"; }
std::string type_name_(int)  { return "Int"; }
std::string type_name_(std::string) { return "String"; }
std::string type_name_(FieldEnum) { return "Enum"; }


} // namespace internal


template class FieldValue_<1,1,FieldEnum>;
template class FieldValue_<1,1,int>;
template class FieldValue_<0,1,int>;
template class FieldValue_<1,1,double>;
template class FieldValue_<0,1,double>;

template class FieldValue_<2,1,double>;
template class FieldValue_<3,1,double>;

template class FieldValue_<2,2,double>;
template class FieldValue_<3,3,double>;

template class FieldValue_<1,1,std::string>;
template class FieldValue_<0,1,std::string>;

template class FieldValue_<2,1,std::string>;
template class FieldValue_<3,1,std::string>;

template class FieldValue_<2,2,std::string>;
template class FieldValue_<3,3,std::string>;

template class FieldValue<1>;
template class FieldValue<2>;
template class FieldValue<3>;
