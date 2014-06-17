/*
 * field_add_potential.impl.hh
 *
 *  Created on: Jan 22, 2013
 *      Author: jb
 */

#ifndef FIELD_ADD_POTENTIAL_IMPL_HH_
#define FIELD_ADD_POTENTIAL_IMPL_HH_

#include "fields/field_add_potential.hh"


template <int spacedim, class Value>
FieldAddPotential<spacedim, Value>::FieldAddPotential(const arma::vec::fixed<spacedim+1> &potential, const Input::AbstractRecord &rec, unsigned int n_comp)
: FieldAlgorithmBase<spacedim, Value>(n_comp),
  inner_field_( FieldAlgorithmBase<spacedim, Value>::function_factory(rec, this->value_.n_rows()) )
{
    grad_=potential.subvec(0,spacedim-1);
    zero_level_=potential[spacedim];
}



/**
 * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
 */
template <int spacedim, class Value>
typename Value::return_type const & FieldAddPotential<spacedim, Value>::value(const Point &p, const ElementAccessor<spacedim> &elm)
{
    this->r_value_ = inner_field_->value(p,elm);

    double potential=arma::dot(grad_ , p) + zero_level_;
    for(unsigned int row=0; row < this->value_.n_rows(); row++)
        for(unsigned int col=0; col < this->value_.n_cols(); col++)
            this->value_(row,col) += potential;

    return this->r_value_;
}



/**
 * Returns std::vector of scalar values in several points at once.
 */
template <int spacedim, class Value>
void FieldAddPotential<spacedim, Value>::value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list)
{
    ASSERT_EQUAL( point_list.size(), value_list.size() );
    inner_field_->value_list(point_list, elm, value_list);
    for(unsigned int i=0; i< point_list.size(); i++) {
        double potential= arma::dot(grad_ , point_list[i]) + zero_level_;
        Value envelope(value_list[i]);

        for(unsigned int row=0; row < this->value_.n_rows(); row++)
            for(unsigned int col=0; col < this->value_.n_cols(); col++)
                envelope(row,col) += potential;
    }
}

template <int spacedim, class Value>
bool FieldAddPotential<spacedim,Value>::set_time (double time)
{
    ASSERT(inner_field_, "Null data pointer.\n");
    //DBGMSG("set_time for field_add_potential\n");
    return inner_field_->set_time(time);
}


template <int spacedim, class Value>
FieldAddPotential<spacedim, Value>::~FieldAddPotential() {
}


#endif /* FIELD_ADD_POTENTIAL_IMPL_HH_ */
