/*
 * field_add_potential.hh
 *
 *  Created on: Jan 19, 2013
 *      Author: jb
 */

#ifndef FIELD_ADD_POTENTIAL_HH_
#define FIELD_ADD_POTENTIAL_HH_

#include <armadillo>

/**
 * This field is meant to be used to implement two possibilities for initialization of pressure fields in
 * Darcy flows. You can either use pressure of piezo-metric  head which are  related by adding gravity potential.
 * For various reasons we use pressure as the primary variable, so if the user enters piezo-head we need to add the potential to
 * the field he/she has provided. This is done by this class. Unfortunately it introduce one more level of indirection,
 * namely one more virtual call for getting the field value.
 *
 * - The field can not be initialized form the input.
 * - We allow only Scalar Value with element_type double.
 */
template <int spacedim, class Value>
class FieldAddPotential : public FieldBase<spacedim, Value> {
public:

    /**
     *
     */
    FieldAddPotential( const arma::vec::fixed<spacedim+1> &potential_grad, const Input::AbstractRecord &rec, unsigned int n_comp=0);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point<spacedim> &p, ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point<spacedim> >  &point_list, ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldAddPotential();

private:
    /// Field to which we add linear potential.
    FieldBase<spacedim, Value> *inner_field_;
    /// Potential gradient.
    arma::vec::fixed<spacedim> grad_;
    /// Potential constant term.
    double zero_level_;
};


#endif /* FIELD_ADD_POTENTIAL_HH_ */
