/*
 * field_add_potential.hh
 *
 *  Created on: Jan 19, 2013
 *      Author: jb
 */

#ifndef FIELD_ADD_POTENTIAL_HH_
#define FIELD_ADD_POTENTIAL_HH_

#include <armadillo>
#include <memory>

#include "fields/field.hh"
#include "flow/old_bcd.hh"


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
class FieldAddPotential : public FieldAlgorithmBase<spacedim, Value> {
public:
	typedef typename Field<spacedim, Value>::FactoryBase FactoryBaseType;
    typedef typename Space<spacedim>::Point Point;
    /**
     *
     */
    FieldAddPotential( const arma::vec::fixed<spacedim+1> &potential_grad, const Input::AbstractRecord &rec, unsigned int n_comp=0);


    /**
     * Factory class (descendant of @p Field<...>::FactoryBase) that is necessary
     * for setting pressure values are piezometric head values.
     */
    class FieldFactory : public FactoryBaseType {
    public:
    	/// Constructor.
    	FieldFactory(arma::vec::fixed<spacedim+1> potential, std::string field_name)
    	: potential_(potential),
    	  field_name_(field_name)
    	{}

    	virtual typename Field<spacedim,Value>::FieldBasePtr create_field(Input::Record rec, const FieldCommon &field) {
       		OldBcdInput *old_bcd = OldBcdInput::instance();
       		old_bcd->read_flow_record(rec, field);
       		auto field_ptr = old_bcd->flow_pressure;

       		Input::AbstractRecord field_a_rec;
        	if (! field_ptr && rec.opt_val(field_name_, field_a_rec)) {
        		return std::make_shared< FieldAddPotential<3, FieldValue<3>::Scalar > >( potential_, field_a_rec);
        	} else {
        		return field_ptr;
        	}
    	}

    	arma::vec::fixed<spacedim+1> potential_;
    	std::string field_name_;
    };


    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    /**
     * Update time and possibly update data.
     */
    bool set_time(const TimeStep &time) override;
    
    
    virtual ~FieldAddPotential();

private:
    /// Field to which we add linear potential.
    std::shared_ptr< FieldAlgorithmBase<spacedim, Value> > inner_field_;
    /// Potential gradient.
    arma::vec::fixed<spacedim> grad_;
    /// Potential constant term.
    double zero_level_;
};


#endif /* FIELD_ADD_POTENTIAL_HH_ */
