/*
 * field_constant.hh
 *
 *  Created on: Dec 15, 2012
 *      Author: jb
 */


#ifndef FIELD_CONSTANT_HH_
#define FIELD_CONSTANT_HH_

#include "system/system.hh"
#include "fields/field_base.hh"
#include "mesh/point.hh"


/**
 * Class representing spatially constant fields.
 *
 */
template <int spacedim, class Value>
class FieldConstant : public FieldBase<spacedim, Value>
{
public:
    typedef typename FieldBase<spacedim, Value>::Point Point;

    /**
     * Default constructor, optionally we need number of components @p n_comp in the case of Vector valued fields.
     */
    FieldConstant(unsigned int n_comp=0);


    static Input::Type::Record input_type;

    /**
     * Return Record for initialization of FieldConstant that is derived from AbstractRecord given by @p a_type
     * and the individual elements of the possible Value (vector, tensor) have Input::Type @p eit.
     */
    static Input::Type::Record get_input_type(Input::Type::AbstractRecord &a_type, typename Value::ElementInputType *eit);

    /**
     * Smart setter from the given value to return.
     */
    FieldConstant<spacedim, Value> &set_value(const typename Value::return_type &val);

    /**
     * This method initialize actual value of the field given from the given Input::Record @p rec.
     */
    virtual void init_from_input(const Input::Record &rec);



    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldConstant();

private:

};


#endif /* FIELD_CONSTANT_HH_ */
