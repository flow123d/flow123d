/*
 * field_constant.hh
 *
 *  Created on: Dec 15, 2012
 *      Author: jb
 */


#ifndef FIELD_FE_HH_
#define FIELD_FE_HH_

#include "petscmat.h"
#include "system/system.hh"
#include "fields/field_algo_base.hh"
#include "mesh/point.hh"
#include "fem/dofhandler.hh"
#include "fem/mapping.hh"
#include "input/factory.hh"



/**
 * Class representing fields given by finite element approximation.
 *
 */
template <int spacedim, class Value>
class FieldFE : public FieldAlgorithmBase<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    /**
     * Default constructor, optionally we need number of components @p n_comp in the case of Vector valued fields.
     */
    FieldFE(unsigned int n_comp=0);

    /**
     * Setter for the finite element data. The mappings are required for computation of local coordinates.
     * @param dh   Dof handler.
     * @param map1 1D mapping.
     * @param map2 2D mapping.
     * @param map3 3D mapping.
     * @param data Vector of dof values.
     */
    void set_fe_data(const DOFHandlerMultiDim *dh,
    		Mapping<1,3> *map1,
    		Mapping<2,3> *map2,
    		Mapping<3,3> *map3,
    		const Vec *data);

    /**
     * Returns one value in one given point. ResultType can be used to avoid some costly calculation if the result is trivial.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm);

    /**
     * Returns std::vector of scalar values in several points at once.
     */
    virtual void value_list (const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list);


    virtual ~FieldFE();

private:

    const DOFHandlerMultiDim *dh_;
    double *data_;
    const Vec *data_vec_;
    unsigned int *dof_indices;

    Mapping<1,3> *map1_;
    Mapping<2,3> *map2_;
    Mapping<3,3> *map3_;

    /// Registrar of class to factory
    static const int registrar;
};


#endif /* FIELD_FE_HH_ */
