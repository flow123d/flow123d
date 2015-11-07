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
 * @file    field.hh
 * @brief   
 */

#ifndef FIELD_HH_
#define FIELD_HH_

#include <memory>
using namespace std;

#include <boost/circular_buffer.hpp>

#include "system/exceptions.hh"
#include "input/accessors_forward.hh"
#include "tools/time_marks.hh"
#include "tools/time_governor.hh"

#include "fields/field_common.hh"
#include "fields/field_algo_base.hh"
#include "fields/field_flag.hh"
#include "io/output_time.hh"


namespace IT=Input::Type;

/**
 * @brief Class template representing a field with values dependent on: point, element, and region.
 *
 * By "field" we mean a mapping of a a pair (Point, Time) to a @p Value, where
 * Point is from @p spacedim  dimensional ambient space, Time is real number (set by @p set_time method),
 * and @p Value type representing range of the field, which can be: real scalar, integer scalar (a discrete value),
 * real vector of fixed (compile time) size, real vector of runtime size, or a matrix of fixed dimensions.
 * Extensions to vectors or matrices of integers, or to variable tensors are possible. For vector and matrix values
 * we use classes provided by Armadillo library for linear algebra.
 * The @p Value template parameter should FieldValue<> template, usual choices are:
 * FieldValue<spacedim>::Scalar, FieldValue<spacedim>::Integer, FieldValue<spacedim>::Enum,
 * FieldValue<spacedim>::VectorFixed, FieldValue<spacedim>::TensorFixed
 * deprecated choices: FieldValue<spacedim>::Vector, FieldValue<spacedim>::VectorEnum.
 *
 * This class assign particular fields (instances of descendants of FiledBase) to the regions. It keeps a table of pointers to fields for every possible bulk
 * region index (very same functionality, but for boundary regions is provided by @p BCField class). This class has interface very similar to  FiledBase, however
 * key methods @p value, and @p value_list are not virtual in this class by contrast these methods are inlined to minimize overhead for
 * simplest fields like FieldConstant.
 *
 * TODO: currently it works only for spacedim==3 since we have only mesh in 3D ambient space.
 *
 */
template<int spacedim, class Value>
class Field : public FieldCommon {
public:

    typedef FieldAlgorithmBase<spacedim, Value> FieldBaseType;
    typedef std::shared_ptr< FieldBaseType > FieldBasePtr;
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef Value ValueType;

    static const unsigned int space_dim = spacedim;


    /**
     * Factory class that creates an instance of FieldBase for field
     * with name @p field_name based on data in field descriptor @p rec.
     *
     * Default implementation in method @p create_field just reads key given by
     * @p field_name and creates instance using @p FieldBase<...>::function_factory.
     * Function should return empty SharedField (that is shared_ptr to FieldBase).
     *
     * Implementation of these descendants is necessary:
     * 1) for backward compatibility with old BCD input files
     * 2) for setting pressure values are piezometric head values
     */
    /**
     * Note for future:
     * We pass through parameter @p field information about field that holds the factory which are necessary
     * for interpreting user input and create particular field instance. It would be clearer to pass these information
     * when the factory is assigned to a field. Moreover some information may not be set to field at all but directly passed
     * to the factory.
     */
    class FactoryBase {
    public:
    	/**
    	 * Default method that creates an instance of FieldBase for field.
    	 *
    	 * Reads key given by @p field_name and creates the field instance using
    	 * @p FieldBase<...>::function_factory.
    	 */
    	virtual FieldBasePtr create_field(Input::Record rec, const FieldCommon &field);
    };

    /**
     * Default constructor.
     *
     */
    Field();

    Field(const string &name, bool bc = false);

    /**
     * Copy constructor. Keeps shared history, declaration data, mesh.
     */
    Field(const Field &other);

    /**
     * Assignment operator. Same properties as copy constructor.
     *
     * Question: do we really need this, isn't copy constructor enough?
     * Answer: It is necessary in (usual) case when Field instance is created as the class member
     * but is filled later by assignment possibly from other class.
     */
    Field &operator=(const Field &other);


    /**
     * Returns reference to input type of particular field instance, this is static member @p input_type of the corresponding FieldBase class
     * (i.e. with same template parameters). However, for fields returning "Enum" we have to create whole unique Input::Type hierarchy using following method
     * @p meka_input_tree.
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    const IT::Instance &get_input_type() override;

    IT::Record &get_multifield_input_type() override;


    /**
     * By this method you can allow that the field need not to be set on regions (and times) where the given @p control_field is
     * FieldConstant and has value in given @p value_list. We check this in the set_time method. Through this mechanism we
     * can switch of e.g. boundary data fields according to the type of the boundary condition.
     */
    auto disable_where(
    		const Field<spacedim, typename FieldValue<spacedim>::Enum > &control_field,
    		const vector<FieldEnum> &value_list) -> Field &;



    /**
     * Set mesh pointer and resize region arrays.
     *
     * Implements abstract method.
     */
    void set_mesh(const Mesh &mesh) override;


    /**
     * Direct read access to the table of Field pointers on regions.
     */
    //boost::shared_ptr< FieldBaseType > operator[] (Region reg);

    /**
     * Implementation of @p FieldCommonBase::is_constant().
     */
    bool is_constant(Region reg) override;


    /**
     * Assigns given @p field to all regions in given region set @p domain.
     * Field is added to the history with given time and possibly used in the next call of the set_time method.
     * Caller is responsible for correct construction of given field.
     *
     * Use this method only if necessary.
     *
     * Default time simplify setting steady fields.
     */
    void set_field(const RegionSet &domain, FieldBasePtr field, double time=0.0);

    /**
     * Same as before but the field is first created using FieldBase::function_factory(), from
     * given abstract record accessor @p a_rec.
     */
    void set_field(const RegionSet &domain, const Input::AbstractRecord &a_rec, double time=0.0);

    /**
     * Check that whole field list is set, possibly use default values for unset regions
     * and call set_time for every field in the field list.
     *
     * Returns true if the field has been changed.
     */
    bool set_time(const TimeStep &time) override;

    /**
     * Check that other has same type and assign from it.
     */
    void copy_from(const FieldCommon & other) override;

    /**
     * Implementation of FieldCommonBase::output().
     */
    void output(std::shared_ptr<OutputTime> stream) override;


    /**
     * Returns true, if field is currently set to a time in which it is discontinuous.
     */
    //bool is_jump_time();


    /**
     * Special field values spatially constant. Could allow optimization of tensor multiplication and
     * tensor or vector addition. field_result_ should be set in constructor and in set_time method of particular Field implementation.
     * We return value @p result_none, if the field is not initialized on the region of the given element accessor @p elm.
     */
    inline FieldResult field_result( ElementAccessor<spacedim> &elm) const;

    /**
     * Returns one value in one given point @p on an element given by ElementAccessor @p elm.
     * It returns reference to he actual value in order to avoid temporaries for vector and tensor values.
     */
    virtual typename Value::return_type const &value(const Point &p, const ElementAccessor<spacedim> &elm) const;

    /**
     * Returns std::vector of scalar values in several points at once. The base class implements
     * trivial implementation using the @p value(,,) method. This is not optimal as it involves lot of virtual calls,
     * but this overhead can be negligible for more complex fields as Python of Formula.
     */
    virtual void value_list(const std::vector< Point >  &point_list, const  ElementAccessor<spacedim> &elm,
                       std::vector<typename Value::return_type>  &value_list) const;

    /**
     * Add a new factory for creating Field algorithms on individual regions.
     * The last factory is tried first, the last one is always the default implementation
     * Field<...>::FactoryBase.
     *
     * The Field<...> object keeps a list of such factories. When the instance of a new field algorithm
     * has to be created from the input field descriptor, we pass through the list of factories backward
     * and let factories to create the field algorithm instance from the actual input field descriptor.
     * The first instance (non-null pointer) is used.
     */
    void add_factory(std::shared_ptr<FactoryBase> factory);

protected:

    /**
     * Read input into @p regions_history_ possibly pop some old values from the
     * history queue to keep its size less then @p history_length_limit_.
     */
    void update_history(const TimeStep &time);



    /**
     *  Check that whole field list (@p region_fields_) is set, possibly use default values for unset regions.
     */
    void check_initialized_region_fields_();

    /**************** Shared data **************/

    /// Pair: time, pointer to FieldBase instance
    typedef pair<double, FieldBasePtr> HistoryPoint;
    /// Nearest history of one region.
    typedef boost::circular_buffer<HistoryPoint> RegionHistory;

    struct SharedData {

        /**
         *  History for every region. Shared among copies.
         */
         std::vector< RegionHistory >  region_history_;
    };

    /**************** Data per copy **************/

    std::shared_ptr<SharedData> data_;

	/**
	 * If this pointer is set, turn off check of initialization in the
	 * @p set_time method on the regions where the method @p get_constant_enum_value
	 * of the control field returns value from @p no_check_values_. This
	 * field is private copy, its set_time method is called from the
	 * set_Time method of actual object.
	 */
    typedef Field<spacedim, typename FieldValue<spacedim>::Enum > ControlField;
	std::shared_ptr<ControlField>  no_check_control_field_;

    /**
     * Table with pointers to fields on individual regions.
     */
    std::vector< FieldBasePtr > region_fields_;

    std::vector<std::shared_ptr<FactoryBase> >  factories_;



    template<int dim, class Val>
    friend class MultiField;

};







/****************************************************************************************
 * Inlined methods of Field< ... >
 */

template<int spacedim, class Value>
inline typename Value::return_type const & Field<spacedim,Value>::value(const Point &p, const ElementAccessor<spacedim> &elm) const
{

    ASSERT(this->set_time_result_ != TimeStatus::unknown, "Unknown time status.\n");
    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, idx: %d, field: %s\n", elm.region().id(), elm.region_idx().idx(), name().c_str());
    return region_fields_[elm.region_idx().idx()]->value(p,elm);
}



template<int spacedim, class Value>
inline void Field<spacedim,Value>::value_list(const std::vector< Point >  &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list) const
{
    ASSERT(this->set_time_result_ != TimeStatus::unknown, "Unknown time status.\n");
    ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
    ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, field: %s\n", elm.region().id(), name().c_str());

    region_fields_[elm.region_idx().idx()]->value_list(point_list,elm, value_list);
}







#endif /* FIELD_HH_ */
