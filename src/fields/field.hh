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

#include <stdio.h>                                     // for sprintf
#include <string.h>                                    // for memcpy
#include <algorithm>                                   // for find, min
#include <boost/circular_buffer.hpp>
#include <memory>                                      // for dynamic_pointe...
#include <new>                                         // for operator new[]
#include <ostream>                                     // for basic_ostream:...
#include <string>                                      // for basic_string
#include <utility>                                     // for pair
#include <vector>                                      // for vector
#include <armadillo>
#include "fields/field_algo_base.hh"                   // for FieldAlgorithm...
#include "fields/field_common.hh"                      // for FieldCommon::T...
#include "fields/field_values.hh"                      // for FieldValue<>::...
#include "fields/field_value_cache.hh"                 // for FieldValueCache
#include "input/accessors.hh"                          // for ExcTypeMismatch
#include "input/accessors_impl.hh"                     // for Record::opt_val
#include "input/factory_impl.hh"                       // for Factory::create
#include "input/input_exception.hh"                    // for FieldCommon::E...
#include "input/storage.hh"                            // for ExcStorageType...
#include "input/type_base.hh"                          // for Array
#include "input/type_generic.hh"                       // for Instance
#include "input/type_record.hh"                        // for Record::ExcRec...
#include "input/input_exception.hh"                    // for Input::Exception
#include "io/output_time.hh"                           // for OutputTime
#include "mesh/elements.h"                             // for Element::dim
#include "mesh/region.hh"                              // for RegionDB::ExcU...
#include "system/asserts.hh"                           // for Assert, ASSERT
#include "system/exc_common.hh"                        // for ExcAssertMsg
#include "system/exceptions.hh"                        // for ExcAssertMsg::...
#include "system/global_defs.h"                        // for OLD_ASSERT, msg
#include "tools/time_governor.hh"                      // for TimeStep

class Mesh;
class Observe;
class EvalPoints;
class BulkPoint;
class EdgePoint;
class CouplingPoint;
class BoundaryPoint;
class FieldSet;
template <int spacedim> class ElementAccessor;
template <int spacedim, class Value> class FieldFE;
namespace detail
{
    template< typename CALLABLE, typename TUPLE, int INDEX >
    struct model_cache_item;
}

using namespace std;
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
 * FieldValue<spacedim>::VectorFixed, FieldValue<spacedim>::TensorFixed.
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

    	/**
    	 * Check if Input::Record accessor contains data of field given by input_name.
    	 *
    	 * Returns true when ever the method create_field returns non-null pointer, otherwise returns false.
    	 */
    	virtual bool is_active_field_descriptor(const Input::Record &in_rec, const std::string &input_name);
    };

    /**
     * Default constructor.
     *
     */
    Field();

    Field(const string &name, bool bc = false);

    /**
     * Constructor that must be used for create of MultiField components.
     *
     * Set parameters @p component_index_, @p shared_->input_name_ and @p name_.
     * Parameter name_ of Field is consisted of component name and MultiField name.
     */
    Field(unsigned int component_index, string input_name, string name = "", bool bc = false);

    /**
     * Copy constructor. Keeps shared history, declaration data, mesh.
     */
    Field(const Field &other);

    /**
     * Assignment operator. Same properties as copy constructor, but class member name_ is not copied.
     *
     * Question: do we really need this, isn't copy constructor enough?
     * Answer: It is necessary in (usual) case when Field instance is created as the class member
     * but is filled later by assignment possibly from other class.
     * TODO: operator can be merged with copy constructor, but we must provide to set correct value
     * of name in method copy_from
     */
    Field &operator=(const Field &other);


    /// Return appropriate value to BulkPoint in FieldValueCache
    typename Value::return_type operator() (BulkPoint &p);


    /// Return appropriate value to EdgePoint in FieldValueCache
    typename Value::return_type operator() (EdgePoint &p);


    /// Return appropriate value to CouplingPoint in FieldValueCache
    typename Value::return_type operator() (CouplingPoint &p);


    /// Return appropriate value to BoundaryPoint in FieldValueCache
    typename Value::return_type operator() (BoundaryPoint &p);


    /**
     * Returns reference to input type of particular field instance, this is static member @p input_type of the corresponding FieldBase class
     * (i.e. with same template parameters). However, for fields returning "Enum" we have to create whole unique Input::Type hierarchy using following method
     * @p make_input_tree.
     * every instance since every such field use different Selection for initialization, even if all returns just unsigned int.
     */
    IT::Instance get_input_type() override;

    IT::Array get_multifield_input_type() override;


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
    //std::shared_ptr< FieldBaseType > operator[] (Region reg);

    /**
     * Implementation of @p FieldCommonBase::is_constant().
     * See also Field<>::field_result which provide better information about special field values.
     */
    bool is_constant(Region reg) override;

    /**
     * Assigns given @p field to all regions in region set given by @p region_set_names.
     * Field is added to the history with given time and possibly used in the next call of the set_time method.
     * Caller is responsible for correct construction of given field.
     *
     * Use this method only if necessary.
     */
    void set(FieldBasePtr field, double time, std::vector<std::string> region_set_names = {"ALL"});

    /**
     * Same as before but the field is first created using FieldBase::function_factory(), from
     * given abstract record accessor @p a_rec.
     */
    void set(const Input::AbstractRecord &a_rec, double time, std::vector<std::string> region_set_names = {"ALL"});

    /**
     * Check that whole field list is set, possibly use default values for unset regions
     * and call set_time for every field in the field list.
     *
     * Returns true if the field has been changed.
     */
    bool set_time(const TimeStep &time, LimitSide limit_side) override;

    /**
     * Check that other has same type and assign from it.
     */
    void copy_from(const FieldCommon & other) override;

    /**
     * Implementation of FieldCommonBase::output().
     */
    void field_output(std::shared_ptr<OutputTime> stream, OutputTime::DiscreteSpaceFlags type) override;

    /**
     * Implementation of FieldCommonBase::observe_output().
     */
    void observe_output(std::shared_ptr<Observe> observe) override;

    /**
     * Returns true, if field is currently set to a time in which it is discontinuous.
     */
    //bool is_jump_time();


    /**
     * @brief Indicates special field states.
     *
     * Return possible values from the enum @p FieldResult, see description there.
     * The initial state is @p field_none, if the field is correctly set on all regions of the @p region_set given as parameter
     * we return state @p field_other
     * - Special field values spatially constant. Could allow optimization of tensor multiplication and
     * tensor or vector addition. field_result_ should be set in constructor and in set_time method of particular Field implementation.
     * We return value @p result_none, if the field is not initialized on the region of the given element accessor @p elm.
     * Other possible results are: result_zeros, result_eye, result_ones, result_constant, result_other
     * see @p FieldResult for explanation.
     */
    FieldResult field_result( RegionSet region_set) const override;

    /**
     * Return value of input type attribute 'field_value_shape' that is appended to the
     * input type of this field in FieldSet::make_field_descriptor_type and also to the output field selection
     * created in EquationOutput::make_output_type.
     * This attribute is used by GeoMop to have semantics of the input and output field data.
     *
     * Attribute value is a valid JSON (and/or flow style YAML) with keys:
     * 'subfields' - value True for multifields, False or not presented for single value fields
     * 'shape' - [ NRows, Ncols] ... given by FieldValue
     * 'type' - <element type> (Double or Integer) ... given by FieldValue
     * 'limit' - bounds of the field values.
     *
     */
    std::string get_value_attribute() const override;

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
    virtual void value_list(const Armor::array &point_list, const  ElementAccessor<spacedim> &elm,
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

    void set_input_list(const Input::Array &list, const TimeGovernor &tg) override;

    /**
     * Interpolate given field into output discrete @p space_type and store the values
     * into storage of output time @p stream for postponed output.
     */
    void compute_field_data(OutputTime::DiscreteSpaceFlags space_type, std::shared_ptr<OutputTime> stream);

    /// Implements FieldCommon::cache_allocate
    void cache_reallocate(const ElementCacheMap &cache_map, unsigned int region_idx) const override;

    /// Implements FieldCommon::cache_update
    void cache_update(ElementCacheMap &cache_map, unsigned int region_patch_idx) const override;

    /// Implements FieldCommon::value_cache
    FieldValueCache<double> * value_cache() override;

    /// Implements FieldCommon::value_cache
    const FieldValueCache<double> * value_cache() const override;

    /**
     * Implementation of FieldCommon::set_dependency().
     */
    std::vector<const FieldCommon *> set_dependency(FieldSet &field_set, unsigned int i_reg) const override;

protected:

    /// Return item of @p value_cache_ given by i_cache_point.
    typename Value::return_type operator[] (unsigned int i_cache_point) const;

    /**
     * Read input into @p regions_history_ possibly pop some old values from the
     * history queue to keep its size less then @p history_length_limit_.
     */
    void update_history(const TimeStep &time);

    /// Fills acutally the data cache with field values, used in @p compute_field_data
    void fill_data_cache(OutputTime::DiscreteSpace space_type,
                         std::shared_ptr<OutputTime> stream,
                         std::shared_ptr<ElementDataCache<typename Value::element_type>> data_cache);

    /**
     *  Check that whole field list (@p region_fields_) is set, possibly use default values for unset regions.
     */
    void check_initialized_region_fields_();

    /**
     * Check that the field is in fact FieldFE set on all bulk regions, return shared pointer to that FieldFE or NULL
     * if the Field is not FieldFE.
     */
    std::shared_ptr< FieldFE<spacedim, Value> > get_field_fe();

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

    /**
     * Field value data cache
     *
     * Data is ordered like three dimensional table. The highest level is determinated by subsets,
     * those data ranges are holds in subset_starts. Data block size of each subset is determined
     * by number of eval_points (of subset) and maximal number of stored elements.
     * The table is allocated to hold all subsets, but only those marked in used_subsets are updated.
     * Order of subsets is same as in eval_points.
     */
    mutable FieldValueCache<typename Value::element_type> value_cache_;



    template<int dim, class Val>
    friend class MultiField;

    template< typename CALLABLE, typename TUPLE, int INDEX >
    friend struct detail::model_cache_item;

};







/****************************************************************************************
 * Inlined methods of Field< ... >
 */

template<int spacedim, class Value>
inline typename Value::return_type const & Field<spacedim,Value>::value(const Point &p, const ElementAccessor<spacedim> &elm) const
{

    ASSERT(this->set_time_result_ != TimeStatus::unknown)(this->name()).error("Unknown time status.\n");
	OLD_ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
	OLD_ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, idx: %d, field: %s\n", elm.region().id(), elm.region_idx().idx(), name().c_str());
    return region_fields_[elm.region_idx().idx()]->value(p,elm);
}



template<int spacedim, class Value>
inline void Field<spacedim,Value>::value_list(const Armor::array &point_list, const ElementAccessor<spacedim> &elm,
                   std::vector<typename Value::return_type>  &value_list) const
{
    ASSERT(this->set_time_result_ != TimeStatus::unknown)(this->name()).error("Unknown time status.\n");
	OLD_ASSERT(elm.region_idx().idx() < region_fields_.size(), "Region idx %u out of range %lu, field: %s\n",
           elm.region_idx().idx(), (unsigned long int) region_fields_.size(), name().c_str());
	OLD_ASSERT( region_fields_[elm.region_idx().idx()] ,
    		"Null field ptr on region id: %d, field: %s\n", elm.region().id(), name().c_str());
    ASSERT_DBG(point_list.n_rows() == spacedim && point_list.n_cols() == 1).error("Invalid point size.\n");

    region_fields_[elm.region_idx().idx()]->value_list(point_list,elm, value_list);
}






#endif /* FIELD_HH_ */
