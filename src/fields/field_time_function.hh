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
 * @file    field_time_function.hh
 * @brief
 */

#ifndef FIELD_TIME_FUNCTION_HH_
#define FIELD_TIME_FUNCTION_HH_

#include <string.h>                        // for memcpy
#include <boost/exception/info.hpp>        // for error_info::error_info<Tag...
#include <new>                             // for operator new[]
#include <ostream>                         // for operator<<
#include <string>                          // for operator<<, string
#include <utility>                         // for pair
#include <armadillo>
#include "fields/field_algo_base.hh"       // for FieldAlgorithmBase
#include "fields/field_constant.hh"        // for FieldConstant
#include "fields/field_values.hh"          // for FieldValue<>::TensorFixed
#include "input/accessors.hh"              // for ExcAccessorForNullStorage
#include "input/accessors_impl.hh"         // for Record::val
#include "input/storage.hh"                // for ExcStorageTypeMismatch
#include "input/type_record.hh"            // for Record::ExcRecordKeyNotFound
#include "tools/time_governor.hh"          // for TimeStep
class UnitSI;

/**
 * Class representing spatially fields defined by time-dependent function.
 *
 */
template <int spacedim, class Value>
class FieldTimeFunction : public FieldConstant<spacedim, Value>
{
public:
    typedef typename FieldAlgorithmBase<spacedim, Value>::Point Point;
    typedef FieldAlgorithmBase<spacedim, Value> FactoryBaseType;

    /**
     * Return Record for initialization of FieldTimeFunction that is derived from Abstract given by @p a_type
     * and the individual elements of the possible Value (vector, tensor) in given times have Input::Type @p eit.
     */
    static const Input::Type::Record & get_input_type();


    /**
     * Default constructor, optionally we need number of components @p n_comp in the case of Vector valued fields.
     */
    FieldTimeFunction(unsigned int n_comp=0);

    /**
     * This method initialize actual value of the field given from the given Input::Record @p rec.
     *
     * TODO: after removing support for vector valued FieldConstant we can merge this method
     * with FieldConstant::init_from_input and move initizaliation of FieldConstant value
     * to set_time method.
     */
    void init_from_input(const Input::Record &rec, const struct FieldAlgoBaseInitData& init_data) override;

    /**
     * Set time and init value_.
     */
    bool set_time(const TimeStep &time) override;

private:
    /// Registrar of class to factory
    static const int registrar;

    /// Accessor to Input::Record
    Input::Record in_rec_;

    /**
     * Initialization data of field. Necessary for check limits.
     *
     * TODO: Temporary solution will be replaced with shared_ptr to field data in FieldAlgoBase
     */
    std::string field_name_;
	UnitSI &unit_si_;
	std::pair<double, double> limits_;


};

#endif /* FIELD_TIME_FUNCTION_HH_ */
