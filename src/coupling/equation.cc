/*!
 *
 * Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
 *
 * Please make a following refer to Flow123d on your project site if you use the program for any purpose,
 * especially for academic research:
 * Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
 *
 * This program is free software; you can redistribute it and/or modify it under the terms
 * of the GNU General Public License version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
 *
 *
 * $Id$
 * $Revision$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 * @file
 * @brief Abstract base class for equation clasess.
 *
 *  @author Jan Brezina
 */

#include <petscmat.h>
#include "time_governor.hh"

#include "equation.hh"
#include "system/system.hh"
#include "input/accessors.hh"

#include <boost/foreach.hpp>
#include "fields/field_base.hh"

/*****************************************************************************************
 * Implementation of EqBase
 */

EquationBase::EquationBase()
: mesh_(NULL), time_(NULL),
  equation_mark_type_(TimeGovernor::marks().new_mark_type()), //creating mark type for new equation
  input_record_()
{}



EquationBase::EquationBase(Mesh &mesh, const  Input::Record in_rec)
: mesh_(&mesh),
  time_(NULL),
  equation_mark_type_(TimeGovernor::marks().new_mark_type()), //creating mark type for new equation
  input_record_(in_rec)
{}


/*****************************************************************************************
 * Implementation of EqDataBase
 */

namespace IT=Input::Type;


EqDataBase::EqDataBase(const std::string& eq_name)
: mesh_(NULL),
  equation_name_(eq_name)
{
    if (equation_name_ == "") xprintf(PrgErr, "You have to provide non-empty equation name when constructing EqDataBase.\n");
}



void EqDataBase::add_field( FieldCommonBase *field, const string &name, const string &desc, Input::Type::Default d_val) {
    field->set_name( name );
    field->set_desc( desc );
    field->set_default( d_val );
    field_list.push_back(field);
}



IT::Record EqDataBase::generic_input_type(bool bc_regions) {
    string rec_name, description;

    if (bc_regions) {
        rec_name = equation_name_ + "_BoundaryData";
        description = "Record to set BOUNDARY fields of the equation '" + equation_name_ + "'.\n"
                "The fields are set only on the domain specified by one of the keys: 'region', 'rid', 'r_set'\n"
                "and after the time given by the key 'time'. The field setting can be overridden by\n"
                " any " + rec_name + " record that comes later in the boundary data array.";
    } else {
        rec_name = equation_name_ + "_BulkData";
        description = "Record to set BULK fields of the equation '" + equation_name_ + "'.\n"
                "The fields are set only on the domain specified by one of the keys: 'region', 'rid', 'r_set'\n"
                "and after the time given by the key 'time'. The field setting can be overridden by\n"
                " any " + rec_name + " record that comes later in the bulk data array.";
    }
    IT::Record rec = IT::Record(rec_name, description)
                     .declare_key("r_set", IT::String(), "Name of region set where to set fields.")
                     .declare_key("region", IT::String(), "Label of the region where to set fields. ")
                     .declare_key("rid", IT::Integer(0), "ID of the region where to set fields." )
                     .declare_key("time", IT::Double(0.0), IT::Default("0.0"),
                             "Apply field setting in this record after this time.\n"
                             "These times have to form an increasing sequence.");

    BOOST_FOREACH(FieldCommonBase * field, field_list)
        if (bc_regions == field->is_bc()) {
            if (field->is_enum_valued())
                rec.declare_key(field->name(), field->make_input_tree(), field->get_default(), field->desc() );
            else
                rec.declare_key(field->name(), field->get_input_type(), field->get_default(), field->desc() );
        }

    return rec;
}



IT::Record EqDataBase::boundary_input_type() {
    return generic_input_type(true);
}


IT::Record EqDataBase::bulk_input_type() {
    return generic_input_type(false);
}



void EqDataBase::set_time(const TimeGovernor &time) {
    /*
     * - read records from arrays until we reach greater time then actual
     * - update fields (delete the previous, use make factory for the new one.
     */
    set_time(time, boundary_input_array_, boundary_it_, true);
    set_time(time, bulk_input_array_, bulk_it_, false);
}



void EqDataBase::set_time(const TimeGovernor &time, Input::Array &list, Input::Iterator<Input::Record> &it, bool bc_regions) {
    // read input up to given time
    if (list.size() != 0) {
        while( it != list.end() && time.ge( it->val<double>("time") ) ) {
            if (bc_regions)
            {
                read_boundary_list_item(*it);
            }
            else read_bulk_list_item(*it);
            ++it;
        }
    }
    // check validity of fields and set current time
    BOOST_FOREACH(FieldCommonBase * field, field_list) 
    {
      if (bc_regions == field->is_bc())
        field->set_time( time.t() );
    }
}



void EqDataBase::set_mesh(Mesh *mesh) {
    mesh_=mesh;
    BOOST_FOREACH(FieldCommonBase * field, field_list) field->set_mesh( mesh_ );
}


void EqDataBase::check_times(Input::Array &list) {
    double time,last_time=0.0;
    if (list.size() == 0) return;
    for( Input::Iterator<Input::Record> it = list.begin<Input::Record>(); it != list.end(); ++it) {
        time = it->val<double>("time");
        if (time < last_time) xprintf(UsrErr, "Time %f in bulk data of equation '%s' is smaller then the previous time %f.\n",
                time, equation_name_.c_str(), last_time );
        last_time=time;
    }
}

void EqDataBase::init_from_input(Input::Array bulk_list, Input::Array bc_list) {
    if (mesh_ == NULL) xprintf(PrgErr, "The mesh pointer wasn't set in the EqData of equation '%s'.\n", equation_name_.c_str());

    bulk_input_array_ = bulk_list;
    boundary_input_array_ = bc_list;

    if (bulk_input_array_.size() !=0) {
        check_times(bulk_input_array_);
        bulk_it_ = bulk_input_array_.begin<Input::Record>();
    }
    if (boundary_input_array_.size() !=0) {
        check_times(boundary_input_array_);
        boundary_it_ = boundary_input_array_.begin<Input::Record>();
    }
}




RegionSet EqDataBase::read_boundary_list_item(Input::Record rec) {
    return read_list_item(rec, true);
}



RegionSet EqDataBase::read_bulk_list_item(Input::Record rec) {
    return read_list_item(rec, false);
}




RegionSet EqDataBase::read_list_item(Input::Record rec, bool bc_regions) {
    RegionSet domain;
    std::string name;
    unsigned int id;
    if (rec.opt_val("r_set", name)) {     
        domain = mesh_->region_db().get_region_set(name);

    } else if (rec.opt_val("region", name)) {
        domain.push_back( mesh_->region_db().find_label(name) );    // try find region by label
        if (! domain[0].is_valid() ) xprintf(UsrErr, "Unknown region with label: '%s'\n", name.c_str());

    } else if (rec.opt_val("rid", id)) {
        domain.push_back( mesh_->region_db().find_id(id) );         // try find region by ID
        if (! domain[0].is_valid() ) xprintf(UsrErr, "Unknown region with id: '%d'\n", id);

    } else {
        if (bc_regions) {
            xprintf(UsrErr, "Missing domain specification in BOUNDARY record, item %d, equation %s\n", boundary_it_.idx(), equation_name_.c_str());
        } else {
            xprintf(UsrErr, "Missing domain specification in BULK record, item %d, equation %s\n", bulk_it_.idx(), equation_name_.c_str());
        }
    }
    
    // init all fields on this region
    BOOST_FOREACH(FieldCommonBase * field, field_list) {
        if (bc_regions == field->is_bc()) {
            Input::Iterator<Input::AbstractRecord> field_it = rec.find<Input::AbstractRecord>(field->name());
            if (field_it) {
                field->set_from_input(domain, *field_it);
                //DBGMSG("set mesh field: %s\n", field->name().c_str() );
                field->set_mesh(mesh_);
            }
        }
    }

    return domain;
}


EqDataBase::~EqDataBase() {}

/*****************************************************************************************
 * Implementation of EquationNothing
 */

EquationNothing::EquationNothing(Mesh &mesh)
: EquationBase(mesh, Input::Record() )
{}
