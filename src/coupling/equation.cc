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
: mesh_(NULL), mat_base(NULL), time_(NULL),
  equation_mark_type_(TimeGovernor::marks().new_mark_type()), //creating mark type for new equation
  input_record_()
{}



EquationBase::EquationBase(Mesh &mesh, MaterialDatabase &mat_base,const  Input::Record in_rec)
: mesh_(&mesh),
  mat_base(&mat_base),
  time_(NULL),
  equation_mark_type_(TimeGovernor::marks().new_mark_type()), //creating mark type for new equation
  input_record_(in_rec)
{}


/*****************************************************************************************
 * Implementation of EqDataBase
 */

namespace IT=Input::Type;


EqDataBase::EqDataBase(const std::string& eq_name)
: equation_name_(eq_name)
{}



void EqDataBase::add_field( FieldCommonBase *field, const string &name, const string &desc, Input::Type::Default d_val) {
    field->set_name( name );
    field->set_desc( desc );
    field->set_default( d_val );
    field_list.push_back(field);
}



IT::Record EqDataBase::generic_input_type(const string &eq_class_name, const string &desc, bool bc_regions) {
    string rec_name = eq_class_name + (bc_regions ? "_BoundaryData" : "_BulkData");
    IT::Record rec = IT::Record(rec_name, desc)
                     .declare_key("region", IT::String(), "Region label")
                     .declare_key("rid", IT::Integer(0), "Region ID (alternative to region label)" );

    BOOST_FOREACH(FieldCommonBase * field, field_list)
        if (bc_regions == field->is_bc()) {
            if (field->is_enum_valued())
                rec.declare_key(field->name(), field->make_input_tree(), field->get_default(), field->desc() );
            else
                rec.declare_key(field->name(), field->get_input_type(), field->get_default(), field->desc() );
        }


    // intentionally we do not call finish here in order to allow adding keys after the generic ones
    // finish should be called at global level through Lazyhttp://stackoverflow.com/questions/4786649/are-variadic-macros-nonstandardTypes

    return rec;
}



void EqDataBase::init_from_input(Input::Array bulk_list, Input::Array bc_list) {
    for(Input::Iterator<Input::Record> it=bulk_list.begin<Input::Record>(); it != bulk_list.end(); ++it)
        init_from_input_one_region(*it, false);
    for(Input::Iterator<Input::Record> it=bc_list.begin<Input::Record>(); it != bc_list.end(); ++it)
        init_from_input_one_region(*it, true);
}



Region EqDataBase::init_from_input_one_region(Input::Record rec, bool bc_regions) {
    Input::Iterator<string> it = rec.find<string>("region");
    Region reg;

    // get the region
    if (it) {
        // try find region by label
        reg = Region::db().find_label(*it);
        if (! reg.is_valid() ) xprintf(UsrErr, "Unknown region with label: '%s'\n", (*it).c_str());
    } else {
        // try find region by ID
        Input::Iterator<unsigned int> id_it = rec.find<unsigned int>("rid");
        reg = Region::db().find_id(*id_it);
        if (! reg.is_valid() ) xprintf(UsrErr, "Unknown region with id: '%d'\n", *id_it);
    }

    // init all fields on this region
    BOOST_FOREACH(FieldCommonBase * field, field_list) {
        if (bc_regions == field->is_bc()) {
            Input::Iterator<Input::AbstractRecord> field_it = rec.find<Input::AbstractRecord>(field->name());
            DBGMSG("%s %d\n", field->name().c_str(), int(field_it) );
            if (field_it) {
                field->init_from_input(reg,*field_it);
            }
        }
    }

    return reg;
}


/*****************************************************************************************
 * Implementation of EquationNothing
 */

EquationNothing::EquationNothing(Mesh &mesh, MaterialDatabase &mat_base)
: EquationBase(mesh, mat_base, Input::Record() )
{}
