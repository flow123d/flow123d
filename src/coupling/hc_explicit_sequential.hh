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
 * @brief
 *
 *  @author Jan Brezina
 */

#ifndef HC_EXPLICIT_SEQUENTIAL_HH_
#define HC_EXPLICIT_SEQUENTIAL_HH_

#include "input/input_type.hh"
#include "input/accessors.hh"

#include "coupling/equation.hh"

class DarcyFlowMH;
class DarcyFlowMHOutput;
class Mesh;
class EquationBase;
class AdvectionProcessBase;
class MaterialDatabase;


/**
 * TODO: should be derived from EquationBase in order to chain couplings
 */
class CouplingBase {
public:
    static Input::Type::AbstractRecord & get_input_type();

};


/**
 * @brief Class for solution of steady or unsteady flow with sequentially coupled explicit transport.
 *
 */
class HC_ExplicitSequential : public CouplingBase {
public:
    static const Input::Type::Record & get_input_type();

    HC_ExplicitSequential(Input::Record in_record);
    void run_simulation();
    ~HC_ExplicitSequential();

private:

    static const int registrar;

    /// mesh common to darcy flow and transport
    Mesh *mesh;

    /// steady or unsteady water flow simulator based on MH scheme
    std::shared_ptr<DarcyFlowMH> water;

    /// explicit transport with chemistry through operator splitting
    std::shared_ptr<AdvectionProcessBase> secondary_eq;

};

#endif /* HC_EXPLICIT_SEQUENTIAL_HH_ */
