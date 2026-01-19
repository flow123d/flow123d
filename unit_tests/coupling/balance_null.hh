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
 * @file    balance_null.hh
 * @brief   
 */

#ifndef BALANCE_NULL_HH_
#define BALANCE_NULL_HH_


#include "coupling/balance.hh" // for LongIdx






class BalanceNull : public Balance {
public:

    /**
     * Constructor.
     * @param file_prefix  Prefix of output file name.
     * @param mesh         Mesh.
     */
    BalanceNull(const std::string &file_prefix, const Mesh *mesh)
    : Balance(file_prefix, mesh) {}

    /**
     * Destructor.
     */
    ~BalanceNull()
    {}


    /**
     * This method must be called before assembling the matrix for computing mass.
     * It actually erases the matrix.
     */
    void start_mass_assembly(FMT_UNUSED unsigned int quantity_idx)
    {}

    /// Variant of the start_mass_assembly() method for a set of quantities.
    inline void start_mass_assembly(FMT_UNUSED std::vector<unsigned int> q_idx_vec)
    {}

	/**
	 * This method must be called before assembling the matrix and vector for fluxes.
	 * It actually erases the matrix and vector.
	 */
	void start_flux_assembly(FMT_UNUSED unsigned int quantity_idx)
    {}

	/// Variant of the start_flux_assembly() method for a set of quantities.
	inline void start_flux_assembly(FMT_UNUSED std::vector<unsigned int> q_idx_vec)
    {}

    /**
     * This method must be called before assembling the matrix and vectors for sources.
     * It actually erases the matrix and vectors.
     */
    void start_source_assembly(FMT_UNUSED unsigned int quantity_idx)
    {}

    /// Variant of the start_source_assembly() method for a set of quantities.
    inline void start_source_assembly(FMT_UNUSED std::vector<unsigned int> q_idx_vec)
    {}

	
    /// Adds elements into matrix for computing mass.
    void add_mass_values(FMT_UNUSED unsigned int quantity_idx,
        FMT_UNUSED const DHCellAccessor &dh_cell,
	    FMT_UNUSED const LocDofVec &loc_dof_indices,
	    FMT_UNUSED const std::vector<double> &mat_values,
	    FMT_UNUSED double vec_value)
    {}

    /// Adds elements into matrix for computing (outgoing) flux
    void add_flux_values(FMT_UNUSED unsigned int quantity_idx,
        FMT_UNUSED const DHCellSide &side,
        FMT_UNUSED const LocDofVec &loc_dof_indices,
        FMT_UNUSED const std::vector<double> &mat_values,
        FMT_UNUSED double vec_value)
    {}

    /// Adds elements into matrix and vector for computing source.
    void add_source_values(FMT_UNUSED unsigned int quantity_idx,
        FMT_UNUSED unsigned int region_idx,
	    FMT_UNUSED const LocDofVec &loc_dof_indices,
	    FMT_UNUSED const std::vector<double> &mult_mat_values,
	    FMT_UNUSED const std::vector<double> &add_mat_values)
    {}


    /// This method must be called after assembling the matrix for computing mass.
    void finish_mass_assembly(FMT_UNUSED unsigned int quantity_idx)
    {}

    /// Variant of the finish_mass_assembly() method for a set of quantities.
    inline void finish_mass_assembly(FMT_UNUSED std::vector<unsigned int> q_idx_vec)
    {}

    /// This method must be called after assembling the matrix for computing flux.
    void finish_flux_assembly(FMT_UNUSED unsigned int quantity_idx)
    {}

    /// Variant of the finish_flux_assembly() method for a set of quantities.
    inline void finish_flux_assembly(FMT_UNUSED std::vector<unsigned int> q_idx_vec)
    {}

    /// This method must be called after assembling the matrix and vectors for computing source.
    void finish_source_assembly(FMT_UNUSED unsigned int quantity_idx)
    {}

    /// Variant of the finish_source_assembly() method for a set of quantities.
    inline void finish_source_assembly(FMT_UNUSED std::vector<unsigned int> q_idx_vec)
    {}

};





#endif // BALANCE_NULL_HH_
