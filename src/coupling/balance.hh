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
 * @file    balance.hh
 * @brief   
 */

#ifndef BALANCE_HH_
#define BALANCE_HH_


#include <fstream>              // for ofstream
#include <string>               // for string
#include <vector>               // for vector
#include <unordered_map>        // for unordered_map
#include "mesh/accessors.hh"    // for SideIter
#include "tools/unit_si.hh"    // for UnitSI
#include "input/accessors.hh"   // for Record
#include "petscmat.h"           // for Mat, _p_Mat
#include "petscvec.h"           // for Vec, _p_Vec
#include "system/file_path.hh"  // for FilePath
#include "tools/time_marks.hh"  // for TimeMark, TimeMark::Type
#include "system/index_types.hh" // for LongIdx

class Mesh;
class TimeGovernor;
class DOFHandlerMultiDim;
class DHCellSide;
class DHCellAccessor;
namespace Input {
	namespace Type {
		class Record;
		class Selection;
	}
}





/**
 * Design of balance class - serves as storage and writer.
 * Equations themselves call methods of Balance that add/modify mass, source and flux
 * of various quantities and generate output.
 *
 * One instance of Balance can handle several conservative quantities of the same type
 * (e.g. mass of several substances or their phases).
 *
 * The mass, flux and source are calculated as follows:
 *
 * 	m(q,r) =  ( M'(q) * solution + mv(q) )[r]
 * 	f(q,r) = -( R' * ( F(q) * solution + fv(q) ) )[r]
 * 	s(q,r) =  ( S'(q) * solution + sv(q) )[r]
 *
 * where M' stands for matrix transpose,
 *
 * 	m(q,r)...mass of q-th substance in region r
 * 	f(q,r)...incoming flux of q-th substance in region r
 * 	s(q,r)...source of q-th substance in region r
 *
 * and
 *
 * 	M(q)...region_mass_matrix_		n_dofs x n_bulk_regions
 * 	F(q)...be_flux_matrix_			n_boundary_edges x n_dofs
 * 	S(q)...region_source_matrix_	n_dofs x n_bulk_regions
 * 	SV(q)..region_source_rhs_		n_dofs x n_bulk_regions
 *  mv(q)..region_mass_vec_         n_bulk_regions
 * 	fv(q)..be_flux_vec_				n_boundary_edges
 * 	sv(q)..region_source_vec_    	n_bulk_regions
 * 	R......region_be_matrix_		n_boundary_edges x n_boundary_regions
 * 
 * Remark: Matrix F and the vector fv are such that F*solution+fv produces _outcoming_ fluxes per boundary edge.
 * However we write to output _incoming_ flux due to users' convention and consistently with input interface.
 *
 * Note that it holds:
 *
 * 	sv(q) = column sum of SV(q)
 *
 * Except for that, we also provide information on positive/negative flux and source:
 *
 * 	fp(q,r) = ( R' * EFP(q) )[r],	EFP(q)[e] = max{ 0, ( F(q) * solution + fv(q) )[e] }
 * 	fn(q,r) = ( R' * EFN(q) )[r],	EFN(q)[e] = min{ 0, ( F(q) * solution + fv(q) )[e] }
 * 	sp(q,r) = sum_{i in DOFS } max{ 0, ( S(q)[i,r] * solution[i] + SV(q)[i,r] ) }
 * 	sn(q,r) = sum_{i in DOFS } min{ 0, ( S(q)[i,r] * solution[i] + SV(q)[i,r] ) }
 *
 * where
 *
 * 	fp(q,r)...positive (inward) flux of q-th quantity in region r
 * 	fn(q,r)...negative (outward) flux of q-th quantity in region r
 * 	sp(q,r)...positive source (spring) of q-th quantity in region r
 * 	sn(q,r)...negative source (sink) of q-th quantity in region r
 *
 * Remark: The matrix R is needed only for calculation of signed fluxes.
 * The reason is that to determine sign, we decompose flux to sum of local contributions
 * per each boundary element and check its sign. It is not possible to decompose flux
 * using shape functions, since their normal derivatives may have any sign.
 *
 *
 *
 * Output values (if not relevant, zero is supplied):
 *
 * #time 	bulk_region     quantity 0    0       0        mass source source_in source_out 0               0                 0
 * #time 	boundary_region quantity flux flux_in flux_out 0    0      0         0          0               0                 0
 * #time 	ALL             quantity flux flux_in flux_out mass source source_in source_out integrated_flux integrated_source error
 *
 * error = current_mass - (initial_mass + integrated_source - integrated_flux)
 *
 */
class Balance {
public:

	/**
	 * Class for storing internal data about conservative quantities handled by the Balance object.
	 * In the future we may store additional data or support definition of derived quantities
	 * (e.g. linear combinations of quantities).
	 */
	class Quantity {
	public:

		Quantity(const unsigned int index, const string &name)
			: name_(name),
			  index_(index)
		{}

		/// Name of quantity (for output).
		string name_;

		/// Internal index within list of quantities.
		const unsigned int index_;

	};

	/**
	 * Possible formats of output file.
	 */
	enum OutputFormat
	{
		legacy,//!< legacy
		txt,   //!< csv
		gnuplot//!< gnuplot
	};

	/// Main balance input record type.
	static const Input::Type::Record & get_input_type();

	/// Input selection for file format.
	static const Input::Type::Selection & get_format_selection_input_type();

	/// Set global variable to output balance files into YAML format (in addition to the table format).
	static void set_yaml_output();

	/**
	 * Constructor.
     * @param file_prefix  Prefix of output file name.
     * @param mesh         Mesh.
	 */
	Balance(const std::string &file_prefix, const Mesh *mesh);

	/**
	 * Destructor.
	 */
	~Balance();


    /**
     * Initialize the balance object according to the input.
     * The balance output time marks are set according to the already existing output time marks of the same equation.
     * So, this method must be called after Output::
     *
     * @param in_rec       Input record of balance.
     * @param tg           TimeGovernor of the equation. We need just equation mark type.
     *
     */
    void init_from_input(const Input::Record &in_rec, TimeGovernor &tg);

	/// Setter for units of conserved quantities.
	void units(const UnitSI &unit);

	/// Getter for cumulative_.
	inline bool cumulative() const { return cumulative_; }


	/**
	 * Define a single conservative quantity.
	 * @param name Name of the quantity.
	 * @return     Quantity's internal index.
	 */
	unsigned int add_quantity(const string &name);

	/**
	 * Define a set of conservative quantities.
	 * @param names List of quantities' names.
	 * @return      List of quantities' indices.
	 */
	std::vector<unsigned int> add_quantities(const std::vector<string> &names);

	/**
	 * Allocates matrices and vectors for balance based on DofHandler.
	 * @param dh DofHandler of the corresponding linear system (and solution).
	 * @param max_dofs_per_boundary Number of dofs contributing to one boundary edge.
	 */
	void allocate(const std::shared_ptr<DOFHandlerMultiDim>& dh,
			unsigned int max_dofs_per_boundary);

	/**
	 * Allocates matrices and vectors for balance.
	 * @param n_loc_dofs            Number of solution dofs on the local process.
	 * @param max_dofs_per_boundary Number of dofs contributing to one boundary edge.
	 */
	void allocate(unsigned int n_loc_dofs,
			unsigned int max_dofs_per_boundary);

    /// Returns true if the current time step is marked for the balance output.
    bool is_current();

	/**
	 * This method must be called before assembling the matrix for computing mass.
	 * It actually erases the matrix.
	 */
	void start_mass_assembly(unsigned int quantity_idx);

	/// Variant of the start_mass_assembly() method for a set of quantities.
	inline void start_mass_assembly(std::vector<unsigned int> q_idx_vec)
	{
		for (auto idx : q_idx_vec)
			start_mass_assembly(idx);
	}

	/**
	 * This method must be called before assembling the matrix and vector for fluxes.
	 * It actually erases the matrix and vector.
	 */
	void start_flux_assembly(unsigned int quantity_idx);

	/// Variant of the start_flux_assembly() method for a set of quantities.
	inline void start_flux_assembly(std::vector<unsigned int> q_idx_vec)
	{
		for (auto idx : q_idx_vec)
			start_flux_assembly(idx);
	}

	/**
	 * This method must be called before assembling the matrix and vectors for sources.
	 * It actually erases the matrix and vectors.
	 */
	void start_source_assembly(unsigned int quantity_idx);

	/// Variant of the start_source_assembly() method for a set of quantities.
	inline void start_source_assembly(std::vector<unsigned int> q_idx_vec)
	{
		for (auto idx : q_idx_vec)
			start_source_assembly(idx);
	}

	
	/**
	 * Adds elements into matrix for computing mass.
	 * The mass matrix M is in format [n_dofs x n_bulk_regions] for each quantity.
	 * The mass vector mv is in format [n_bulk_regions] for each quantity.
	 * The region mass is later computed using transpose multiplication (M'*u + mv)[r].
	 * See class @p Balance description above for details.
	 * @param quantity_idx  Index of quantity.
	 * @param dh_cell       Dofhandler cell accessor.
	 * @param loc_dof_indices   Local dof indices (to the solution vector) to be added.
	 * @param mat_values    Values to be added into matrix M.
     * @param vec_value     Value to be added into vector mv.
	 */
	void add_mass_values(unsigned int quantity_idx,
			const DHCellAccessor &dh_cell,
			const LocDofVec &loc_dof_indices,
			const std::vector<double> &mat_values,
			double vec_value);

	/**
	 * Adds elements into matrix for computing (outgoing) flux.
	 * The flux matrix F is in format [n_boundary_edges x n_dofs] for each quantity.
	 * The flux vector fv is in format [n_boundary_edges] for each quantity.
	 * See class @p Balance description above for details.
	 * @param quantity_idx  Index of quantity.
	 * @param side          DHCellSide iterator.
	 * @param loc_dof_indices   Local dof indices (to the solution vector) to be added.
	 * @param mat_values    Values to be added into matrix F.
     * @param vec_value     Value to be added into vector fv.
	 */
	void add_flux_values(unsigned int quantity_idx,
			const DHCellSide &side,
			const LocDofVec &loc_dof_indices,
			const std::vector<double> &mat_values,
			double vec_value);

    /**
	 * Adds elements into matrix and vector for computing source.
	 * The source region matrix S is in format [n_dofs x n_bulk_regions] for each quantity (multiplies solution).
	 * The source region rhs matrix SV is in format [n_dofs x n_bulk_regions] for each quantity (addition).
	 * The region source is then computed as (S'(q) * solution + sv(q))[r], while sv(q) being column sum of SV(q).
	 * The source term for a specific region and quantity is then
	 * @param quantity_idx  Index of quantity.
	 * @param region_idx    Index of bulk region.
	 * @param loc_dof_indices   Local dof indices (to the solution vector) to be added.
	 * @param mult_mat_values   Values to be added into matrix S.
     * @param add_mat_values    Values to be added into matrix SV.
	 */
	void add_source_values(unsigned int quantity_idx,
			unsigned int region_idx,
			const LocDofVec &loc_dof_indices,
			const std::vector<double> &mult_mat_values,
            const std::vector<double> &add_mat_values);


	/// This method must be called after assembling the matrix for computing mass.
	void finish_mass_assembly(unsigned int quantity_idx);

	/// Variant of the finish_mass_assembly() method for a set of quantities.
	inline void finish_mass_assembly(std::vector<unsigned int> q_idx_vec)
	{
		for (auto idx : q_idx_vec)
			finish_mass_assembly(idx);
	}

	/// This method must be called after assembling the matrix and vector for computing flux.
	void finish_flux_assembly(unsigned int quantity_idx);

	/// Variant of the finish_flux_assembly() method for a set of quantities.
	inline void finish_flux_assembly(std::vector<unsigned int> q_idx_vec)
	{
		for (auto idx : q_idx_vec)
			finish_flux_assembly(idx);
	}

	/// This method must be called after assembling the matrix and vectors for computing source.
	void finish_source_assembly(unsigned int quantity_idx);

	/// Variant of the finish_source_assembly() method for a set of quantities.
	inline void finish_source_assembly(std::vector<unsigned int> q_idx_vec)
	{
		for (auto idx : q_idx_vec)
			finish_source_assembly(idx);
	}

	/**
	 * Updates cumulative quantities for balance.
	 * This method can be called in substeps even if no output is generated.
	 * It calculates the sum of source and sum of (incoming) flux over time interval.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 */
	void calculate_cumulative(unsigned int quantity_idx,
			const Vec &solution);

	/**
	 * Calculates actual mass and save it to given vector.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 * @param output_array	Vector of output masses per region.
	 */
	void calculate_mass(unsigned int quantity_idx,
			const Vec &solution,
			vector<double> &output_array);

	/**
	 * Calculates actual mass, incoming flux and source.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 */
	void calculate_instant(unsigned int quantity_idx,
			const Vec &solution);

    /**
	* Adds provided values to the cumulative sources.
	* @param quantity_idx Index of quantity.
	* @param sources Sources per region.
	* @param dt Actual time step.
	*/
	void add_cumulative_source(unsigned int quantity_idx, double source);

	/// Perform output to file for given time instant.
	void output();

private:
	/// Size of column in output (used if delimiter is space)
	static const unsigned int output_column_width = 20;
	/**
	 * Postponed allocation and initialization to allow calling setters in arbitrary order.
	 * In particular we need to perform adding of output times after the output time marks are set.
	 * On the other hand we need to read the input before we make the allocation.
	 *
	 * The initialization is done during the first call of any start_*_assembly method.
	 */
	void lazy_initialize();

	/// Perform output in old format (for compatibility)
	void output_legacy(double time);

	/// Perform output in csv format
	void output_csv(double time, char delimiter, const std::string& comment_string, unsigned int repeat = 0);

	/// Perform output in yaml format
	void output_yaml(double time);

	/// Return part of output represented by zero values. Count of zero values is given by cnt parameter.
	std::string csv_zero_vals(unsigned int cnt, char delimiter);

	/// Print output header
	void format_csv_output_header(char delimiter, const std::string& comment_string);

	/// Format string value of csv output. Wrap string into quotes and if delimiter is space, align text to column.
	std::string format_csv_val(std::string val, char delimiter, bool initial = false);

	/// Format double value of csv output. If delimiter is space, align text to column.
	std::string format_csv_val(double val, char delimiter, bool initial = false);

    /** Computes unique id of local boundary edge from local element index and element side index
     * 
     * @param side is a side of locally owned element
    */
    inline LongIdx get_boundary_edge_uid(SideIter side)
    { return 4*side->elem_idx() + side->side_idx();}    // 4 is maximum of sides per element

	//**********************************************

	static bool do_yaml_output_;

	/// Allocation parameters. Set by the allocate method used in the lazy_initialize.
	unsigned int n_loc_dofs_par_;
	unsigned int n_loc_dofs_seq_;
    unsigned int max_dofs_per_boundary_;


    /// Save prefix passed in in constructor.
    std::string file_prefix_;

    /// File path for output_ stream.
    FilePath balance_output_file_;

    /// Handle for file for output in given OutputFormat of balance and total fluxes over individual regions and region sets.
    std::ofstream output_;

    // The same as the previous case, but for output in YAML format.
    std::ofstream output_yaml_;

    /// Format of output file.
    OutputFormat output_format_;

    /// Names of conserved quantities.
    std::vector<Quantity> quantities_;

    const Mesh *mesh_;

    /// Units of conserved quantities.
    UnitSI units_;


    /// Matrices for calculation of mass (n_dofs x n_bulk_regions).
    Mat *region_mass_matrix_;

    /// Matrices for calculation of flux (n_boundary_edges x n_dofs).
    Mat *be_flux_matrix_;

    /// Matrices for calculation of source (n_dofs x n_bulk_regions).
    Mat *region_source_matrix_;

    /// Matrices for calculation of signed source (n_dofs x n_bulk_regions).
    Mat *region_source_rhs_;

    /// Vectors for calculation of flux (n_boundary_edges).
    Vec *be_flux_vec_;
    
    /// Vectors for calculation of mass (n_bulk_regions).
    Vec *region_mass_vec_;

    /** Maps unique identifier of (local bulk element idx, side idx) returned by @p get_boundary_edge_uid(side)
     * to local boundary edge.
     * Example usage:
     *     be_id = be_id_map_[get_boundary_edge_uid(side)]
     */
    std::unordered_map<LongIdx, unsigned int> be_id_map_;
    
    /// Maps local boundary edge to its region boundary index.
    std::vector<unsigned int> be_regions_;

    /// Offset for local part of vector of boundary edges.
    int be_offset_;


    // Vectors storing mass and balances of fluxes and volumes.
    // substance, phase, region
    std::vector<std::vector<double> > fluxes_;
    std::vector<std::vector<double> > fluxes_in_;
    std::vector<std::vector<double> > fluxes_out_;
    std::vector<std::vector<double> > masses_;
    std::vector<std::vector<double> > sources_in_;
    std::vector<std::vector<double> > sources_out_;

    // Sums of the above vectors over phases and regions
    std::vector<double> sum_fluxes_;
    std::vector<double> sum_fluxes_in_;
    std::vector<double> sum_fluxes_out_;
    std::vector<double> sum_masses_;
    std::vector<double> sum_sources_;
    std::vector<double> sum_sources_in_;
    std::vector<double> sum_sources_out_;
    std::vector<double> initial_mass_;

	// time integrated quantities
    std::vector<double> integrated_sources_;
    std::vector<double> integrated_fluxes_;
    std::vector<double> increment_fluxes_;
    std::vector<double> increment_sources_;

	/// time of last calculated balance
	double last_time_;

	/// TimeMark type for balance output of particular equation.
	TimeMark::Type balance_output_type_;

	/// TimeMark type for output of particular equation.
	TimeMark::Type output_mark_type_;

	/// true before calculating the mass at initial time, otherwise false
	bool initial_;

	/// if true then cumulative balance is computed
	bool cumulative_;

	/// true before allocating necessary internal structures (Petsc matrices etc.)
	bool allocation_done_;

	/// If the balance is on. Balance is off in the case of no balance output time marks.
    bool balance_on_;

    /// Add output time marks to balance output time marks.
    bool add_output_times_;


	/// MPI rank.
	int rank_;

	/// hold count of line printed into output_
	unsigned int output_line_counter_;

	/// marks whether YAML output has printed header
	bool output_yaml_header_;

    /// Record for current balance
    Input::Record input_record_;

    const TimeGovernor *time_;


};





#endif // BALANCE_HH_
