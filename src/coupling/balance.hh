#ifndef BALANCE_HH_
#define BALANCE_HH_


#include "la/distribution.hh"
#include "transport/substance.hh"
#include "petscmat.h"
#include "fields/unit_si.hh"

class RegionDB;





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
 * 	m(q,c,r) = ( M'(q,c) * solution )[r]
 * 	f(q,c,r) = ( R' * ( F(q,c) * solution + fv(q,c) ) )[r]
 * 	s(q,c,r) = ( S'(q,c) * solution + sv(q,c) )[r]
 *
 * where M' stands for matrix transpose,
 *
 * 	m(q,c,r)...mass of q-th substance's c-th component in region r
 * 	f(q,c,r)...flux of q-th substance's c-th component in region r
 * 	s(q,c,r)...source of q-th substance's c-th component in region r
 *
 * and
 *
 * 	M(q,c)...region_mass_matrix_[q*quantities_.size() + c]		n_dofs x n_bulk_regions
 * 	F(q,c)...be_flux_matrix_[q*quantities_.size() + c]			n_boundary_edges x n_dofs
 * 	S(q,c)...region_source_matrix_[q*quantities_.size() + c]	n_dofs x n_bulk_regions
 * 	SV(q,c)..region_source_rhs_[q*quantities_.size() + c]		n_dofs x n_bulk_regions
 * 	fv(q,c)..be_flux_vec_[q*quantities_.size() + c]				n_boundary_edges
 * 	sv(q,c)..region_source_vec_[q*quantities_.size() + c]    	n_bulk_regions
 * 	R........region_be_matrix_									n_boundary_edges x n_boundary_regions
 *
 * Note that it holds:
 *
 * 	sv(q,c) = column sum of SV(q,c)
 *
 * Except for that, we also provide information on positive/negative flux and source:
 *
 * 	fp(q,c,r) = ( R' * EFP(q,c) )[r],	EFP(q,c)[e] = max{ 0, ( F(q,c) * solution + fv(q,c) )[e] }
 * 	fn(q,c,r) = ( R' * EFN(q,c) )[r],	EFN(q,c)[e] = min{ 0, ( F(q,c) * solution + fv(q,c) )[e] }
 * 	sp(q,c,r) = sum_{i in DOFS } max{ 0, ( S(q,c)[i,r] * solution[i] + SV(q,c)[i,r] ) }
 * 	sn(q,c,r) = sum_{i in DOFS } min{ 0, ( S(q,c)[i,r] * solution[i] + SV(q,c)[i,r] ) }
 *
 * where
 *
 * 	fp(q,c,r)...positive (outward) flux of q-th quantity's c-th component in region r
 * 	fn(q,c,r)...negative (inward) flux of q-th quantity's c-th component in region r
 * 	sp(q,c,r)...positive source (spring) of q-th quantity's c-th component in region r
 * 	sn(q,c,r)...negative source (sink) of q-th quantity's c-th component in region r
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

	/**
	 * Constructor.
	 * @param file_prefix  Prefix of output file name.
	 * @param mesh         Mesh.
	 * @param in_rec       Input record of balance.
	 */
	Balance(const std::string &file_prefix,
			const Mesh *mesh,
			const Input::Record &in_rec);

	~Balance();

	/// Setter for units of conserved quantities.
	void units(const UnitSI &unit) { units_ = unit; units_.undef(false); }

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
	 * Allocates matrices and vectors for balance.
	 * @param el_ds        Distribution of elements.
	 * @param el_4_loc     Local to global element numbering.
	 * @param n_loc_dofs            Number of solution dofs on the local process.
	 * @param max_dofs_per_boundary Number of dofs contributing to one boundary edge.
	 */
	void allocate(const Distribution *el_ds,
			const int *el_4_loc,
			unsigned int n_loc_dofs,
			unsigned int max_dofs_per_boundary);

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
	 * @param quantity_idx  Index of quantity.
	 * @param region_idx    Index of bulk region.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void add_mass_matrix_values(unsigned int quantity_idx,
			unsigned int region_idx,
			const std::vector<int> &dof_indices,
			const std::vector<double> &values);

	/**
	 * Adds elements into matrix for computing flux.
	 * @param quantity_idx  Index of quantity.
	 * @param elem_idx      Local index of boundary edge.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void add_flux_matrix_values(unsigned int quantity_idx,
			unsigned int elem_idx,
			const std::vector<int> &dof_indices,
			const std::vector<double> &values);

	/**
	 * Adds elements into matrix for computing source.
	 * @param quantity_idx  Index of quantity.
	 * @param region_idx    Index of bulk region.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void add_source_matrix_values(unsigned int quantity_idx,
			unsigned int region_idx,
			const std::vector<int> &dof_indices,
			const std::vector<double> &values);

	/**
	 * Adds element into vector for computing flux.
	 * @param quantity_idx  Index of quantity.
	 * @param elem_idx      Local index of boundary edge.
	 * @param value         Value to be added.
	 */
	void add_flux_vec_value(unsigned int quantity_idx,
			unsigned int elem_idx,
			double value);

	/**
	 * Adds elements into vector for computing source.
	 * @param quantity_idx  Index of quantity.
	 * @param region_idx    Index of bulk region.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void add_source_rhs_values(unsigned int quantity_idx,
			unsigned int region_idx,
			const std::vector<int> &dof_values,
			const std::vector<double> &values);

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
	 * It calculates the sum of flux and source over time interval.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 * @param dt            Actual time step.
	 */
	void calculate_cumulative_sources(unsigned int quantity_idx,
			const Vec &solution,
			double dt);

	/**
	 * Updates cumulative quantities for balance.
	 * This method can be called in substeps even if no output is generated.
	 * It calculates the sum of flux and source over time interval.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 * @param dt            Actual time step.
	 */
	void calculate_cumulative_fluxes(unsigned int quantity_idx,
			const Vec &solution,
			double dt);

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
	 * Calculates actual mass and save it to internal vector.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 */
	void calculate_mass(unsigned int quantity_idx,
			const Vec &solution)
	{
		calculate_mass(quantity_idx, solution, masses_[quantity_idx]);
	}

	/**
	 * Calculates actual flux.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 */
	void calculate_flux(unsigned int quantity_idx,
			const Vec &solution);

	/**
	 * Calculates actual source.
	 * @param quantity_idx  Index of quantity.
	 * @param solution      Solution vector.
	 */
	void calculate_source(unsigned int quantity_idx,
			const Vec &solution);

	/**
	* Adds provided values to the cumulative sources.
	* @param quantity_idx Index of quantity.
	* @param sources Sources per region.
	* @param dt Actual time step.
	*/
	void add_cumulative_source(unsigned int quantity_idx, double source);

	/// Perform output to file for given time instant.
	void output(double time);

private:
	/// Size of column in output (used if delimiter is space)
	static const unsigned int output_column_width = 20;

	/// Perform output in old format (for compatibility)
	void output_legacy(double time);

	/// Perform output in csv format
	void output_csv(double time, char delimiter, const std::string& comment_string, unsigned int repeat = 0);

	/// Return part of output represented by zero values. Count of zero values is given by cnt parameter.
	std::string csv_zero_vals(unsigned int cnt, char delimiter);

	/// Print output header
	void format_csv_output_header(char delimiter, const std::string& comment_string);

	/// Format string value of csv output. Wrap string into quotes and if delimiter is space, align text to column.
	std::string format_csv_val(std::string val, char delimiter, bool initial = false);

	/// Format double value of csv output. If delimiter is space, align text to column.
	std::string format_csv_val(double val, char delimiter, bool initial = false);


    /// Handle for file for output of balance and total fluxes over individual regions and region sets.
    ofstream output_;

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

    /// Vectors for calculation of source (n_bulk_regions).
    Vec *region_source_vec_;

    /**
     * Auxiliary matrix for transfer of quantities between boundary edges and regions
     * (n_boundary_edges x n_boundary_regions).
     */
    Mat region_be_matrix_;

    /// auxiliary vectors for summation of matrix columns
    Vec ones_, ones_be_;

    /// Number of boundary region for each local boundary edge.
    std::vector<unsigned int> be_regions_;

    /// Offset for local part of vector of boundary edges.
    int be_offset_;


    // Vectors storing mass and balances of fluxes and volumes.
    // substance, phase, region
    std::vector<std::vector<double> > fluxes_;
    std::vector<std::vector<double> > fluxes_in_;
    std::vector<std::vector<double> > fluxes_out_;
    std::vector<std::vector<double> > masses_;
    std::vector<std::vector<double> > sources_;
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

	/// initial time
	double initial_time_;

	/// time of last calculated balance
	double last_time_;

	/// true before calculating the mass at initial time, otherwise false
	bool initial_;

	/// if true then cumulative balance is computed
	bool cumulative_;

	/// true before allocating necessary internal structures (Petsc matrices etc.)
	bool allocation_done_;

	/// MPI rank.
	int rank_;

	/// hold count of line printed into output_
	unsigned int output_line_counter_;

};





#endif // BALANCE_HH_
