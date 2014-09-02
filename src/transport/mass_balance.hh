#ifndef MASS_BALANCE_HH_
#define MASS_BALANCE_HH_


#include "la/distribution.hh"
#include "transport/substance.hh"
#include "petscmat.h"

/**
 * Interface class for equation which implements methods required for mass balance.
 */

class EquationForMassBalance {
public:

	// TODO: Think if we really need TimeIntegrationScheme.
	// Currently the cummulative quantities are calculated by the method calculate()
	// in the same way for both explicit and implicit methods.
    enum TimeIntegrationScheme {
    	none,
    	explicit_euler,
    	implicit_euler,
    	crank_nicholson
    };

	virtual ~EquationForMassBalance() {};

    /// Returns number of trnasported substances.
    virtual unsigned int n_substances() = 0;

    /// Returns reference to the vector of substnace names.
    virtual SubstanceList &substances() = 0;

    /// Returns the time integration scheme of the equation.
	virtual TimeIntegrationScheme time_scheme() = 0;


protected:

    /**
     * Calculates the total flux through boundaries of all regions, and additionally positive and negative fluxes.
     * The actual calculation depends on the numerical scheme, so each descendant of TransportBase implements this method.
     * @param bcd_balance       bcd_balance[i][j] is the calculated total flux
     *                          of @p ith substance through boundary of @p jth region.
     * @param bcd_plus_balance  bcd_plus_balance[i][j] is the total positive flux
     *                          of @p ith substance through boundary of @p jth region.
     * @param bcd_minus_balance bcd_minus_balance[i][j] is the total negative flux
     *                          of @p ith substance through boundary of @p jth region.
     */
    virtual void calc_fluxes(vector<vector<double> > &bcd_balance, vector<vector<double> > &bcd_plus_balance, vector<vector<double> > &bcd_minus_balance) = 0;

    /**
     * Calculates the substance mass and sources on all regions.
     * The actual calculation depends on the numerical scheme, so each descendant of TransportBase implements this method.
     * @param mass        mass[i][j] is the calculated mass of @p ith
     *                    substance on @p jth region.
     * @param src_balance src_balance[i][j] is the source mass
     *                    of @p ith substance on @p jth region.
     */
    virtual void calc_elem_sources(vector<vector<double> > &mass, vector<vector<double> > &src_balance) = 0;

    /// Returns the region database.
    virtual const RegionDB *region_db() = 0;

    friend class MassBalance;

};


/**
 * @brief Class for calculation and writing the balance of mass, volume sources and fluxes.
 *
 * At each time instant we calculate:
 * - Flux rate through all boundary regions
 * - Volume sources rate on all bulk regions
 * - Mass on all bulk regions
 *
 * In addition, the quantities are summed over all regions and integrated in time,
 * so that the cumulative (integrated) quantities should satisfy:
 *
 *     INITIAL_MASS + INTEGRATED_SOURCES - INTEGRATED_FLUXES = CURRENT_MASS.
 *
 * All quantities are written to the file "mass_balance.txt" in the output directory.
 */
class MassBalance {
public:

    MassBalance(EquationForMassBalance *eq, const Input::Record &in_rec);

    ~MassBalance();

    /**
     * @brief Write computed fields to file.
     */
    void output(double time);

    /**
     * Calculate mass balance: flux through boundary, mass and volume sources
     */
    void calculate(double time);


    static Input::Type::Record input_type;


protected:

    /// Pointer to the class which implements calculation of mass, sources and fluxes.
    EquationForMassBalance *equation_;

    /// Handle for output file for output of balance and total fluxes over individual regions and region sets.
    FILE *balance_output_file;

    // Vectors storing mass and balances of fluxes and volumes.
    vector<vector<double> > bcd_balance;
    vector<vector<double> > bcd_plus_balance;
    vector<vector<double> > bcd_minus_balance;
    vector<vector<double> > mass;
    vector<vector<double> > src_balance;

	vector<double> bcd_total_balance;
	vector<double> bcd_total_inflow;
	vector<double> bcd_total_outflow;
	vector<double> mass_total;
	vector<double> src_total_balance;

	vector<double> initial_mass;
	vector<double> integrated_sources;
	vector<double> integrated_fluxes;

	/// initial time
	double initial_time;

	/// time of last calculated balance
	double last_time;

	/// true before calculating the mass at initial time, otherwise false
	bool initial;

	/// if true then cumulative balance is computed
	bool cumulative;

};





/**
 * New design of balance class - serves as storage and writer.
 * Equations themselves call methods of Balance that add/modify mass, source and flux
 * of various quantities and generate output.
 *
 * Another new feature is that quantities can have more components,
 * e.g. velocity.{x,y,z} or concentration.{mobile,immobile}.
 * Each quantity has an implicit component "default".
 *
 * The mass, flux and source are calculated as follows:
 *
 * 	m(q,c,r) = ( M'(q,c) * solution )[r]
 * 	f(q,c,r) = ( R' * ( F(q,c) * solution + fv(q,c) ) )[r]
 * 	s(q,c,r) = ( S'(q,c) * solution + sv(q,c) )[r]
 *
 * where
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
 *
 *
 * Output values (if not relevant, zero is supplied):
 *
 * #time 	bulk_region     	quantity 	component 	0    	0       	0        	mass 	source 	0               	0                 	0
 * #time 	boundary_region 	quantity 	component 	flux 	flux_in 	flux_out 	0    	0      	0               	0                 	0
 * #time 	ALL             	quantity 	component 	flux 	flux_in 	flux_out 	mass 	source 	integrated_flux 	integrated_source 	error
 *
 * error = current_mass - (initial_mass + integrated_source - integrated_flux)
 *
 */
class Balance {
public:

	enum OutputFormat
	{
		legacy,
		csv,
		gnuplot
	};

	static Input::Type::Selection format_selection_input_type;

	/**
	 * Constructor.
	 * @param quantities   Names of conserved quantities.
	 * @param components   Names of components of each quantity.
	 * @param elem_regions Vector of region numbers for each boundary edge.
	 * @param region_db    Region database.
	 * @param cumulative   If true, cumulative sums will be calculated.
	 * @param file         Name of output file.
	 */
	Balance(const std::vector<std::string> &quantities,
			const std::vector<std::string> &components,
			const std::vector<unsigned int> &elem_regions,
			const RegionDB *region_db,
			const bool cumulative,
			const std::string file);

	~Balance();

	/// Getter for cumulative_.
	inline bool cumulative() const { return cumulative_; }

	/**
	 * Allocates matrices and vectors for balance.
	 * @param n_loc_dofs            Number of solution dofs on the local process.
	 * @param max_dofs_per_boundary Number of dofs contributing to one boundary edge.
	 */
	void allocate_matrices(unsigned int n_loc_dofs, unsigned int max_dofs_per_boundary);

	/**
	 * This method must be called before assembling the matrix for computing mass.
	 * It actually erases the matrix.
	 */
	void start_mass_assembly(unsigned int quantity_idx, unsigned int component_idx);

	/**
	 * This method must be called before assembling the matrix and vector for fluxes.
	 * It actually erases the matrix and vector.
	 */
	void start_flux_assembly(unsigned int quantity_idx, unsigned int component_idx);

	/**
	 * This method must be called before assembling the matrix and vectors for sources.
	 * It actually erases the matrix and vectors.
	 */
	void start_source_assembly(unsigned int quantity_idx, unsigned int component_idx);

	/**
	 * Adds elements into matrix for computing mass.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param region_idx    Index of bulk region.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void set_mass_matrix_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	/**
	 * Adds elements into matrix for computing flux.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param elem_idx      Local index of boundary edge.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void set_flux_matrix_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int elem_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	/**
	 * Adds elements into matrix for computing source.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param region_idx    Index of bulk region.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void set_source_matrix_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	/**
	 * Adds element into vector for computing flux.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param elem_idx      Local index of boundary edge.
	 * @param value         Value to be added.
	 */
	void set_flux_vec_value(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int elem_idx,
			double value);

	/**
	 * Adds elements into vector for computing source.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param region_idx    Index of bulk region.
	 * @param n_dofs        Number of dofs to be added.
	 * @param dof_indices   Dof indices to be added.
	 * @param values        Values to be added.
	 */
	void set_source_rhs_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_values,
			double *values);

	/// This method must be called after assembling the matrix for computing mass.
	void finish_mass_assembly(unsigned int quantity_idx, unsigned int component_idx);

	/// This method must be called after assembling the matrix and vector for computing flux.
	void finish_flux_assembly(unsigned int quantity_idx, unsigned int component_idx);

	/// This method must be called after assembling the matrix and vectors for computing source.
	void finish_source_assembly(unsigned int quantity_idx, unsigned int component_idx);

	/**
	 * Updates cumulative quantities for balance.
	 * This method can be called in substeps even if no output is generated.
	 * It calculates the sum of flux and source over time interval.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param solution      Solution vector.
	 * @param dt            Actual time step.
	 */
	void calculate_cumulative_sources(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution,
			double dt);

	/**
	 * Updates cumulative quantities for balance.
	 * This method can be called in substeps even if no output is generated.
	 * It calculates the sum of flux and source over time interval.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param solution      Solution vector.
	 * @param dt            Actual time step.
	 */
	void calculate_cumulative_fluxes(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution,
			double dt);

	/**
	 * Calculates actual mass.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param solution      Solution vector.
	 */
	void calculate_mass(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution);

	/**
	 * Calculates actual flux.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param solution      Solution vector.
	 */
	void calculate_flux(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution);

	/**
	 * Calculates actual source.
	 * @param quantity_idx  Index of quantity.
	 * @param component_idx Index of component.
	 * @param solution      Solution vector.
	 */
	void calculate_source(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution);

	/// Perform output to file for given time instant.
	void output(double time);

private:

	/// Perform output in old format (for compatibility)
	void output_legacy(double time);


    /// Handle for file for output of balance and total fluxes over individual regions and region sets.
    ofstream output_;

    /// Format of output file.
    OutputFormat output_format_;

    /// Names of conserved quantities.
    const std::vector<std::string> quantities_;

    /// Names of components.
    const std::vector<std::string> components_;

    /// Database of bulk and boundary regions.
    const RegionDB &regions_;


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
    const std::vector<unsigned int> be_regions_;

    /// Offset for local part of vector of boundary edges.
    int be_offset_;


    // Vectors storing mass and balances of fluxes and volumes.
    // substance, phase, region
    std::vector<std::vector<std::vector<double> > > fluxes_;
    std::vector<std::vector<std::vector<double> > > fluxes_in_;
    std::vector<std::vector<std::vector<double> > > fluxes_out_;
    std::vector<std::vector<std::vector<double> > > masses_;
    std::vector<std::vector<std::vector<double> > > sources_;
    std::vector<std::vector<std::vector<double> > > sources_in_;
    std::vector<std::vector<std::vector<double> > > sources_out_;

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

	/// initial time
	double initial_time_;

	/// time of last calculated balance
	double last_time_;

	/// true before calculating the mass at initial time, otherwise false
	bool initial_;

	/// if true then cumulative balance is computed
	bool cumulative_;

	int rank_;

};





#endif // MASS_BALANCE_HH_
