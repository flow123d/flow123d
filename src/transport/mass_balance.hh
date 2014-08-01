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

	Balance(const std::vector<std::string> &quantities,
			const std::vector<std::string> &components,
			const RegionDB *region_db,
			const bool cumulative,
			const std::string file);

	~Balance();

	inline bool cumulative() const
	{
		return cumulative_;
	}

	void allocate_matrices(unsigned int lsize);

	void start_mass_assembly(unsigned int quantity_idx,
			unsigned int component_idx)
	{
		MatZeroEntries(region_mass_matrix_[quantity_idx*components_.size()+component_idx]);
	}

	void start_flux_assembly(unsigned int quantity_idx,
			unsigned int component_idx)
	{
		MatZeroEntries(region_flux_matrix_[quantity_idx*components_.size()+component_idx]);
		MatZeroEntries(region_flux_rhs_[quantity_idx*components_.size()+component_idx]);
		VecZeroEntries(region_flux_vec_[quantity_idx*components_.size()+component_idx]);
	}

	void start_source_assembly(unsigned int quantity_idx,
			unsigned int component_idx)
	{
		MatZeroEntries(region_source_matrix_[quantity_idx*components_.size()+component_idx]);
		MatZeroEntries(region_source_rhs_[quantity_idx*components_.size()+component_idx]);
		VecZeroEntries(region_source_vec_[quantity_idx*components_.size()+component_idx]);
	}


	void set_mass_matrix_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	void set_flux_matrix_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	void set_source_matrix_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	void set_flux_rhs_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_indices,
			double *values);

	void set_source_rhs_values(unsigned int quantity_idx,
			unsigned int component_idx,
			unsigned int region_idx,
			int n_dofs,
			int *dof_values,
			double *values);

	void finish_mass_assembly(unsigned int quantity_idx,
			unsigned int component_idx)
	{
		MatAssemblyBegin(region_mass_matrix_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(region_mass_matrix_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
	}

	void finish_flux_assembly(unsigned int quantity_idx,
			unsigned int component_idx)
	{
		MatAssemblyBegin(region_flux_matrix_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(region_flux_matrix_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(region_flux_rhs_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(region_flux_rhs_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatMultTranspose(region_flux_rhs_[quantity_idx*components_.size()+component_idx], ones_, region_flux_vec_[quantity_idx*components_.size()+component_idx]);
	}

	void finish_source_assembly(unsigned int quantity_idx,
			unsigned int component_idx)
	{
		MatAssemblyBegin(region_source_matrix_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(region_source_matrix_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyBegin(region_source_rhs_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(region_source_rhs_[quantity_idx*components_.size()+component_idx], MAT_FINAL_ASSEMBLY);
		MatMultTranspose(region_source_rhs_[quantity_idx*components_.size()+component_idx], ones_, region_source_vec_[quantity_idx*components_.size()+component_idx]);
	}


	void calculate_cumulative(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution,
			double dt);

	void calculate_mass(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution);

	void calculate_flux(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution);

	void calculate_source(unsigned int quantity_idx,
			unsigned int component_idx,
			const Vec &solution);

	void output(double time);

private:


	void output_legacy(double time);


    /// Handle for file for output of balance and total fluxes over individual regions and region sets.
    ofstream output_;

    OutputFormat output_format_;

    const std::vector<std::string> quantities_;

    const std::vector<std::string> components_;

    const RegionDB &regions_;

    Mat *region_mass_matrix_;
    Mat *region_flux_matrix_;
    Mat *region_flux_rhs_;
    Mat *region_source_matrix_;
    Mat *region_source_rhs_;

    /// auxiliary vector for summation of matrix columns
    Vec ones_;
    Vec *region_flux_vec_;
    Vec *region_source_vec_;


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
