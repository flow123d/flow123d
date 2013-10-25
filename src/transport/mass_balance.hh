#ifndef MASS_BALANCE_HH_
#define MASS_BALANCE_HH_



/**
 * Interface class for equation which implements methods required for mass balance.
 */

class EquationForMassBalance {
public:

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
    virtual vector<string> &substance_names() = 0;

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

    MassBalance(EquationForMassBalance *eq, const char *output_f_name);

    ~MassBalance();

    /**
     * @brief Write computed fields to file.
     */
    void output(double time);

    /**
     * Calculate mass balance: flux through boundary, mass and volume sources
     */
    void calculate(double time);


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

};








#endif // MASS_BALANCE_HH_
