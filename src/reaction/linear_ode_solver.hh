#ifndef LINEAR_ODE_SOLVER_H_
#define LINEAR_ODE_SOLVER_H_

#include "armadillo"
#include "input/accessors.hh"


/// @brief Base class for linear ODE solver.
/** This class represents an interface to a solver of a system of linear ordinary differential 
 *  equations with constant coefficients.
 */
class LinearODESolverBase
{
public:
    /**
     * Abstract record for the linear ODE solver.
     */
    static Input::Type::AbstractRecord input_type;
    
    LinearODESolverBase();
    virtual ~LinearODESolverBase();
    
    void set_system_matrix(const arma::mat &matrix);  ///< Sets the matrix of ODE system.
    void set_step(double step);                 ///< Sets the step of the numerical method.
    
    /// Updates solution of the ODEs system.
    virtual void update_solution(arma::vec &init_vector, arma::vec &output_vec) = 0;
    
    /// Updates solution of the system with different initial vectors.
    /**
     * Column initial and output vectors are grouped in the matrices.
     * Parameter @p mask can be used to skip some of the vectors.
     */
    virtual void update_solution(arma::mat &init_vectors, arma::mat &output_vecs, 
                                 const std::vector<unsigned int> &mask = std::vector<unsigned int>(0)) = 0;
                                 
protected:
    arma::mat system_matrix_;     ///< the matrix of ODE system
    arma::vec rhs_;               ///< the vector of RHS values (not used currently)
    double step_;           ///< the step of the numerical method
    bool step_changed_;     ///< flag is true if the step has been changed
};


/** @brief Template class of the linear ODE solver.
 * 
 * It provides a common method @p update_solution which can compute the same system of ODEs with 
 * different initial vectors at once.
 * 
 * This class represents the Curiously Recurring Template Pattern (CRTP). Therefore, the method update_solution
 * of the template called inside @p update_solution using static_cast is not a virtual one.
 * 
 */
template<class Method>
class LinearODESolver : public LinearODESolverBase
{   
public:
    LinearODESolver(){};
    virtual ~LinearODESolver(){};
    
    /// Updates solution of the system with different initial vectors.
    /**
     * Column initial and output vectors are grouped in the matrices.
     * Parameter @p mask can be used to skip some of the vectors.
     */
    virtual void update_solution(arma::mat &init_vectors, arma::mat &output_vecs, 
                         const std::vector<unsigned int> &mask = std::vector<unsigned int>(0)) override;
    
private:
};

template<class Method>
void LinearODESolver<Method>::update_solution(arma::mat& init_vectors, arma::mat& output_vecs, const std::vector< unsigned int > &mask)
{  
    ASSERT(0,"Not implemented yet.");
    ASSERT_EQUAL(init_vectors.n_cols, output_vecs.n_cols);
    ASSERT_EQUAL(init_vectors.n_rows, output_vecs.n_rows);
    
    for(unsigned int j=0; j < init_vectors.n_cols; j++)
    {
        //static_cast<Method*>(this)->update_solution(init_vectors.col(j), output_vecs.col(j));
    }
}





/*********************** ANALYTIC SOLUTION ************************/
/** @brief This class implements the analytic solution of a system of linear ODEs with constant matrix.
 *
 * The analytic solution can be obtained in the special case when decrease of one quantity is a fraction
 * of increase of other quantity. The fractions are derived from the matrix.
 * 
 * It is used in first order reactions and decays.
 */
class LinearODEAnalytic : public LinearODESolver<LinearODEAnalytic>
{
public:
    /**
     * Input record for class LinearODE_analytic.
     */
    static Input::Type::Record input_type;
    
    ///Default constructor is possible because the input record is not needed.
    LinearODEAnalytic(){};
    
    /// Constructor from the input data.
    LinearODEAnalytic(Input::Record in_rec);

    /// Destructor.
    ~LinearODEAnalytic(void);
    
    void update_solution(arma::vec &init_vector, arma::vec &output_vec) override;
    
protected:
    /**
     *   Computes the standard fundamental matrix.
     */
    void compute_matrix();
    
    /// The solution is computed only by a matrix multiplication (standard fundamental matrix).
    arma::mat solution_matrix_;
};

#endif // LINEAR_ODE_SOLVER_H_