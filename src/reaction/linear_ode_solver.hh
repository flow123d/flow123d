#ifndef LINEAR_ODE_SOLVER_H_
#define LINEAR_ODE_SOLVER_H_

#include "armadillo"
#include "input/accessors.hh"

using namespace arma;

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
    
    void set_system_matrix(const mat &matrix);  ///< Sets the matrix of ODE system.
    void set_step(double step);                 ///< Sets the step of the numerical method.
    
    /// Updates solution of the system.
    virtual void update_solution(vec &init_vector, vec &output_vec) = 0;
    
    /// Updates solution of the system with different initial vectors.
    /**
     * Column initial and output vectors are grouped in the matrices.
     * Parameter @p mask can be used to skip some of the vectors.
     */
    virtual void update_solution(mat &init_vectors, mat &output_vecs, 
                                 const std::vector<unsigned int> &mask = std::vector<unsigned int>(0)) = 0;
                                 
protected:
    mat system_matrix_;     ///< the matrix of ODE system
    vec rhs_;               ///< the vector of RHS values (not used currently)
    double step_;           ///< the step of the numerical method
    bool step_changed_;     ///< flag is true if the step has been changed
};


/** @brief Template class of the linear ODE solver.
 * 
 * It provides a common method @p update_solution which can compute the same system of ODEs with 
 * different initial vectors at once.
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
    virtual void update_solution(mat &init_vectors, mat &output_vecs, 
                         const std::vector<unsigned int> &mask = std::vector<unsigned int>(0)) override;
    
private:
};

// template<class Method>
// LinearODESolver<Method>::LinearODESolver()
// {
// }


// template<class Method>
// void LinearODESolver<Method>::update_solution(vec& init_vector, vec& output_vec)
// {
//     static_cast<Method*>(this)->update_solution(init_vector,output_vec);
//     step_changed_ = false;
// }

template<class Method>
void LinearODESolver<Method>::update_solution(mat& init_vectors, mat& output_vecs, const std::vector< unsigned int > &mask)
{  
}





/*********************** ANALYTIC SOLUTION ************************/
/** @brief This class implements the analytic solution of a system of linear ODEs with constant matrix.
 *
 * The analytic solution can be obtained in special case when decrease of one quantity is a fraction
 * of increase of other quantity. The fractions are passed through the bifurcation matrix.
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
    
    void update_solution(vec &init_vector, vec &output_vec) override;
//     void update_solution(mat &init_vectors, mat &output_vecs, 
//                          const std::vector<unsigned int> &mask = std::vector<unsigned int>(0)) override;
    
    void set_bifurcation_matrix(const mat &bifurcation);
    
protected:
    /**
     *   Computes the standard fundamental matrix.
     */
    void compute_matrix();
    
    /// The bifurcation matrix. 
    mat bifurcation_matrix_;
    
    /// The solution is computed only by matrix multiplication (standard fundamental matrix).
    mat solution_matrix_;
};

#endif // LINEAR_ODE_SOLVER_H_