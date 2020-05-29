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
 * @file    hm_iterative.hh
 * @brief   
 * @author  Jan Stebel
 */

#ifndef HM_ITERATIVE_HH_
#define HM_ITERATIVE_HH_

#include <memory>
#include <string>
#include <vector>
#include "input/input_type_forward.hh"
#include "input/accessors_forward.hh"
#include "coupling/equation.hh"
#include "flow/darcy_flow_interface.hh"
#include "mechanics/elasticity.hh"

class Mesh;
class FieldCommon;
class RichardsLMH;


class IterativeCoupling {
public:

    static const Input::Type::Record &record_template() {
        return it::Record("Coupling_Iterative_AUX",
                "Record with data for iterative coupling.\n")
            .declare_key( "max_it", it::Integer(0), it::Default("100"),
                    "Maximal count of HM iterations." )
            .declare_key( "min_it", it::Integer(0), it::Default("1"),
                    "Minimal count of HM iterations." )
            .declare_key( "a_tol", it::Double(0), it::Default("0"),
                    "Absolute tolerance for difference in HM iteration." )
            .declare_key( "r_tol", it::Double(0), it::Default("1e-7"),
                    "Relative tolerance for difference in HM iteration." )
            .close();
    }

    IterativeCoupling(Input::Record in_record)
    : it(0)
    {
        min_it_ = in_record.val<unsigned int>("min_it");
        max_it_ = in_record.val<unsigned int>("max_it");
        a_tol_ = in_record.val<double>("a_tol");
        r_tol_ = in_record.val<double>("r_tol");
    }

    void solve_step()
    {
        it = 0;
        double abs_error = std::numeric_limits<double>::max();
        double rel_error = std::numeric_limits<double>::max();

        while ( it < min_it_ || (abs_error > a_tol_ && rel_error > r_tol_ && it < max_it_) )
        {
            it++;
            solve_iteration();
            compute_iteration_error(abs_error, rel_error);
            update_after_iteration();
        }
        update_after_converged();
    }

    unsigned int iteration()
    { return it; }

protected:

    /// Solve equations and update data (fields).
    virtual void solve_iteration() = 0;

    /// Compute absolute and relative error in the solution.
    virtual void compute_iteration_error(double &abs_error, double &rel_error) = 0;

    /// Save data (e.g. solution fields) for the next iteration.
    virtual void update_after_iteration() = 0;

    /// Save data after iterations have finished.
    virtual void update_after_converged() = 0;


    /// Minimal number of iterations to perform.
    unsigned int min_it_;
    
    /// Maximal number of iterations.
    unsigned int max_it_;
    
    /// Absolute tolerance for difference between two succeeding iterations.
    double a_tol_;
    
    /// Relative tolerance for difference between two succeeding iterations.
    double r_tol_;

private:

    /// Iteration index.
    unsigned int it;

};


/**
 * @brief Class for solution of fully coupled flow and mechanics using fixed-stress iterative splitting.
 * 
 * Flow and mechanics are solved separately and within each iteration the coupling terms are updated.
 * Here we use the fixed-stress splitting [see Mikelic&Wheeler, Comput. Geosci. 17(3), 2013] which uses
 * a tuning parameter "beta" to speed up the convergence.
 */
class HM_Iterative : public DarcyFlowInterface, public IterativeCoupling {
public:
    
    class EqData : public FieldSet
    {
    public:
        EqData();
        
        void initialize(Mesh &mesh);
        
        Field<3, FieldValue<3>::Scalar> alpha;   ///< Biot coefficient.
        Field<3, FieldValue<3>::Scalar> density; ///< Density of fluid.
        Field<3, FieldValue<3>::Scalar> gravity; ///< Standard gravity.
        Field<3, FieldValue<3>::Scalar> beta;
        
        /// Potential -alpha*pressure whose gradient is passed to mechanics as additional load.
        Field<3, FieldValue<3>::Scalar> pressure_potential;
        Field<3, FieldValue<3>::Scalar> flow_source;
        
        /// FieldFE for pressure_potential field.
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > potential_ptr_;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > beta_ptr_;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > flow_source_ptr_;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > old_pressure_ptr_;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > old_iter_pressure_ptr_;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > div_u_ptr_;
        std::shared_ptr<FieldFE<3, FieldValue<3>::Scalar> > old_div_u_ptr_;
    };
    
    /// Define input record.
    static const Input::Type::Record & get_input_type();

    HM_Iterative(Mesh &mesh, Input::Record in_record);
    void initialize() override;
    void zero_time_step() override;
    void update_solution() override;
    double last_t() override;
    ~HM_Iterative();

private:
    
    void update_potential();
    
    void update_flow_fields();

    void solve_iteration() override;

    void update_after_iteration() override;

    void update_after_converged() override;
    
    void compute_iteration_error(double &abs_error, double &rel_error) override;
    
    static const int registrar;

    /// steady or unsteady water flow simulator based on MH scheme
    std::shared_ptr<RichardsLMH> flow_;

    /// solute transport with chemistry through operator splitting
    std::shared_ptr<Elasticity> mechanics_;
    
    EqData data_;

    /// Tuning parameter for iterative splitting.
    double beta_;
    
};

#endif /* HC_EXPLICIT_SEQUENTIAL_HH_ */
