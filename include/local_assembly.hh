/*
 * local_assembly.hh
 *
 *  Created on: Nov 1, 2011
 *      Author: jb
 */

#ifndef LOCAL_ASSEMBLY_HH_
#define LOCAL_ASSEMBLY_HH_

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/lapack_full_matrix.h>

#include <base/tensor_base.h>
#include <base/quadrature_lib.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_dg_vector.h>
#include <fe/fe_raviart_thomas.h>
#include <base/polynomials_raviart_thomas.h>
#include <fe/fe_face.h>
#include <fe/fe_dgp_monomial.h>

#include <fe/fe_values.h>

#include <spatial_functions.hh>

using namespace dealii;

/**
 * Keep solution vector, z-coordinates of its support points and
 * precomputed nonlinear functions.
 */
template <int dim>
class Solution {
public:
    Solution(RichardsData<dim> *rd)
    : richards_data(rd) {}

    /**
     * Compute z-coordinates of support points of solution values.
     */
    void reinit(DoFHandler<dim> &dh, bool save_old_solution) {
        z_coord.reinit(dh.n_dofs());
        phead.reinit(dh.n_dofs());
        old_phead.reinit(dh.n_dofs());

        lambda.reinit(dh.n_dofs());
        lambda_diff.reinit(dh.n_dofs());
        lambda_new.reinit(dh.n_dofs());

        conductivity_new.reinit(dh.n_dofs());
        conductivity_new_diff.reinit(dh.n_dofs());
        conductivity_old.reinit(dh.n_dofs());

        saturation.reinit(dh.n_dofs());
        capacity.reinit(dh.n_dofs());
        old_saturation.reinit(dh.n_dofs());
        residual.reinit(dh.n_dofs());
        head_diff.reinit(dh.n_dofs());
        head_second_diff.reinit(dh.n_dofs());

        out_aux.reinit(dh.n_dofs());
        el_sum.reinit(dh.n_dofs());
        //sat_diff.reinit(dh.n_dofs());
        //time_lambda.reinit(dh.n_dofs());
        //sat_diff = 0.0;

        typename DoFHandler<dim>::active_cell_iterator
            cell = dh.begin_active(),
            endc = dh.end();
        std::vector<unsigned int> face_dof_indices(1);

        for (; cell != endc; ++cell) {
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
                cell->face(face_no)->get_dof_indices(face_dof_indices);
                z_coord(face_dof_indices[0]) = cell->face(face_no)->barycenter() [dim-1];
            }
        }

        head_diff= 0;
        dt_ =1;
        cond_type = trapezoid;
        //cond_type = mid_point_quad;
        //cond_type = mid_point;
    }

    void revert() {
        phead=old_phead;
        saturation=old_saturation;
    }
    /**
     * Save old values and update the new.
     */
    void timestep_update(double dt) {
        double p;
        dt_last = dt_;
        dt_=dt;

        for(unsigned int i=0;i<phead.size();i++) {
            head_diff(i) = (phead(i) - old_phead(i)) /dt_;
            old_phead(i) = phead(i);
            old_saturation(i) = saturation(i);
            conductivity_old(i) = conductivity_new(i);

        }
        last_lmb_change=false;
    }

    void update(double s_param) {
        double p, p_scaled;
        double cap_crit;
        double h_crit_1 = richards_data->cap_arg_max_/2 ;
        double h_crit_0 = richards_data->cap_arg_max_/2 ;
        richards_data->fq(h_crit_0,cap_crit);
        cap_crit = cap_crit / richards_data->fk(h_crit_0);
        double h_crit = richards_data->cap_arg_max_/4;



        int  n_lmb_change = 0;
        for(unsigned int i=0;i<phead.size();i++) {
            p = pressure(i);
            saturation(i)=richards_data->fq(p,capacity(i));

            // lambda
            /*// linear lambda
            if (p>0) {
                lambda(i)=1;
                lambda_diff(i)=0;
                lambda_s_diff(i)=0;
            } else {
                p_scaled = p / (richards_data->cap_arg_max_/2.0);
                if (p_scaled > 1) {
                    lambda(i) = 1 - s_param * 0.5;
                    lambda_diff(i)=0;
                    lambda_s_diff(i)=-0.5;
                } else {
                    lambda(i) = 1 - s_param * 0.5 * p_scaled;
                    lambda_diff(i)= - s_param * 0.5/ (richards_data->cap_arg_max_/2.0);
                    lambda_s_diff(i)= - 0.5 * p_scaled;
                }
            } */
            // constant lambda
            //lambda(i)=0.5;
            //lambda_diff(i)=0;
            //lambda_s_diff(i)=0;


            // linear time lambda betwean old saturation point and actual small capacity ( betwean these points we have problem)
            double hd = pressure(i) - h_crit;
            double interpol= hd/(hd-old_pressure(i));
            interpol = max(interpol, 0.0);
            interpol = min(interpol, 1.0);

            lambda(i) = 0.5+0.5*interpol;
            //if (old_pressure(i) > 0.0) lambda(i) =1;
            //else lambda(i)=0.5;

/*

            if (p> h_crit_0) {
                time_lambda(i)=0;
                //lambda_diff(i)=0;
                //lambda_s_diff(i)=0;
            } else {
                p_scaled = (p-h_crit_0) / (h_crit_1- h_crit_0);
                if (p_scaled > 1) {
                    time_lambda(i) = 1;
                    //lambda_diff(i)=0;
                    //lambda_s_diff(i)=-0.5;
                } else {
                    time_lambda(i) = 1.0/(1.0+ pow((1-p_scaled),1.0)*pow(p_scaled,-1.0)); // second exp. is near zero
                    //time_lambda(i) = time_lambda(i) / (2 - time_lambda(i));
                    //lambda_diff(i)= - s_param * 0.5/ (richards_data->cap_arg_max_/2.0);
                    //lambda_s_diff(i)= - 0.5 * p_scaled;
                }
            }

//            lambda(i) = 1.0/(1.0+time_lambda(i));
            time_lambda(i) =1;*/

            // lambda cap/cond
/*
            if (cond_type != mid_point_quad) {
                s_param = 1;
                double c_diff;
                double cond = richards_data->fk(p, c_diff);
                if (p > 0) {
                    lambda(i) = 1;
                    lambda_diff(i) = 0;
                    //lambda_s_diff(i) = 0;
                } else {
                    p_scaled = p / h_crit_0;
                    if (p_scaled > 1) {
                        lambda(i) = 1 - s_param * 0.5;
                        lambda_diff(i) = 0;
                      //  lambda_s_diff(i) = -0.5;
                    } else {
                        lambda(i) = 1 - s_param * 0.5 * (capacity(i) / cond) / dt_; // / (cap_crit);
                        lambda_diff(i) = -s_param * 0.5 * (richards_data->fqq(p) / cond - capacity(i) / cond / cond * c_diff) / dt_; // / (cap_crit);
                        //lambda_s_diff(i) = 0.0; //- 0.5 * (capacity(i)/cond) / (cap_over_k__crit);
                    }
                }
                if (lambda(i) < 0.5) {
                    lambda(i) = 0.5;
                    lambda_diff(i) = 0;
                }
            }*/
            //lambda(i) = lambda_new(i);
            /*
            double c_diff;
            double cond=richards_data->fk(p, c_diff);
            if (p>0) {
                lambda(i)=1;
                lambda_diff(i)=0;
                lambda_s_diff(i)=0;
            } else {
                p_scaled = p / h_crit;
                if (p_scaled > 1) {
                    lambda(i) = 1 - s_param * 0.5;
                    lambda_diff(i)=0;
                    lambda_s_diff(i)=-0.5;
                } else {
                    lambda(i) =  (1-0.5*s_param) + 0.5*s_param/
                            (1.0+square(
                                       1.0/(1.0- square(p_scaled)) - 1.0
                                      )
                             );
                    lambda_diff(i)= 0.0;//- s_param * 0.5 * (richards_data->fqq(p) /cond - capacity(i)/cond /cond * c_diff);
                    lambda_s_diff(i)= -0.5 + 0.5/
                            (1.0+square(
                                       1.0/(1.0- square(p_scaled)) - 1.0
                                      )
                             );
                }
            }*/

            //p=lambda(i)*p + (1-lambda(i))*old_pressure(i);

            if (cond_type == mid_point) {
                p=lambda(i)*pressure(i) + (1-lambda(i))*old_pressure(i);
            } else {
                p=pressure(i);
            }
            conductivity_new(i)=richards_data->fk(p, conductivity_new_diff(i));
            //if (cond_type == mid_point_quad) {
            //    conductivity_old(i)=richards_data->fk(pressure(i));
           // }

        }
    }

    void lambda_update() {
        for(unsigned int i=0;i<phead.size();i++) {
            lambda(i) =1;
          /*
          double div =  capacity(i) * el_sum(i) * (phead(i) - old_phead(i));

          if (fabs(div) < 1e-24) {
              div = (div >= 0 ? 1 : -1) * 1.0e-24;
          }

          {
              lambda(i)=0.75 - 0.25*lambda_new(i) / div;
              lambda(i) = max(0.5, lambda(i));
              lambda(i) = min(1.0, lambda(i));
          }*/

        }
    }

    inline double pressure(unsigned int i)
        { return richards_data->pressure(phead(i), z_coord(i)); }

    inline double old_pressure(unsigned int i)
        { return richards_data->pressure(old_phead(i), z_coord(i)); }

    inline double square(const double x) const
        {return x*x;}

    // smooth approximation of max(a,b)
    inline double my_max(double a, double b) {
        return 0.5*(sqrt(a*a+b*b+(a-b)*(a-b)) +a+b);
    }

    // derivative of approximative maximum by a, this is same as derivative by b, but you has to change variables
    // i.e. my_max_diff_a(a,b) is derivative by a, my_max_diff_a(b,a) is derivative by b,
    inline double my_max_diff_a(double a, double b) {
        return 0.5*0.5/sqrt((a*a+b*b+(a-b)*(a-b)))*(2*a +2*a -2*b)  + 0.5;
    }

    inline double conduct_new(unsigned int i) {
        if (cond_type == trapezoid) {
            return conductivity_new(i);
        } else {
            return conductivity_midpoint(i);
        }
    }

    inline double conduct_old(unsigned int i) {
        if (cond_type == trapezoid) {
            return conductivity_old(i);
        } else {
            return conductivity_midpoint(i);
        }
    }

    inline double conduct_new_diff(unsigned int i) {
        return conductivity_new_diff(i);
    }

    void compute_head_second_diff() {
        for(unsigned int i=0;i<phead.size();i++) {
            head_second_diff(i) = ((phead(i) - old_phead(i)) /dt_ - head_diff(i) )/ (dt_ + dt_last) * (1.0/6 *dt_*dt_);
        }
    }

    double conductivity_new_midpoint(unsigned int i) {
        return richards_data->fk(pressure(i));
    }

    Vector<double> z_coord;
    Vector<double> phead;
    Vector<double> old_phead;

    Vector<double> lambda;
    Vector<double> lambda_diff;
    Vector<double> lambda_new;

    Vector<double> saturation;
    Vector<double> capacity;
    Vector<double> old_saturation;
    Vector<double> residual;
    Vector<double> head_diff;
    Vector<double> head_second_diff;
    Vector<double> out_aux;
    Vector<double> el_sum;
    //Vector<double> sat_diff;
    //Vector<double> time_lambda;


    RichardsData<dim> *richards_data;
    double dt_, dt_last;
    bool last_lmb_change;


    typedef enum { trapezoid=0, mid_point_approx=1, mid_point=2, mid_point_quad=3 } CondType;
    CondType cond_type;

private:
    inline double conductivity_midpoint(unsigned int i) {
        if (cond_type == mid_point_approx) {
            return lambda(i)*conductivity_new(i) + (1-lambda(i))* conductivity_old(i);
        } else {
            return conductivity_new(i);
        }
    }
    Vector<double> conductivity_new;
    Vector<double> conductivity_old;
    Vector<double> conductivity_new_diff;


};

template <int dim>
class LocalAssembly {
public:
    //typedef FE_DGPMonomial<dim> FE_PostPressure;
    typedef FE_DGQ<dim> FE_PostPressure;
        //FE_Q<dim>;
        //FE_DGP<dim>;
    const static unsigned int post_order=1; // Q1=1,DGQ1=1,DGP2=2

    enum block_index_names {
        velocity_bl=0,
        half_flux_bl=1,
        pressure_bl=2,
        saturation_bl=3,
        p_error_bl=4,
        q_error_bl=5,
        post_p_bl=6,
        post_old_p_bl=7,
        post_residual=8,
        post_lambda=9,
        post_aux=10
    };

    LocalAssembly(unsigned int order);

    void set_data(Solution<dim> *sol)
    {
        solution = sol;
    }

    void set_dt(const double in_dt, const double in_t) {
        dt = in_dt;
        time = in_t;
    }

    void reinit(typename DoFHandler<dim>::active_cell_iterator cell);
    void output_evaluate(double & bc_flux_total,double &volume_total);
    void compute_add_variation(double &flux_var, double &head_var);
    void set_lambda();
    double get_p_error() {
        return  output_vector( out_idx[p_error_bl][0] ); }
    double get_q_error() { return  output_vector( out_idx[q_error_bl][0] ); }

    void make_A();
    void make_C();
    void apply_bc(typename DoFHandler<dim>::active_cell_iterator cell, std::map<unsigned int, double> &boundary_values);
    FullMatrix<double> &get_matrix(bool symmetric);
    Vector<double> &get_func();
    void compute_steady_flux();
    Vector<double> &get_s_diff();
    void update_lambda();


    Vector<double> &get_output_vector() {return output_vector;}

    double compute_p_error();
    double compute_q_error();

    void make_interpolation(FE_DGQ<dim> &fe);
    void make_interpolation(FE_DGPMonomial<dim> &fe);
    void out_interpolate(Vector<double> &trace_phead, unsigned int out_block);
    //void out_interpolate_pressure(FE_DGP<dim> &fe, Vector<double> &trace_phead) {};
    //void out_interpolate_pressure(FE_DGPMonomial<dim> &fe, Vector<double> &trace_phead);
    //void out_interpolate_pressure(FE_Q<dim> &fe, Vector<double> &trace_phead) {};

    //RichardsData<dim> *richards_data;
    Solution<dim> *solution;

private:
    typename DoFHandler<dim>::active_cell_iterator dh_cell;

    QGauss<dim> quadrature_formula;
    QGauss<dim> error_quadrature;
    QGauss < dim - 1 > face_quadrature_formula;

    FE_DGVector<PolynomialsRaviartThomas<dim>, dim> velocity_fe;
    FE_FaceQ<dim> trace_fe;
    FE_DGQ<dim> pressure_fe;
    FE_PostPressure post_pressure_fe;
public:
    FESystem<dim> output_fe;
private:
    //FEValues<dim> pressure_fe_values;
    FEValues<dim> velocity_fe_values;
    FEValues<dim> error_fe_values;

    FEFaceValues<dim> trace_fe_face_values;
    FEFaceValues<dim> velocity_fe_face_values;

    const FEValuesExtractors::Vector vec_extr;
    const FEValuesExtractors::Scalar scal_extr;


    FullMatrix<double>  local_matrix_00,    // block A
                        local_matrix_01,    // block -B = C * ones ; need not to be assembled
                        local_matrix_02,    // block C, +/- one on diagonal
                        local_matrix_22,    // final local matrix
                        inv_local_matrix_00; // Ct A- C - alpha_i alpha_j / alpha_tot (schur complement add)

    // auxiliary local vectors
    Vector<double>  alphas, lumping_weights;
    Vector<double> local_velocity;
    Vector<double> local_phead, lambda, local_old_phead,a_p_diff,cond_diff,velocit, steady_flux;

    // vectors provided outside
    Vector<double> local_func,     // local RHS for global system
                   output_vector;   // local par of output global vector

    // auxiliary local values
    double alphas_total,
           element_volume;

    std::map<unsigned int, double> RT_boundary_values; // Dirichlet boundary values

    std::vector<unsigned int> local_dof_indices;
    std::vector< std::vector<unsigned int>  > out_idx;

    //std::vector<unsigned int> face_dofs;

    //const unsigned int z_component;
    //const double density, gravity;
    double dt, time, out_cond, lambda_el;

    std::vector<Tensor < 2, dim> > k_inverse_values;

    LAPACKFullMatrix<double> post_interpolation;

};

template <int dim>
LocalAssembly<dim>::LocalAssembly(unsigned int order)
:
quadrature_formula(order+2),
error_quadrature(order+post_order+1),
face_quadrature_formula(order + 2),

velocity_fe(order,mapping_piola),
trace_fe(order),
pressure_fe(order),
post_pressure_fe(order+post_order),
output_fe (velocity_fe,2,   // flux RT0
           pressure_fe,4,    // pressure, saturation, flux_error, pressure error, Q0
           post_pressure_fe,5),    // postprocessed pressure + difference with analytical

//pressure_fe_values(pressure_fe, quadrature_formula,
//        update_values ),

velocity_fe_values(velocity_fe, quadrature_formula,
        update_values | update_gradients | update_quadrature_points | update_JxW_values),

error_fe_values(post_pressure_fe, error_quadrature,
        update_values | update_gradients | update_quadrature_points | update_JxW_values),

trace_fe_face_values(trace_fe, face_quadrature_formula,
        update_values | update_normal_vectors | update_quadrature_points | update_JxW_values),

velocity_fe_face_values(velocity_fe, face_quadrature_formula,
        update_values),

vec_extr (0),
scal_extr (0),


local_matrix_00(velocity_fe.dofs_per_cell, velocity_fe.dofs_per_cell),
inv_local_matrix_00(velocity_fe.dofs_per_cell, velocity_fe.dofs_per_cell),
local_matrix_01(velocity_fe.dofs_per_cell, pressure_fe.dofs_per_cell),
local_matrix_02(velocity_fe.dofs_per_cell, trace_fe.dofs_per_cell),
local_matrix_22(trace_fe.dofs_per_cell, trace_fe.dofs_per_cell),
local_dof_indices(trace_fe.dofs_per_cell),

local_func(trace_fe.dofs_per_cell),
alphas(velocity_fe.dofs_per_cell),
lumping_weights(velocity_fe.dofs_per_cell),
local_velocity(velocity_fe.dofs_per_cell),
steady_flux(velocity_fe.dofs_per_cell),
local_phead(trace_fe.dofs_per_cell),
local_old_phead(trace_fe.dofs_per_cell),
lambda(trace_fe.dofs_per_cell),
a_p_diff(trace_fe.dofs_per_cell),
cond_diff(trace_fe.dofs_per_cell),
velocit(trace_fe.dofs_per_cell),
output_vector(output_fe.dofs_per_cell),
out_idx(output_fe.n_blocks()),


//density(1.0),
//gravity(1.0),
//z_component(dim-1),

k_inverse_values(quadrature_formula.size())

{
    /// make table that gives output local dof index for given block and index in block
    std::pair<unsigned int,unsigned int> block_index;
    // first allocate block subvectors
    for(unsigned int global_i=0;global_i<output_fe.dofs_per_cell; global_i++) {
        block_index=output_fe.system_to_block_index(global_i);
        out_idx[block_index.first].push_back(0);
    }
    for(unsigned int global_i=0;global_i<output_fe.dofs_per_cell; global_i++) {
        block_index=output_fe.system_to_block_index(global_i);
        out_idx[block_index.first][block_index.second]= global_i;
        cout <<"(block, b.size index):" << block_index.first << " "<<out_idx[block_index.first].size()<<" " << block_index.second<<endl;
    }

    make_interpolation(post_pressure_fe);

    cout << "n_comp: " << output_fe.n_components() << endl;
}

template <int dim>
FullMatrix<double> &LocalAssembly<dim>::get_matrix(bool symmetry)
{
    local_matrix_22=0;
    double cond_new=0, cond_old=0;
    Vector<double> grad_p_new(local_phead);
    Vector<double> grad_p_old(local_old_phead);
    Vector<double> lambda_diff(local_old_phead);
    grad_p_new=0;
    grad_p_old=0;

    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points

    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {
        unsigned int idx = local_dof_indices[i];

        lambda(i) = solution->lambda(idx);
        lambda_diff(i) = solution->lambda_diff(idx);
        local_phead(i) = solution->phead(idx);
        local_old_phead(i) = solution->old_phead(idx);
    }

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        unsigned int idx = local_dof_indices[i];
        cond_new+=lumping_weights(i)*solution->conduct_new(idx);
        cond_old+=lumping_weights(i)*solution->conduct_old(idx);

        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            grad_p_new(i) += inv_local_matrix_00(i,j) * local_phead(j);
            grad_p_old(i) += inv_local_matrix_00(i,j) * local_old_phead(j);
        }

        if (!symmetry) {
            if (solution->cond_type == solution->trapezoid) {
                cond_diff(i)= lambda_el * lumping_weights(i) * solution->conduct_new_diff(idx);
                velocit(i) = grad_p_new(i);
            } else if (solution->cond_type == solution->mid_point_approx) {
                cond_diff(i)= lumping_weights(i) * solution->conduct_new_diff(idx) * solution->lambda(idx)
                        + lambda_diff(i) * (solution->conduct_new(local_dof_indices[i]) - solution->conduct_old(local_dof_indices[i]));
                velocit(i) = lambda_el * grad_p_new(i) + (1-lambda_el) * grad_p_old(i);

            } else if (solution->cond_type == solution->mid_point) {
                cond_diff(i)= lumping_weights(i) * solution->conduct_new_diff(idx) *
                        (solution->lambda_diff(idx)*(local_phead(i) - local_old_phead(i)) + lambda(i));
                velocit(i) = lambda_el * grad_p_new(i) + (1-lambda_el) * grad_p_old(i);

            }

        } else cond_diff(i) = 0;
    }

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {

        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            local_matrix_22(i,j) = cond_new * lambda_el * inv_local_matrix_00(i, j)
                                   + lumping_weights(j)*solution->lambda_diff(local_dof_indices[j]) * (cond_new * grad_p_new(i) - cond_old * grad_p_old(i))
                                   +  cond_diff(j) * velocit(i);
        }

        unsigned int idx = local_dof_indices[i];
        local_matrix_22(i, i) +=
                solution->capacity(idx) *lumping_weights(i) * element_volume / dt;
    }
    //inv_local_matrix_00.print_formatted(cout);
    //alphas.print(cout);
    //inv_local_matrix_00.print_formatted(cout);
    //local_matrix_22.print_formatted(cout);

    //local_matrix_22.add(1.0, inv_local_matrix_00);


    //cout << "local matrix" << endl;
    //local_matrix_22.print_formatted(cout);

    //cout << "local rhs" << endl;
    //local_rhs_2.print(cout);
    //Assert(false, ExcNotImplemented());
    return local_matrix_22;
}

template <int dim>
Vector<double> &LocalAssembly<dim>::get_func() {


    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points
    compute_steady_flux();

    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {
        unsigned int idx = local_dof_indices[i];

        local_func(i) =
                 (solution->saturation(idx) - solution->old_saturation(idx))
                  * lumping_weights(i) * element_volume / dt + steady_flux(i);
    }

    return local_func;
}

template <int dim>
void LocalAssembly<dim>::update_lambda() {

    double cond=0;
    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {

        local_phead(i) =0;
        local_old_phead(i)=0;
        for(unsigned int j=0; j < inv_local_matrix_00.size(0); j++) {
            unsigned int idx = local_dof_indices[j];
            local_phead(i) += inv_local_matrix_00(i,j) * solution->phead(idx);
        }

        // new mid point rule
        cond+=lumping_weights(i)*solution->conductivity_new_midpoint(local_dof_indices[i]);
    }

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        // new midpoint rule
        solution->el_sum(local_dof_indices[i]) += lumping_weights(i) * element_volume;
        solution->lambda_new(local_dof_indices[i]) +=cond * local_phead(i);
        //cout << local_dof_indices[i] << " " << solution->lambda_new(local_dof_indices[i]) << " " << cond * local_phead(i) << endl;
    }
}

template <int dim>
void LocalAssembly<dim>::compute_steady_flux() {

    double cond_new=0, cond_old=0, cond=0;
    for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i) {

        local_phead(i) =0;
        local_old_phead(i)=0;
        for(unsigned int j=0; j < inv_local_matrix_00.size(0); j++) {
            unsigned int idx = local_dof_indices[j];
            local_phead(i) += inv_local_matrix_00(i,j) * solution->phead(idx);
            local_old_phead(i) += inv_local_matrix_00(i,j) * solution->old_phead(idx);
        }
        cond_new+=lumping_weights(i)*solution->conduct_new(local_dof_indices[i]);
        cond_old+=lumping_weights(i)*solution->conduct_old(local_dof_indices[i]);
   }

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        steady_flux(i)= cond_new * lambda_el *local_phead(i) + cond_old * (1-lambda_el) * local_old_phead(i);
    }
}

template <int dim>
Vector<double> &LocalAssembly<dim>::get_s_diff() {

    double conductivity = 0;
    double cond_diff = 0;
    Vector<double> &lambda_s_diff = velocit;
    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points

    local_func = 0;
    /*
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no; // assume that local cell indices are counted as local faces
        unsigned int idx = local_dof_indices[i];

        lambda(i) = solution->lambda(idx);
        lambda_s_diff(i) = solution->lambda_s_diff(idx);
        local_phead(i) = solution->phead(idx);
        local_old_phead(i) = solution->old_phead(idx);

        conductivity += lumping_weights(i) * (solution->conductivity_mid(idx));
        cond_diff += lumping_weights(i) * solution->conductivity_mid_diff(idx) * lambda_s_diff(i) * (local_phead(i) - local_old_phead(i));
        local_func(i) = 0;
    }

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            local_func(i) += cond_diff *  inv_local_matrix_00(i, j) *
                           ( lambda(i)*local_phead (j) + (1-lambda(i))*local_old_phead(j) )
                           +
                           conductivity *  inv_local_matrix_00(i, j)  *
                           ( lambda_s_diff(i)* ( local_phead (j) - local_old_phead(j) ) );
        }
    }*/
    return local_func;
}


/**
 * REinit Fe_values, compute local matrices, vectors, values
 */
template <int dim>
void LocalAssembly<dim>::reinit(typename DoFHandler<dim>::active_cell_iterator cell)
{

    dh_cell = cell;
    typename Triangulation<dim,dim>::active_cell_iterator tria_cell(cell);
    //pressure_fe_values.reinit(tria_cell);
    velocity_fe_values.reinit(tria_cell);
    error_fe_values.reinit(tria_cell);  // TODO: move into error calculation

    // reinit local trace pressure indices
    std::vector<unsigned int> local_face_indices(1);
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        cell->face(face_no)->get_dof_indices(local_face_indices);
        local_dof_indices[face_no]=local_face_indices[0]; // assume that local cell indices are counted as local faces
    }

    element_volume = dh_cell->measure();

    make_A();
    make_C();

    // lm_00 inverse and compute lumping weights
    inv_local_matrix_00.invert(local_matrix_00);
    // scale 00 matrix by 02 matrix
    alphas_total = 0;
    for (unsigned int i = 0; i < local_matrix_00.size(0); ++i) {
        alphas(i) = 0;
        for (unsigned int j = 0; j < local_matrix_00.size(0); ++j) {
            inv_local_matrix_00(i, j) *= local_matrix_02(i, i) * local_matrix_02(j, j);
            alphas(i) += inv_local_matrix_00(i, j);
        }
        alphas_total += alphas(i);
    }
    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            inv_local_matrix_00(i, j) += -alphas(i) * alphas(j) / alphas_total;
        }
    }



    lumping_weights = alphas;
    lumping_weights.scale(1 / alphas_total);

    lambda_el=0;
    for (unsigned int i = 0; i < local_matrix_00.size(0); ++i)
        lambda_el+= lumping_weights(i)*solution->lambda(local_dof_indices[i]);
    //cout << "lmb: " << lambda_el << endl;

    //inv_local_matrix_00.print_formatted(cout);
    //    cout << lumping_weights <<endl;
}

template <int dim>
void LocalAssembly<dim>::make_A() {
    local_matrix_00=0;
    //local_matrix_01=0;
    solution->richards_data->k_inverse.value_list(velocity_fe_values.get_quadrature_points(), k_inverse_values);

    // assembly velocity and cell pressure blocks
    // cycle over quadrature points
    for (unsigned int q = 0; q < quadrature_formula.size(); ++q) {
        for (unsigned int i = 0; i < velocity_fe.dofs_per_cell; ++i) {
            const Tensor < 1, dim > phi_i_u = velocity_fe_values[vec_extr].value(i, q);
            const double div_phi_i_u = velocity_fe_values[vec_extr].divergence(i, q);

            // matrix_00, velocity * velocity
            for (unsigned int j = 0; j < velocity_fe.dofs_per_cell; ++j) {
                const Tensor < 1, dim> phi_j_u = velocity_fe_values[vec_extr].value(j, q);


                local_matrix_00(i,j) += (phi_i_u * k_inverse_values[q] *  phi_j_u) * velocity_fe_values.JxW(q);
            }
/*
            // matrix 01, velocity * element pressure
            for (unsigned int j = 0; j < pressure_fe.dofs_per_cell; ++j) {
                const double phi_j_p = pressure_fe_values[scal_extr].value(j, q);

                local_matrix_01(i,j) += - div_phi_i_u * phi_j_p * velocity_fe_values.JxW(q);
            }
 */
        }
    }
}

/**
 * Actually this do not perform any useful calculation for RT0 elements since C block is just identity matrix, up to the
 * face orientation (can be obtained   in cheaper way)
 * However for other discretization of velocity field as ABF, this is block is notrivial. Then we should also modify
 * calculation of local Schur complements.
 */
template <int dim>
void LocalAssembly<dim>::make_C() {
    local_matrix_02=0;
    // assembly 02 matrix
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        trace_fe_face_values.reinit(dh_cell, face_no);
        velocity_fe_face_values.reinit(dh_cell, face_no);

        // 02 matrix should be just identity matrix up to face orientation
        for (unsigned int q = 0; q < face_quadrature_formula.size(); ++q) {
            for (unsigned int i = 0; i < velocity_fe.dofs_per_cell; ++i) {
                const Tensor<1, dim> phi_i_u = velocity_fe_face_values[vec_extr].value(i, q);

                // matrix 02, velocity * trace pressure
                for (unsigned int j = 0; j < trace_fe.dofs_per_cell; ++j) {
                    const double phi_j_p = trace_fe_face_values.shape_value(j, q);

                    local_matrix_02(i, j) += trace_fe_face_values.normal_vector(q) * phi_i_u * phi_j_p
                            * trace_fe_face_values.JxW(q);
                    Assert( i==j || local_matrix_02(i,j) == 0, ExcDimensionMismatch(i,j));
                }
            }
        }
    }
}


template <int dim>
void LocalAssembly<dim>::apply_bc(typename DoFHandler<dim>::active_cell_iterator cell, std::map<unsigned int, double> &boundary_values)
{
    vector<unsigned int> local_face_indices(1);

    // apply BC conditions
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        trace_fe_face_values.reinit(cell, face_no);
        velocity_fe_face_values.reinit(cell,face_no);

        // apply boundary conditions
        unsigned int b_ind = cell->face(face_no)->boundary_indicator();
        if (b_ind == 255) continue;

        unsigned int n_dofs = cell->get_fe().n_dofs_per_face();
        Assert (  n_dofs== 1, ExcDimensionMismatch (n_dofs, 1) );

        cell->face(face_no)->get_dof_indices(local_face_indices);
        unsigned int idx = local_face_indices[0];


        switch ( solution->richards_data->bc_type(b_ind) ) {
             case RichardsData<dim>::Dirichlet:
             // use boundary values map to eliminate Dirichlet boundary trace pressures
             boundary_values[idx] =
                     solution->richards_data->bc_func(b_ind)->value(
                                 cell->face(face_no)->barycenter()
                             );

             /*  cout << "BC(z,i,v):" << cell->face(face_no)->barycenter()[dim-1] << " "
                  << idx << " "
                  << solution->richards_data->bc_func(b_ind)->value(
                          dh_cell->face(face_no)->barycenter()
                      )
                      <<endl;
            */
             break;

             case RichardsData<dim>::Neuman:
             // add BC to RHS - not necessary for homogeneous Neuman
             /*
             for (unsigned int j = 0; j < trace_fe.dofs_per_cell; ++j) {
                 const double phi_j_p = trace_fe_face_values[scal_extr].value(j, q);

                 local_rhs_02(j) += phi_j_p * neuman_bc* trace_fe_face_values.JxW(q);

             } */
             break;

             default:
                     Assert (false, ExcNotImplemented());
             break;

        }
    }
}


/**
 * After reinit by solution dof handler. This should reconstruct local velocity and element pressure from solution.
 */
template <int dim>
void LocalAssembly<dim>::output_evaluate(double & bc_flux_total, double & volume_total)
{

    Vector<double> trace_head(trace_fe.dofs_per_cell);
    Vector<double> weight_trace_phead(trace_fe.dofs_per_cell);
    Vector<double> trace_res(trace_fe.dofs_per_cell);
    Vector<double> aux(trace_fe.dofs_per_cell);
    aux=0;
    lambda=0;


    output_vector=0;
    double conductivity =0;
    double el_phead=0;
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no;
        unsigned int idx = local_dof_indices[i];

        local_phead(i) = solution->phead(idx);
        trace_res(i) = solution->residual(idx);
        trace_head(i)= solution->pressure(idx);
        local_old_phead(i) = solution->old_pressure(idx);

        double cap = max(solution->capacity(idx), 1e-3);

        //aux(i) = solution->lambda_new(idx)/solution->el_sum(idx);
                //0.75* solution->pressure(idx) + 0.25*solution->old_pressure(idx) + dt*0.25*solution->lambda_new(idx)/solution->el_sum(idx)/cap;
        //cout << cap << " " << aux(i) << " " << solution->lambda_new(idx) <<endl;

        lambda(i) = solution->lambda(idx); //solution->residual(idx) - (solution->saturation(idx) - solution->old_saturation(idx))/dt * solution->el_sum(idx);
                          //element_volume*lumping_weights(i)*solution->new_sat_diff(idx)/dt;
        el_phead    += trace_head(face_no) * lumping_weights(face_no);

         // to get velocity compatible with saturation we should use  lambda weighted phead
//        conductivity += lumping_weights(i)*solution->conductivity_mid(idx);
        //solution->richards_data->fk(trace_head(i))* lumping_weights(i);

        //cout << "i f l:"<<local_dof_indices[0]<<" "<<face_no<<" "<<local_lambda(face_no)<<
        //        " "<<trace_phead(face_no)<<" "<<weight_trace_phead(face_no)<<endl;
    }


    //cout << "head: " << dh_cell->barycenter()[dim-1] << " "
    //     << output_el_head << " " << output_el_phead << endl;
    //if (solution->cond_type == solution->trapezoid) lambda_el  = 1.0;
    compute_steady_flux();
    Vector<double> half_velocity(local_velocity);
    half_velocity = steady_flux;
    lambda_el  = 1.0;
    compute_steady_flux();

    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no; // assume that local cell indices are counted as local faces

        double sat_old = solution->old_saturation(local_dof_indices[i]);
        double sat_new = solution->saturation(local_dof_indices[i]);

        local_velocity(i) = -( element_volume*lumping_weights(i)*
                    (sat_new - sat_old) / dt + steady_flux(i)) * local_matrix_02(i,i) ; //
        half_velocity(i) = -( element_volume*lumping_weights(i)*
                    (sat_new - sat_old) / dt + half_velocity(i)) * local_matrix_02(i,i) ; //
        //local_velocity(i) = 0;
        //for(unsigned int j=0; j<local_velocity.size(); j++) {
        //    local_velocity(i) = ( lambda(i)*local_phead (j) + (1-lambda(i))*local_old_phead(j) ) * inv_local_matrix_00(i,j) * conductivity;
        // }
        //local_velocity(i)*= -local_matrix_02(i,i); //
        output_vector( out_idx[velocity_bl][i] ) = local_velocity(i);
        output_vector( out_idx[half_flux_bl][i] ) = half_velocity(i);


        output_vector( out_idx[saturation_bl][0] )
          += lumping_weights(i) * sat_new;

        volume_total += element_volume*lumping_weights(i) * sat_new;
        if (dh_cell->face(face_no)->boundary_indicator() != 255) bc_flux_total += half_velocity(i)*local_matrix_02(i,i)*dt;
    }

    output_vector( out_idx[pressure_bl][0] ) = el_phead;

    out_interpolate(trace_head, post_p_bl);
    //trace_phead.add( -1.0, old_trace_phead);
    out_interpolate(trace_res, post_residual);

    //get_func();
    out_interpolate(aux, post_aux);

    out_interpolate(lambda, post_lambda);
    out_interpolate(local_old_phead, post_old_p_bl);

    //cout << "cond diff: " << dh_cell->barycenter()[dim-1] << " "

    //double anal_sol = solution->richards_data->anal_sol->value(dh_cell->barycenter());
    //double p_error = el_phead - anal_sol;
    if (solution->richards_data->has_exact_solution()) {
        double p_error = compute_p_error();

        output_vector(out_idx[p_error_bl][0]) = p_error; //dh_cell->measure() * p_error * p_error;

        //double flux = (local_velocity(2) + local_velocity(3)) / 2.0;
        //double anal_flux = solution->richards_data->anal_flux->value(dh_cell->barycenter());
        //double q_error = flux - anal_flux;
        double q_error = compute_q_error();
        output_vector(out_idx[q_error_bl][0]) = q_error; //dh_cell->measure() * q_error * q_error;
    }

}

template <int dim>
void LocalAssembly<dim>::set_lambda() {

}


// L2 norm squared of pressure error over one element
template <int dim>
double LocalAssembly<dim>::compute_p_error() {

    const std::vector< Point<dim> > &q_points = error_fe_values.get_quadrature_points();

    double error =0;
    for(unsigned int q=0; q < error_quadrature.size(); ++q) {
        double post_pressure = 0;
        for(unsigned int i=0; i< out_idx[post_p_bl].size(); ++i) {
            post_pressure += output_vector( out_idx[post_p_bl][i] ) * error_fe_values[scal_extr].value(i,q);
        }

        post_pressure = solution->richards_data->pressure( post_pressure, q_points[q] );
        double anal_sol = solution->richards_data->anal_sol->value( q_points[q] );
        double diff = post_pressure - anal_sol;
        error += diff * diff * error_fe_values.JxW(q);
    }

    return error;
}

// L2 norm squared of flux error over one element
template <int dim>
double LocalAssembly<dim>::compute_q_error() {

    const std::vector< Point<dim> > &q_points = error_fe_values.get_quadrature_points();

  //  double flux = (local_velocity(2) + local_velocity(3) ) / 2.0;
  //  double anal_flux =  solution->richards_data->anal_flux->value(dh_cell->barycenter());

    Tensor<1,dim> velocity;
    double error =0;
    for(unsigned int q=0; q < error_quadrature.size(); ++q) {
        velocity = 0;
        for(unsigned int i=0; i< out_idx[velocity_bl].size(); ++i) {
            velocity += output_vector( out_idx[velocity_bl][i] ) * velocity_fe_values[vec_extr].value(i,q);
        }
        double anal_flux = solution->richards_data->anal_flux->value( q_points[q] );
        double diff = velocity[1] - anal_flux;
        error += (diff * diff + velocity[0]*velocity[0]) * velocity_fe_values.JxW(q);
    }

    return error;
}


template <int dim>
void LocalAssembly<dim>::make_interpolation(FE_DGQ<dim> &fe) {

    post_interpolation.reinit(GeometryInfo<dim>::faces_per_cell, post_pressure_fe.n_dofs_per_cell());
    /*
    // make interpolation
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
       unsigned i = face_no;    // row , trace point index
       Point<dim> center(true); // initialize by zeros
       for(unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face; ++j) {
           Point<dim> vtx( GeometryInfo<dim>::unit_cell_vertex( GeometryInfo< dim >::face_to_cell_vertices(i,j) ) );
           center += vtx;

       }
       center /= GeometryInfo<dim>::vertices_per_face;
       cout << center << endl;
       cout << "oreeeer: " << post_pressure_fe.get_degree() <<endl;
       for (unsigned int j = 0; j < post_pressure_fe.n_dofs_per_cell(); ++j)  // column DGQ shape function index
            post_interpolation(i,j) = post_pressure_fe.shape_value(j, center);
    }
    post_interpolation.print_formatted(cout);
    post_interpolation.compute_inverse_svd();
    Vector<double> one(4), res(4);
    one=0; one(0)=1;

    post_interpolation.vmult(res,one);
    cout << res <<endl;*/
    // columns - trace points
    // rows - DGQ shape functions
    post_interpolation(0,0) =3 * 0.25;
    post_interpolation(0,1) =-1 * 0.25;
    post_interpolation(0,2) =3 * 0.25;
    post_interpolation(0,3) =-1 * 0.25;
    post_interpolation(1,0) =-1 * 0.25;
    post_interpolation(1,1) =3 * 0.25;
    post_interpolation(1,2) =3 * 0.25;
    post_interpolation(1,3) =-1 * 0.25;
    post_interpolation(2,0) =3 * 0.25;
    post_interpolation(2,1) =-1 * 0.25;
    post_interpolation(2,2) =-1 * 0.25;
    post_interpolation(2,3) =3 * 0.25;
    post_interpolation(3,0) =-1 * 0.25;
    post_interpolation(3,1) =3 * 0.25;
    post_interpolation(3,2) =-1 * 0.25;
    post_interpolation(3,3) =3 * 0.25;

}

/**
 * Interpolation of pressure into post-processed FE space. However, it can be used also for interpolation of different fields.
 */
template <int dim>
void LocalAssembly<dim>::out_interpolate(Vector<double> &trace_phead, unsigned int out_block) {
    // compute interpolation through matrix
        for(unsigned int i = 0; i< post_pressure_fe.n_dofs_per_cell(); ++i) {
            output_vector( out_idx[out_block][i] )=0;
            for(unsigned int j = 0; j < local_dof_indices.size(); ++j) {
                output_vector(out_idx[out_block][i]) += post_interpolation(i,j) * trace_phead(j);
            }
        }
}

/**
 * Interpolation of trace values into DGP, we use Monomial version, since
 * "DGP" use normalized Legendre polynomials 1, sqrt(3)(2x-1), sqrt(5)(6x^2 - 6x +1)
 * and interpolation matrix is not easily computable. Moreover since we perform only local interpolation
 * we have no problems with preconditioning.
 *
 * TODO: Use two more degrees of freedom to take infromation form velocity field.
 * Without this information we get saddle point interpolation on one cell even in 1D solved in 2D.
 */
template <int dim>
void LocalAssembly<dim>::make_interpolation(FE_DGPMonomial<dim> &fe) {

    post_interpolation.reinit(post_pressure_fe.n_dofs_per_cell(), GeometryInfo<dim>::faces_per_cell);
    Assert( post_pressure_fe.n_dofs_per_cell() == 6, ExcDimensionMismatch(post_pressure_fe.n_dofs_per_cell() ,6));
    Assert( GeometryInfo<dim>::faces_per_cell == 4,  ExcDimensionMismatch(GeometryInfo<dim>::faces_per_cell ,4));

    /*
    // make interpolation
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
       unsigned i = face_no;    // row , trace point index
       Point<dim> center(true); // initialize by zeros
       for(unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face; ++j) {
           Point<dim> vtx( GeometryInfo<dim>::unit_cell_vertex( GeometryInfo< dim >::face_to_cell_vertices(i,j) ) );
           center += vtx;

       }
       center /= GeometryInfo<dim>::vertices_per_face;
       cout << center << endl;
       cout << "oreeeer: " << post_pressure_fe.get_degree() <<endl;
       for (unsigned int j = 0; j < post_pressure_fe.n_dofs_per_cell(); ++j)  // column DGQ shape function index
            post_interpolation(i,j) = post_pressure_fe.shape_value(j, center);
    }
    post_interpolation.print_formatted(cout);
    post_interpolation.compute_inverse_svd();
    Vector<double> one(4), res(4);
    one=0; one(0)=1;

    post_interpolation.vmult(res,one);
    cout << res <<endl;*/
    // columns - trace points
    // rows - DGQ shape functions
    post_interpolation(0,0) =3.0 / 4.0;
    post_interpolation(0,1) =-1.0/4;
    post_interpolation(0,2) =3.0/4;
    post_interpolation(0,3) =-1.0/4;
    post_interpolation(1,0) =-2;
    post_interpolation(1,1) =0;
    post_interpolation(1,2) =1;
    post_interpolation(1,3) =1;
    post_interpolation(2,0) =1;
    post_interpolation(2,1) =1;
    post_interpolation(2,2) =-2;
    post_interpolation(2,3) =0;
    post_interpolation(3,0) =0;
    post_interpolation(3,1) =0;
    post_interpolation(3,2) =0;
    post_interpolation(3,3) =0;
    post_interpolation(4,0) =1;
    post_interpolation(4,1) =1;
    post_interpolation(4,2) =-1;
    post_interpolation(4,3) =-1;
    post_interpolation(5,0) =-1;
    post_interpolation(5,1) =-1;
    post_interpolation(5,2) =1;
    post_interpolation(5,3) =1;

    // check
    std::vector< Point<dim> > face_bary(4);
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
       unsigned i = face_no;    // row , trace point index
       Point<dim> center(true); // initialize by zeros
       for(unsigned int j = 0; j < GeometryInfo<dim>::vertices_per_face; ++j) {
           Point<dim> vtx( GeometryInfo<dim>::unit_cell_vertex( GeometryInfo< dim >::face_to_cell_vertices(i,j) ) );
           center += vtx;

       }
       face_bary[i] = center / GeometryInfo<dim>::vertices_per_face;
    }

    for(int i=0;i<4; ++i) {
        for(int k=0; k< 4; ++k) {
        double value=0;
        //cout << face_bary[k] << endl;
        for (unsigned int j = 0; j < post_pressure_fe.n_dofs_per_cell(); ++j) {  // column DGQ shape function index
            value+=post_interpolation(j,i)* post_pressure_fe.shape_value(j, face_bary[k]);
            cout << post_pressure_fe.shape_value(j, face_bary[k]) << " ";
        }

        cout <<":"<< value<< " "<<endl;
    }
    cout << endl;
    }

}

/* Returns local variation of flux and pressure. I.e. sum of absolute value of differences in every direction of reference element. */
// computes  approximation for second order term of Taylor expansion of h(t + dt)
template <int dim>
void LocalAssembly<dim>::compute_add_variation(double &flux_var, double &head_var) {
    head_var += fabs(solution->pressure(local_dof_indices[0]) - solution->pressure(local_dof_indices[1]))
            + fabs(solution->pressure(local_dof_indices[2]) - solution->pressure(local_dof_indices[3]));

    Vector<double> face_area(4);

    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        trace_fe_face_values.reinit(dh_cell, face_no);
        face_area(face_no) = 0;
         // compute face area
         for (unsigned int q = 0; q < face_quadrature_formula.size(); ++q) face_area(face_no)+= trace_fe_face_values.JxW(q);
     }


    flux_var += fabs(local_velocity(0)/face_area(0) - local_velocity(1)/face_area(1)) +
            fabs(local_velocity(2)/face_area(2) - local_velocity(3)/face_area(3));
}

/*
template <int dim>
void LocalAssembly<dim>::out_interpolate_pressure(FE_DGPMonomial<dim> &fe, Vector<double> &trace_phead) {


    // compute interpolation through matrix
    if (fe.get_degree() == 2) {
        for(unsigned int i = 0; i< post_pressure_fe.n_dofs_per_cell(); ++i) {
            output_vector( out_idx[post_p_bl][i] )=0;
            for(unsigned int j = 0; j < local_dof_indices.size(); ++j) {
                output_vector(out_idx[post_p_bl][i]) += post_interpolation(i,j) * trace_phead(j);
            }
            //Point<dim> vtx( dh_cell->vertex(i) );
            output_vector(out_idx[post_p_diff_bl][i])= 0;
            //output_vector(out_idx[post_p_bl][i])
            //        - solution->richards_data->anal_sol->value( vtx );
        }

    }

    // compute error
}
*/


#endif /* LOCAL_ASSEMBLY_HH_ */
