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
        lambda_s_diff.reinit(dh.n_dofs());

        conductivity_mid.reinit(dh.n_dofs());
        conductivity_mid_diff.reinit(dh.n_dofs());
        //old_conductivity.reinit(dh.n_dofs());
        saturation.reinit(dh.n_dofs());
        capacity.reinit(dh.n_dofs());
        old_saturation.reinit(dh.n_dofs());
        residual.reinit(dh.n_dofs());

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
        dt_=dt;

        for(unsigned int i=0;i<phead.size();i++) {
            old_phead(i) = phead(i);
            old_saturation(i) = saturation(i);
        }
        last_lmb_change=false;
    }

    void update(double s_param) {
        double p, p_scaled;

        int  n_lmb_change = 0;
        for(unsigned int i=0;i<phead.size();i++) {
            p = pressure(i);
            saturation(i)=richards_data->fq(p,capacity(i));

            // lambda
            /*
            if (p>0) {
                lambda(i)=1;
                lambda_diff(i)=0;
                lambda_s_diff(i)=0;
            } else {
                p_scaled = p / (richards_data->cap_arg_max_/3.0);
                if (p_scaled > 1) {
                    lambda(i) = 1 - s_param * 0.5;
                    lambda_diff(i)=0;
                    lambda_s_diff(i)=-0.5;
                } else {
                    lambda(i) = 1 - s_param * 0.5 * p_scaled;
                    lambda_diff(i)= - s_param * 0.5/ (richards_data->cap_arg_max_/3.0);
                    lambda_s_diff(i)= - 0.5 * p_scaled;
                }
            }*/
            lambda(i)=1;
            lambda_diff(i)=0;
            lambda_s_diff(i)=0;


            p=lambda(i)*p + (1-lambda(i))*old_pressure(i);
            conductivity_mid(i)=richards_data->fk(p, conductivity_mid_diff(i));
        }
    }

    inline double pressure(unsigned int i)
        { return richards_data->pressure(phead(i), z_coord(i)); }

    inline double old_pressure(unsigned int i)
        { return richards_data->pressure(old_phead(i), z_coord(i)); }

    //inline double phead_mid(unsigned int i)
    //        { return lambda(i)*phead(i)+(1-lambda(i))*old_phead(i); }




    Vector<double> z_coord;
    Vector<double> phead;
    Vector<double> old_phead;

    Vector<double> lambda;
    Vector<double> lambda_diff;
    Vector<double> lambda_s_diff;

    Vector<double> conductivity_mid;
    Vector<double> conductivity_mid_diff;
    //Vector<double> old_conductivity;
    Vector<double> saturation;
    Vector<double> capacity;
    Vector<double> old_saturation;
    Vector<double> residual;

    RichardsData<dim> *richards_data;
    double dt_;
    bool last_lmb_change;

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
        pressure_bl=1,
        saturation_bl=2,
        p_error_bl=3,
        q_error_bl=4,
        post_p_bl=5,
        post_p_diff_bl=6
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
    void output_evaluate();
    void set_lambda();
    double get_p_error() { return  output_vector( out_idx[p_error_bl][0] ); }
    double get_q_error() { return  output_vector( out_idx[q_error_bl][0] ); }

    void make_A();
    void make_C();
    void apply_bc(typename DoFHandler<dim>::active_cell_iterator cell, std::map<unsigned int, double> &boundary_values);
    FullMatrix<double> &get_matrix(bool symmetric);
    Vector<double> &get_func();
    Vector<double> &get_s_diff();


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
    Vector<double> local_phead, lambda, local_old_phead,a_p_diff,cond_diff,velocit;

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
    double dt, time;

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
output_fe (velocity_fe,1,   // flux RT0
           pressure_fe,4,    // pressure, saturation, flux_error, pressure error, Q0
           post_pressure_fe,2),    // postprocessed pressure + difference with analytical

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
}

template <int dim>
FullMatrix<double> &LocalAssembly<dim>::get_matrix(bool symmetry)
{
    local_matrix_22=0;
    double conductivity = 0;

    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points

    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no; // assume that local cell indices are counted as local faces
        unsigned int idx = local_dof_indices[i];

        lambda(i) = solution->lambda(idx);
        local_phead(i) = solution->phead(idx);
        local_old_phead(i) = solution->old_phead(idx);
        conductivity += lumping_weights(i) * solution->conductivity_mid(idx);

        if (!symmetry) {
            cond_diff(i)= lumping_weights(i) * solution->conductivity_mid_diff(idx) *
                          (solution->lambda_diff(idx)*(local_phead(i) - local_old_phead(i)) + lambda(i));
        } else cond_diff(i) = 0;

        velocit(i)=0;
        a_p_diff(i)=0;
        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            if (! symmetry) velocit(i) += inv_local_matrix_00(i, j)  * ( lambda(i)*local_phead (j) + (1-lambda(i))*local_old_phead(j) );
            a_p_diff(i) += inv_local_matrix_00(i,j) * (local_phead(j) - local_old_phead(j));
        }
        lambda(i) = solution->lambda(idx);

    }

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {

        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            local_matrix_22(i,j) = lambda(i) * conductivity * inv_local_matrix_00(i, j)
                                   +  cond_diff(j) * velocit(i);
        }

        unsigned int idx = local_dof_indices[i];
        local_matrix_22(i, i) += lumping_weights(i) * element_volume / dt * solution->capacity(idx)
                + solution->lambda_diff(idx) * conductivity * a_p_diff(i);
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

    double conductivity = 0;

    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points

    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no; // assume that local cell indices are counted as local faces
        unsigned int idx = local_dof_indices[i];

        lambda(i) = solution->lambda(idx);
        local_phead(i) = solution->phead(idx);
        local_old_phead(i) = solution->old_phead(idx);

        conductivity += lumping_weights(i) * (solution->conductivity_mid(idx));
        local_func(i) = (solution->saturation(idx) - solution->old_saturation(idx)) * lumping_weights(i) * element_volume / dt;
    }
    //cout << conductivity << " func: " << local_func << endl;

    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            local_func(i) += conductivity *  inv_local_matrix_00(i, j)  *
                          ( lambda(i)*local_phead (j) + (1-lambda(i))*local_old_phead(j) );
        }
    }
    //cout << conductivity << " func: " << local_func << endl;

    return local_func;
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
    }
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
                                 dh_cell->face(face_no)->barycenter()
                             );
             /*cout << "BC(z,i,v):" << dh_cell->face(face_no)->barycenter()[dim-1] << " "
                  << idx << " "
                  << solution->richards_data->bc_func(b_ind)->value(
                          dh_cell->face(face_no)->barycenter()
                      )<<endl;*/


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
void LocalAssembly<dim>::output_evaluate()
{
    Vector<double> trace_phead(trace_fe.dofs_per_cell);
    Vector<double> trace_head(trace_fe.dofs_per_cell);
    Vector<double> weight_trace_phead(trace_fe.dofs_per_cell);
    Vector<double> trace_res(trace_fe.dofs_per_cell);
    Vector<double> local_velocity(velocity_fe.dofs_per_cell);

    output_vector=0;
    double conductivity =0;
    double el_phead=0;
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no;
        unsigned int idx = local_dof_indices[i];

        trace_phead(i) = solution->phead(idx);
        trace_res(i) = solution->residual(idx);
        trace_head(i)= solution->pressure(idx);
        el_phead    += trace_head(face_no) * lumping_weights(face_no);

        weight_trace_phead(i) =  trace_phead(i); // to get velocity compatible with saturation we should use  lambda weighted phead
        conductivity += solution->richards_data->fk(trace_head(i))* lumping_weights(i);

        //cout << "i f l:"<<local_dof_indices[0]<<" "<<face_no<<" "<<local_lambda(face_no)<<
        //        " "<<trace_phead(face_no)<<" "<<weight_trace_phead(face_no)<<endl;
    }


    //cout << "head: " << dh_cell->barycenter()[dim-1] << " "
    //     << output_el_head << " " << output_el_phead << endl;

    inv_local_matrix_00.vmult(local_velocity, weight_trace_phead);

    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        unsigned int i=face_no; // assume that local cell indices are counted as local faces

        double sat_old = solution->old_saturation(local_dof_indices[i]);
        double sat_new = solution->saturation(local_dof_indices[i]);


        local_velocity(i) = -( element_volume*lumping_weights(i)* (sat_new - sat_old) / dt + conductivity*local_velocity(i) ) * local_matrix_02(i,i); //
        output_vector( out_idx[velocity_bl][i] ) = local_velocity(i);

        output_vector( out_idx[saturation_bl][0] )
          += lumping_weights(i) * sat_new;
    }

    output_vector( out_idx[pressure_bl][0] ) = el_phead;

    out_interpolate(trace_head, post_p_bl);
    //trace_phead.add( -1.0, old_trace_phead);
    out_interpolate(trace_res, post_p_diff_bl);



    //cout << "cond diff: " << dh_cell->barycenter()[dim-1] << " "

    //double anal_sol = solution->richards_data->anal_sol->value(dh_cell->barycenter());
    //double p_error = el_phead - anal_sol;
    if (solution->richards_data->has_exact_solution()) {
        double p_error = compute_p_error();
        output_vector(out_idx[p_error_bl][0]) = p_error; //dh_cell->measure() * p_error * p_error;

        double flux = (local_velocity(2) + local_velocity(3)) / 2.0;
        double anal_flux = solution->richards_data->anal_flux->value(dh_cell->barycenter());
        double q_error = flux - anal_flux;
        output_vector(out_idx[q_error_bl][0]) = dh_cell->measure() * q_error * q_error;
    }
    output_vector(out_idx[p_error_bl][0]) = conductivity;
}

template <int dim>
void LocalAssembly<dim>::set_lambda() {

}


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
