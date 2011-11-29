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

#include <base/tensor_base.h>
#include <base/quadrature_lib.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_dg_vector.h>
#include <fe/fe_raviart_thomas.h>
#include <base/polynomials_raviart_thomas.h>
#include <fe/fe_face.h>

#include <fe/fe_values.h>

#include <spatial_functions.hh>

using namespace dealii;


template <int dim>
class LocalAssembly {
public:
    LocalAssembly(unsigned int order);
    void set_vectors(const Vector<double> &sol,const  Vector<double> &old_sol)
    {
        solution=&sol;
        old_solution=&old_sol;
    }

    void set_data(struct RichardsData<dim> * data)
    {
        richards_data=data;
    }

    void set_dt(const double in_dt, const double in_t) {
        dt = in_dt;
        time = in_t;
    }

    void reinit(typename DoFHandler<dim>::active_cell_iterator cell);
    void output_evaluate();

    void make_A();
    void make_C();
    void apply_bc(std::map<unsigned int, double> &boundary_values);
    FullMatrix<double> &get_matrix() {return local_matrix_22;}
    Vector<double> &get_rhs() {return local_rhs_2;}
    Vector<double> &get_output_velocity() {return local_velocity;}
    double get_output_el_sat() {return output_el_sat;}
    double get_output_el_phead() {return output_el_phead;}
    double get_output_el_head() {return output_el_head;}

    RichardsData<dim> *richards_data;

private:
    typename DoFHandler<dim>::active_cell_iterator dh_cell;

    QGauss<dim> quadrature_formula;
    QGauss < dim - 1 > face_quadrature_formula;
    const unsigned int n_q_points;
    const unsigned int n_face_q_points;

    FE_DGVector<PolynomialsRaviartThomas<dim>, dim> velocity_fe;
    FE_FaceQ<dim> trace_fe;
    FE_DGQ<dim> pressure_fe;
    FE_DGQ<dim> postprocessed_fe;

    FEValues<dim> pressure_fe_values;
    FEValues<dim> velocity_fe_values;
    FEFaceValues<dim> trace_fe_face_values;
    FEFaceValues<dim> velocity_fe_face_values;
    const FEValuesExtractors::Vector vec_extr;
    const FEValuesExtractors::Scalar scal_extr;


    FullMatrix<double>  local_matrix_00,    // block A
                        local_matrix_01,    // block -B = C * ones ; need not to be assembled
                        local_matrix_02,    // block C, +/- one on diagonal
                        local_matrix_22,    // final local matrix
                        inv_local_matrix_00; // Ct A- C - alpha_i alpha_j / alpha_tot (schur complement add)
    Vector<double> local_rhs_2, alphas, lumping_weights;
    Vector<double> local_velocity;
    double output_el_phead;
    double output_el_head;
    double output_el_sat;
    double alphas_total;

    std::map<unsigned int, double> RT_boundary_values; // Dirichlet boundary values

    const Vector<double> * solution;
    const Vector<double> * old_solution;


    std::vector<unsigned int> local_dof_indices;
    //std::vector<unsigned int> face_dofs;

    const unsigned int z_component;
    const double density, gravity;
    double dt, time;

    std::vector<Tensor < 2, dim> > k_inverse_values;
    double conductivity;
    double element_volume;
};

template <int dim>
LocalAssembly<dim>::LocalAssembly(unsigned int order)
:
quadrature_formula(order+2),
face_quadrature_formula(order + 2),
n_q_points(quadrature_formula.size()),
n_face_q_points(face_quadrature_formula.size()),

velocity_fe(order,mapping_piola),
pressure_fe(order),
trace_fe(order),
postprocessed_fe(order+2),

pressure_fe_values(pressure_fe, quadrature_formula,
        update_values ),

velocity_fe_values(velocity_fe, quadrature_formula,
        update_values | update_gradients | update_quadrature_points | update_JxW_values),


//FEValues<dim> sat_fe_values(sat_fe, quadrature_formula,
//        update_values);
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
local_dof_indices(trace_fe.dofs_per_face),

local_rhs_2(trace_fe.dofs_per_cell),
alphas(velocity_fe.dofs_per_cell),
lumping_weights(velocity_fe.dofs_per_cell),
local_velocity(velocity_fe.dofs_per_cell),


density(1.0),
gravity(1.0),
z_component(dim-1),

k_inverse_values(n_q_points)

{
}

template <int dim>
void LocalAssembly<dim>::reinit(typename DoFHandler<dim>::active_cell_iterator cell)
{

    local_matrix_00=0;
    local_matrix_01=0;
    local_matrix_02=0;
    local_matrix_22=0;
    local_rhs_2=0;

    dh_cell = cell;
    typename Triangulation<dim,dim>::active_cell_iterator tria_cell(cell);
    pressure_fe_values.reinit(tria_cell);
    velocity_fe_values.reinit(tria_cell);

    make_A();
    make_C();


    // lm_00 inverse and compute lumping weights
    inv_local_matrix_00.invert(local_matrix_00);
//    local_matrix_00.print_formatted(cout);
//    inv_local_matrix_00.print_formatted(cout);
    // scale 00 matrix by 02 matrix
    alphas_total = 0;
    //local_matrix_02.print_formatted(cout);
    for (unsigned int i = 0; i < local_matrix_00.size(0); ++i) {
        alphas(i) = 0;
        for (unsigned int j = 0; j < local_matrix_00.size(0); ++j) {
            inv_local_matrix_00(i, j) *= local_matrix_02(i, i) * local_matrix_02(j, j);
            alphas(i) += inv_local_matrix_00(i, j);
        }
        alphas_total += alphas(i);
    }

    lumping_weights = alphas;
    lumping_weights.scale(1 / alphas_total);
    //alphas.print(cout);
    //lumping_weights.print(cout);

    element_volume = cell->measure();
    conductivity = 0;

    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points



    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        cell->face(face_no)->get_dof_indices(local_dof_indices);
        unsigned int i=face_no; // assume that local cell indices are counted as local faces
        Point<dim> barycenter = cell->face(face_no)->barycenter();

        double trace_phead = (*solution)(local_dof_indices[0]);
        double trace_pressure = richards_data->pressure( trace_phead, barycenter);
        double trace_pressure_old = richards_data->pressure( (*old_solution)(local_dof_indices[0]) , barycenter);

        //double sat=f.val();
        double capacity;
        //double capacity = 0.0;
//        double sat_new = capacity* p.val(); //f.val();
//        double sat_old = capacity* trace_pressure_old; //richards_data->fq(trace_pressure_old);
        double sat_new = richards_data->fq(trace_pressure, capacity);
        double sat_old = richards_data->fq(trace_pressure_old);

        //cout << "tph sn so c:" << trace_phead << " "<<sat_new<<" "<<sat_old<<" "<<capacity<<endl;

        local_matrix_22(i, i) = lumping_weights(i) * element_volume / dt * capacity;
        local_rhs_2(i) = lumping_weights(i) * element_volume / dt * (capacity * trace_phead - (sat_new - sat_old));

        conductivity += lumping_weights(i) * richards_data->fk(trace_pressure);
/*
        cout << i << "(lw, ev, dt, capcap: " << capacity << "lm: " << local_matrix_22(i, i) << endl;
        cout << "rhs: " << local_rhs_2(i) << "con: " << conductivity << endl;
        cout << p.val() << endl;
        cout << richards_data->inv_fk(p.val()) << endl;
*/
    }

    //conductivity = 0.1;
    // make schur
    // modify inverse matrix

    //inv_local_matrix_00.print_formatted(cout);
    for (unsigned int i = 0; i < local_matrix_22.size(0); ++i) {
        for (unsigned int j = 0; j < local_matrix_22.size(0); ++j) {
            inv_local_matrix_00(i, j) += -alphas(i) * alphas(j) / alphas_total;
            inv_local_matrix_00(i, j) *= conductivity;
        }
    }
    //inv_local_matrix_00.print_formatted(cout);
    //alphas.print(cout);
    //inv_local_matrix_00.print_formatted(cout);
    //local_matrix_22.print_formatted(cout);

    local_matrix_22.add(1.0, inv_local_matrix_00);


    //cout << "local matrix" << endl;
    //local_matrix_22.print_formatted(cout);

    //cout << "local rhs" << endl;
    //local_rhs_2.print(cout);
    //Assert(false, ExcNotImplemented());

}

template <int dim>
void LocalAssembly<dim>::make_A() {
    richards_data->k_inverse.value_list(velocity_fe_values.get_quadrature_points(), k_inverse_values);

    // assembly velocity and cell pressure blocks
    // cycle over quadrature points
    for (unsigned int q = 0; q < n_q_points; ++q) {
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

template <int dim>
void LocalAssembly<dim>::make_C() {
    // assembly 02 matrix
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        trace_fe_face_values.reinit(dh_cell, face_no);
        velocity_fe_face_values.reinit(dh_cell, face_no);

        // 02 matrix should be just identity matrix up to face orientation
        for (unsigned int q = 0; q < n_face_q_points; ++q) {
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
void LocalAssembly<dim>::apply_bc(std::map<unsigned int, double> &boundary_values)
{
    double s;

    // assembly 02 matrix & apply BC conditions
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        trace_fe_face_values.reinit(dh_cell, face_no);
        velocity_fe_face_values.reinit(dh_cell,face_no);

        // apply boundary conditions
        unsigned int b_ind = dh_cell->face(face_no)->boundary_indicator();
        if (b_ind == 255) continue;

        unsigned int n_dofs = dh_cell->get_fe().n_dofs_per_face();
        Assert (  n_dofs== 1, ExcDimensionMismatch (n_dofs, 1) );

        switch ( richards_data->bc_type(b_ind) ) {
             case RichardsData<dim>::Dirichlet:
             // use boundary values map to eliminate Dirichlet boundary trace pressures
             boundary_values[dh_cell->face(face_no)->dof_index(0)] =
                     richards_data->bc_func(b_ind)->value(
                                 dh_cell->face(face_no)->barycenter()
                             );
             //cout << "BC(z,s,v):" << dh_cell->face(face_no)->barycenter()[dim-1] << " "
             //     << s << " "
             //     << -tan( (exp(s)-1) / (exp(s) +1 )) << " "<<endl;
                     //face_bc->value(trace_fe_face_values.quadrature_point(0));

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
    output_el_phead = 0;
    output_el_sat = 0;


    Vector<double> trace_phead(trace_fe.dofs_per_cell);
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        dh_cell->face(face_no)->get_dof_indices(local_dof_indices);
        trace_phead(face_no) = (*solution)(local_dof_indices[0]);
        output_el_phead += trace_phead(face_no) * lumping_weights(face_no);

    }
    output_el_head = richards_data->pressure(output_el_phead,dh_cell->barycenter());
    //output_el_head = trace_phead(3);
    //cout << "head: " << dh_cell->barycenter()[dim-1] << " "
    //     << output_el_head << " " << output_el_phead << endl;

    inv_local_matrix_00.vmult(local_velocity, trace_phead);


    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        dh_cell->face(face_no)->get_dof_indices(local_dof_indices);
        unsigned int i=face_no; // assume that local cell indices are counted as local faces
        Point<dim> barycenter = dh_cell->face(face_no)->barycenter();

        double sat_old = richards_data->fq( richards_data->pressure( (*old_solution)(local_dof_indices[0]) , barycenter) );
        double sat_new = richards_data->fq( richards_data->pressure( trace_phead(i) , barycenter) );

        local_velocity(i)= -( element_volume*lumping_weights(i)* (sat_new - sat_old) / dt + local_velocity(i) ) * local_matrix_02(i,i); //
        output_el_sat += lumping_weights(i) * sat_new;
    }

    //cout << "cond diff: " << dh_cell->barycenter()[dim-1] << " "
}

#endif /* LOCAL_ASSEMBLY_HH_ */
