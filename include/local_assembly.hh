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
    void set_vectors(const Vector<double> &sol,const  Vector<double> &old_sol,const Vector<double> &sat,const Vector<double> &old_sat)
    {
        solution=&sol;
        old_solution=&old_sol;
        saturation=&sat;
        old_saturation=&old_sat;
    }

    void set_data(struct RichardsData<dim> * data)
    {
        richards_data=data;
    }

    void set_dt(const double in_dt) {
        dt = in_dt;
    }

    void reinit(typename DoFHandler<dim>::active_cell_iterator cell);
    FullMatrix<double> &get_matrix() {return local_matrix_22;}
    Vector<double> &get_rhs() {return local_rhs_2;}

private:
    QGauss<dim> quadrature_formula;
    QGauss < dim - 1 > face_quadrature_formula;
    const unsigned int n_q_points;
    const unsigned int n_face_q_points;

    FE_DGVector<PolynomialsRaviartThomas<dim>, dim> velocity_fe;
    FE_FaceQ<dim> trace_fe;
    FE_DGQ<dim> pressure_fe;

    FEValues<dim> pressure_fe_values;
    FEValues<dim> velocity_fe_values;
    FEFaceValues<dim> trace_fe_face_values;
    FEFaceValues<dim> velocity_fe_face_values;
    const FEValuesExtractors::Vector vec_extr;
    const FEValuesExtractors::Scalar scal_extr;


    FullMatrix<double> local_matrix_00, local_matrix_01, local_matrix_02, local_matrix_22, inv_local_matrix_00;
    Vector<double> local_rhs_0,local_rhs_2, alphas, lumping_weights, ones_vec_0;
    double alphas_total;

    std::map<unsigned int, double> RT_boundary_values; // Dirichlet boundary values

    const Vector<double> * solution;
    const Vector<double> * old_solution;
    const Vector<double> * saturation;
    const Vector<double> * old_saturation;
    RichardsData<dim> *richards_data;

    std::vector<unsigned int> local_dof_indices;
    std::vector<unsigned int> face_dofs;

    const unsigned int z_component;
    const double density, gravity;
    double dt;

    std::vector<Tensor < 2, dim> > k_inverse_values;


/*
    const RightHandSide<dim> right_hand_side;
    const PressureBoundaryValues<dim> pressure_boundary_values;


    std::vector<Point<dim> > quadrature_points;
    std::vector<double> rhs_values(n_q_points);
    std::vector<double> pressure_values(n_q_points);
    std::vector<double> old_pressure_values(n_q_points);
    std::vector<double> saturation_old_values(n_q_points);
    std::vector<double> saturation_values(n_q_points);
    std::vector<Tensor <1,dim> > old_u_values(n_q_points);
    std::vector<Tensor <1,dim> > u_values(n_q_points);
    std::vector<double> old_div_u_values(n_q_points);
    std::vector<Tensor <1,dim> > face_u_values(n_face_q_points);
    std::vector<double> face_p_values(n_face_q_points);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    const FEValuesExtractors::Scalar pressure_trace(dim+1);

    std::map<unsigned int, double> RT_boundary_values;
*/
};

template <int dim>
LocalAssembly<dim>::LocalAssembly(unsigned int order)
:
quadrature_formula(order+1),
face_quadrature_formula(order + 1),
n_q_points(quadrature_formula.size()),
n_face_q_points(face_quadrature_formula.size()),

velocity_fe(order,mapping_piola),
pressure_fe(order),
trace_fe(order),

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


local_rhs_0(velocity_fe.dofs_per_cell),
local_rhs_2(trace_fe.dofs_per_cell),
alphas(velocity_fe.dofs_per_cell),
lumping_weights(velocity_fe.dofs_per_cell),
ones_vec_0(velocity_fe.dofs_per_cell),

density(1.0),
gravity(1.0),
z_component(dim-1),

k_inverse_values(n_q_points)

/*
const RightHandSide<dim> right_hand_side;
const PressureBoundaryValues<dim> pressure_boundary_values;


std::vector<Point<dim> > quadrature_points;
std::vector<double> rhs_values(n_q_points);
std::vector<Tensor < 2, dim> > k_inverse_values(n_q_points);
std::vector<double> pressure_values(n_q_points);
std::vector<double> old_pressure_values(n_q_points);
std::vector<double> saturation_old_values(n_q_points);
std::vector<double> saturation_values(n_q_points);
std::vector<Tensor <1,dim> > old_u_values(n_q_points);
std::vector<Tensor <1,dim> > u_values(n_q_points);
std::vector<double> old_div_u_values(n_q_points);
std::vector<Tensor <1,dim> > face_u_values(n_face_q_points);
std::vector<double> face_p_values(n_face_q_points);

const FEValuesExtractors::Vector velocities(0);
const FEValuesExtractors::Scalar pressure(dim);
const FEValuesExtractors::Scalar pressure_trace(dim+1);

std::map<unsigned int, double> RT_boundary_values;
*/
{
    ones_vec_0=1.0;
}

template <int dim>
void LocalAssembly<dim>::reinit(typename DoFHandler<dim>::active_cell_iterator cell)
{

    local_matrix_00=0;
    local_matrix_01=0;
    local_matrix_02=0;
    local_matrix_22=0;

    local_rhs_0=0;
    local_rhs_2=0;

    typename Triangulation<dim,dim>::active_cell_iterator tria_cell(cell);
    pressure_fe_values.reinit(tria_cell);
    velocity_fe_values.reinit(tria_cell);

    richards_data->k_inverse->value_list(velocity_fe_values.get_quadrature_points(), k_inverse_values);

    std::cout << "initialized" << std::endl;
    // assembly velocity and cell pressure blocks
    // cycle over quadrature points
    for (unsigned int q = 0; q < n_q_points; ++q) {
        for (unsigned int i = 0; i < velocity_fe.dofs_per_cell; ++i) {
            const Tensor < 1, dim> phi_i_u = velocity_fe_values[vec_extr].value(i, q);
            const double div_phi_i_u = velocity_fe_values[vec_extr].divergence(i, q);

            // matrix_00, velocity * velocity
            for (unsigned int j = 0; j < velocity_fe.dofs_per_cell; ++j) {
                const Tensor < 1, dim> phi_j_u = velocity_fe_values[vec_extr].value(j, q);

                local_matrix_00(i,j) += (phi_i_u * k_inverse_values[q] *  phi_j_u) * velocity_fe_values.JxW(q);
            }

            // matrix 01, velocity * element pressure
            for (unsigned int j = 0; j < pressure_fe.dofs_per_cell; ++j) {
                const double phi_j_p = pressure_fe_values[scal_extr].value(j, q);

                local_matrix_01(i,j) += - div_phi_i_u * phi_j_p * velocity_fe_values.JxW(q);
            }

            // rhs 0, piezometric head
            local_rhs_0(i) += density * gravity
              * ( (velocity_fe_values.get_quadrature_points())[q])(z_component)
              * div_phi_i_u * velocity_fe_values.JxW(q);

        }
    }
    // lm_00 inverse and compute lumping weights
    inv_local_matrix_00.invert(local_matrix_00);
    inv_local_matrix_00.vmult(alphas, ones_vec_0);
    alphas_total = alphas * ones_vec_0;
    lumping_weights = alphas;
    lumping_weights.scale(1/alphas_total);

    const double element_volume = cell->measure();

    // matrix 22, trace * trace - lumped diagonal matrix for time term
    // not clear how to write this as cycle over quadrature points, namely i it not clear
    // what are base functions for boundary dofs and how to extend this to higher dimensions
    //
    // but form face_fe_values we can extract only values in quadrature points

    fadbad::B<double> p,f;

    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < trace_fe.dofs_per_cell; ++i) {


                   double sat_new = (*saturation)(local_dof_indices[i]);
                   double sat_old = (*old_saturation)(local_dof_indices[i]);
                   double trace_pressure = (*solution)(local_dof_indices[i]);
                   p = trace_pressure;
                   f = richards_data->fq_diff(p);
                   f.diff(0, 1);
                   //double sat=f.val();
                   double capacity = p.d(0);

                   local_matrix_22(i,i)=lumping_weights(i) * element_volume / dt * capacity;
                   local_rhs_2(i)=lumping_weights(i) * element_volume / dt * (
                           capacity * trace_pressure  - (sat_new - sat_old));

    }

    local_matrix_22.add(1.0,inv_local_matrix_00);
    for(unsigned int i=0; i < local_matrix_22.size(0); ++i) {
        for(unsigned int j=0; j < local_matrix_22.size(0); ++j) {
            local_matrix_22(i,j) += - alphas(i)*alphas(j)/alphas_total;
        }
    }



    // assembly boundary terms
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
        trace_fe_face_values.reinit(cell, face_no);

        // 02 matrix should be just identity matrix
        for (unsigned int q = 0; q < n_face_q_points; ++q) {
            for (unsigned int i = 0; i < velocity_fe.dofs_per_cell; ++i) {
                const Tensor < 1, dim> phi_i_u = velocity_fe_values[vec_extr].value(i, q);

                // matrix 02, velocity * trace pressure
                for (unsigned int j = 0; j < trace_fe.dofs_per_cell; ++j) {
                    const double phi_j_p = trace_fe_face_values[scal_extr].value(j, q);

                    local_matrix_02(i,j) += trace_fe_face_values.normal_vector(q) * phi_i_u *
                            phi_j_p * trace_fe_face_values.JxW(q);
                    Assert ( (i==j) || local_matrix_02(i,j)==0.0,
                             ExcDimensionMismatch (i, j));

                }
            }
        }

        // apply boundary conditions
        unsigned int b_ind = cell->face(face_no)->boundary_indicator();
        if (b_ind == 255) continue;
        BoundaryCondition<dim> *face_bc = richards_data->bc[b_ind];

         //cout << "BC:" << b_ind << " "<<face_bc->type() <<endl;
         switch (face_bc->type()) {
             case BoundaryCondition<dim>::Dirichlet:
             /*
                    //! add boundary term
                    //! \f$ \int_{\prtl \Omega} V\cdot N ( h_D +z) \f$
                    quadrature_points = fe_face_values.get_quadrature_points();
                    for (unsigned int q = 0; q < n_face_q_points; ++q) {
                        Point<dim> &q_point=quadrature_points[q];
                        //cout << "val:" << face_bc->value(q_point)<<endl;
                        for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                            local_rhs(i) += -(fe_face_values[velocities].value(i, q) *
                                fe_face_values.normal_vector(q) *
                                (face_bc->value(q_point) + density * gravity * q_point(z_component)) *
                                fe_face_values.JxW(q));

                    }
*/

                    break;

             case BoundaryCondition<dim>::Neuman:
             /*
                    // up to now only homogeneous, and not sure it is OK
                    // set all dofs on boundary to zero
                    //
                    // questions:
                    // 1) Are the basic functions with nodal points on the boundary (which means taht with boundary moments)
                    //    the only functions with non zero normal flux over the boundary?
                    //    Or also interior functions have nonzero normal flux?
                    // 2) Can I use FETools::inteploate  to set Dirichlet DOFs, giving
                    // a vector function in the quadrature points?

                    cell->face(face_no)->get_dof_indices(face_dofs, cell->active_fe_index());

                    for (unsigned int i = 0; i < face_dofs.size(); ++i) {
                        if (fe.n_nonzero_components(fe.face_to_cell_index(i, face_no)) == dim) {
                            // for FESystem with only one RT comonent this
                            // this is enough to identify it, otherwise
                            // one have to use get_nonzero_components
                            // and analyse its contents

                            RT_boundary_values[face_dofs[i]] = 0.0;

                        };
                    }
                    break;
                    */
                default:
                     Assert (false, ExcNotImplemented());
            }


    }


    // apply constrains

}

#endif /* LOCAL_ASSEMBLY_HH_ */
