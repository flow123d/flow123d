#include "output.hh"



///////////////// TODO:
/* - otazka je proc infiltraca pri mixed schamtu trva dele a to i pri jemnejsisch casovych/prostorovych krocich
 * - pro prostorovy krok 0.001 je treba zvysit presnost nelinearniho solveru - dalsi neblaha vlastnost
 *   v pocastu simulace se i pro jemny krok objevily kladne rychlosti
 * - i pro takhle jemnou sit se rychlost chova divne alespon na pocatku, sit stale neni schopna vystihnout 
 *   krivost tlaku blizko bodu nasyceni
 *
 * - mixed schema je nestabilni pro male hodnoty parametru dt*Ks/dx, presneji male hodnoty
 *   tohoto parametru snaze vedou k situaci, kdy skok ve vodivosi je prilis veliky coz vyvola
 *   oscilace rychlosti a tim i tlaku
 *
 *   pro RT0 nejsou oscilace prilis velke, pouze zakmit v tlaku
 *   pro vyssi rad to uplne polozi nelinearni solver
 *
 * - theta time discretization - funguje pro hodnoty theta > 0.5, pro 0.5 dochazi k oscilacim v nasycene elipticke zone
 *   otazka jestli ma pak theta vyznam
 *
 * - dodelat poradne okrajove podminky, zadavani toku, 
 *
 *   napr chybu zpusobenou casovym krokem by slo odovdit z rozdilu implic a explic resp implic a CN schematu
 *   ve vysledku by melo vyjit neco jako dela ten kluk ... viz. simona
 *   .. hm tak ten kluk se jmenuje Kuraz, pro kazdou celociselnou hodnotu tlaku urci rozsah kam az se muze tlak
 *   pohnout aby diference saturaci byla dobre (v dane presnosti) aproximovana kapacitou*diference tlaku
 *   pak voli casovy krok tak, aby tato podminka byla splnena. Jiste je, ze to patrne povede k dobre konvergenci
 *   nelinearniho solveru. Nicmene to nic nerika o casove chybe.
 *
 * - pro prostorovou adaptivitu potrebuju dobry estimator. Kellyho nelze pouzit. Ad hoc napad je diference mezi tlaky, to by se asi pocitalo podobne jako Kelly.
 *   lepsi estimator, ktery by daval skutecnou chybu by byl lepsi.
 *
 * - vykreslovani vtoku v gnuplotu (porovnani oproti jemnemu s1d), druhy indikator
 * - u dalsiho zkoumat jednak chybu a pak pocet linearnich iteraci * pocet dofu
 *
 * - vliv eliminace rychlosti pred vypoctem residua
 * - zkusit adaptivitu v prostoru
 * - rozmyslet kriteria pro volbu casoveho kroku
 * - ma se na zjemnene siti pro ukonceni nelin. solveru pocitat l2 norma vektoru residua
 *   nebo radeji L@ integral z prislusne chyby ve water content
 * - zjistit kde unika voda, asi problem okrajovych podminek
 *
 * 1) pro order > 0 je potreba poradne rozbyslet pocitani saturace z tlaku
 * 2) nejprve odladit 0 order model
 * 3) toky by mely jit pocitat z hodnot bazovych funkci
 * 4) hodilo by se mi, aby konstantni polynomy byly primo v bazi Q_k na elementu i na hranici (pro RT prvky)
 *    pak bych mohl lepe zajistit bilanci
 *
 * adaptivita: viz. tutorial 23, bez Dirichleta pres constrains
 * zacit refine_grid, a zkontrolovat kdy se musi aplikovat constrainy
 */

#define MAX_NUM_OF_DEALII_BOUNDARIES 255

#include <fstream>
#include <iostream>

#include <base/quadrature_lib.h>
#include <base/logstream.h>
//#include <base/function.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/fe_raviart_thomas.h>
#include <numerics/data_out.h>

#include <base/tensor_function.h>
#include <lac/sparse_ilu.h>

 #include <numerics/error_estimator.h>
 #include <grid/grid_refinement.h>
 #include <numerics/solution_transfer.h>

//#include <petsc.h>

using namespace dealii;


//#include <boost/assign.hpp>

//#include <schur_solver.hh> // eliminace rychlosti, zavest jako metodu rozhrani k solveru
//#include "direct_solver.hh"



#include <UMFPACK_system.hh>
#include <time.hh>
#include <spatial_functions.hh>
#include <hydro_functions.hh>
#include <richards_bc.hh>
#include <output.hh>



/**
 *  Richards' equation solver.
 */

static const HydrologyParams h_params=
{
    1.14, //n
    0.1, //alfa
    
    0.01, //Qr
    0.480,  //Qs;
    //1.0,
    0.999,
    //0.99792, //cut_fraction;  simulation is very sensitive for cutting value
                        
    2,  // Ks;             // saturated conductivity
    2   // Kk;             // conductivity at the cut point
};

// nasledujici data jsou z Genuchtenova clanku
static const HydrologyParams sand_stone_params=
{
    10.4, //n
    0.79, //alfa m^-1 (0.0079 cm^-1

    0.153, //Qr
    0.25,  //Qs;
    0.9978, //cut_fraction;

    1.25E-5,  // Ks;             // saturated conductivity m/s (108. cm/day)
    1.25E-5   // Kk;             // conductivity at the cut point
};

static const HydrologyParams silt_loam_GE3_params=
{
    2.06, //n
    0.423, //alfa

    0.131, //Qr
    0.396,  //Qs;
    0.9978, //cut_fraction;

    5.74E-7,  // Ks;             // saturated conductivity
    5.74E-7   // Kk;             // conductivity at the cut point
};


/***
 * rovnice Richards
 * zamyslene vstupy:
 * - hruba sit, pri konstrukci
 * - pocatecni podminka, set .. (i behem vypoctu)
 * - okrajove podminky (std:map pro okrajove indexy hrube site)
 * - hydrologicke funkce  (std:map pro materialove indexy hrube site)
 * - pocatecni a cilovy cas
 *  nepovine:
 * - scaling faktory, anisotropie (Function ..)
 *
 * metody pro pouziti:
 *  1) run - simulace na danem casovem intervalu
 *  2) get_solution, set_solution (pouze tlaky)
 *  3) assembly, .. pouziti nad danou casti globalni matice, v ramci externiho newtonovskeho resice
 *
 */


typedef UMFPACKSystem SolverType;

template <int dim>
class Richards : public Subscriptor {
public:
    Richards(Triangulation<dim> &coarse_tria, unsigned int order);
    //! Add boundary condition bc for the given boundary_index.
    //! boundary without assigned BC are homogeneous Dirichlet, but rather leads to exception
    void add_bc(const unsigned int boundary_index,  BoundaryCondition<dim> *one_bc)
        { Assert( boundary_index < MAX_NUM_OF_DEALII_BOUNDARIES, ExcMessage("invalid index.") );
        bc[boundary_index]=one_bc;
        }

    void add_resistivity(TensorFunction<2,dim> *resist)
        { k_inverse=resist; }

    void run();

private:
    enum block_names {vel_b=0, head_b=1};

    void saturation_update();
    void assemble_system();
    void grid_refine();
    //void compute_errors () const;
    void output_results();
    void output_residuum();

    //! physical data
    SmartPointer< TensorFunction<2,dim> > k_inverse;       //! inverse of (possibly) heterogeneous hydraulic tensor
    //! there should be whole BC descriptor object which returns
    //! BC objects for given index with checking, possibly returning None type of BC
    SmartPointer< BoundaryCondition<dim> > bc[MAX_NUM_OF_DEALII_BOUNDARIES];
    //! TODO: material functions should be packed up into a Soil object
    //! moreover there should be array of Materials for individual mat. indexes
    INV_FK< B<double> > inv_fk_diff;
    FQ< B<double> > fq_diff;
    FC<double> fc;
    FQ<double> fq;
    INV_FK<double> inv_fk;
    //! TODO: precise formulation to include variable density, and gravity field
    double density;
    double gravity;
    int z_component;        // vertical component of the Point<dim>

    //! numerics
    SolverTime time;

    double print_time_step;
    double print_time;
    int print_step;

    int order; // order of FE
    double time_theta; // theta time method
    double weight; // velocity mass lumping weight

    int max_non_lin_iter;
    double nonlin_tol;
    bool newton;



//    std::vector<bool> velocity_component_mask;

    Triangulation<dim> triangulation;
    FESystem<dim> fe;
    FE_DGQ<dim> sat_fe;
    DoFHandler<dim> dof_handler;
    DoFHandler<dim> sat_dh;
    std::vector<unsigned int> component_to_block;

    //! linear algerba objects
    typedef typename SolverType::VectorType VectorType;
    typedef typename SolverType::MatrixType MatrixType;
    SolverType linear_system;

    ConstraintMatrix     hanging_node_constraints;

    //! field vectors
    Vector<double> saturation;
    Vector<double> old_saturation;
    Vector<double> m_mat_crit;
    VectorType solution;
    VectorType old_solution; // previous time step
    VectorType residual;
    Vector<double> theta_diff;

    BCOutput<dim>    bc_out;
};

/**
 *  Construct the solver - discretization space and DOF handler
 */

template <int dim>
Richards<dim>::Richards (Triangulation<dim> &coarse_tria, unsigned int order)
		:
                //bc(MAX_NUM_OF_DEALII_BOUNDARIES),
                inv_fk_diff(h_params),
                fq_diff(h_params),
                fc(h_params),
                fq(h_params),
                inv_fk(h_params),
                time(0.0, 1.0, 0.000001),                  // start, end, dt_init
                order(order),
                fe (FE_RaviartThomas<dim>(order), 1,    // one component is velocity vector
                    FE_DGQ<dim>(order), 1),             // one component is pressure-head, seccond saturation
                sat_fe(order),
		dof_handler (triangulation),             // dof_handler is initialized by empty tria.
                sat_dh(triangulation),
                linear_system(2),                     // three blocks two active
                bc_out("bc_output", order)
{
    density = 1.0;
    gravity = 1.0;
    z_component = dim - 1; // z-coord in Point
    print_time_step = 0.01;

    max_non_lin_iter = 60;
    nonlin_tol = 0.0001;
    //theta=1. / 2.;  // 1 implicit, 0 explicit
    // Crank Nicholson still makes problem with velocity oscilations
    //
    time_theta = 1; //0.5;
    weight=0.0;//0.9; // 1 - mixed. ; 0 - decoupled ?? doesn't work

 

    triangulation.copy_triangulation(coarse_tria);

    dof_handler.distribute_dofs(fe);
    sat_dh.distribute_dofs(sat_fe);
    std::cout << "active cells "<< triangulation.n_active_cells()
	    << " total cells: " << triangulation.n_cells()
            << " dofs: " << dof_handler.n_dofs() << std::endl;

    component_to_block.assign(dim,vel_b);
    component_to_block.push_back(head_b);
    DoFRenumbering::component_wise (dof_handler,component_to_block);
    linear_system.reinit(dof_handler);

    solution.reinit(linear_system.get_rhs());
    old_solution.reinit (linear_system.get_rhs());
    saturation.reinit(solution.block(head_b));
    old_saturation.reinit(solution.block(head_b));
    m_mat_crit.reinit(solution.block(head_b));
    residual.reinit(linear_system.get_rhs());
}


/**
 * Refine grid between time steps, this is far to be optimal, but is simplest
 * since we need not any grid and solution saving.
 * refinement and coarsing shall be done after sucesfull solution of the timestep.
 */
template <int dim>
void Richards<dim>::grid_refine()
{
   Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
   std::vector<bool> component_mask(dim+1, true);
   component_mask[dim]=false;

   // compute refinement indicators
   KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(order+1),
                                      typename FunctionMap<dim>::type(),    // Neuman boundary function
                                      solution,
                                      estimated_error_per_cell,
                                      component_mask);
/*
   { // compute error estimate based on jumps in pressure
    // Main cycle over the cells.
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell != endc; ++cell, ++sat_cell) {
        fe_values.reinit(cell);
        sat_fe_values.reinit(sat_cell);
        local_matrix = 0;
        local_rhs = 0;

    }
    }*/
  /* for(unsigned int i=0; i< estimated_error_per_cell.size(); i++) cout << i << " : " << estimated_error_per_cell(i) <<endl;*/
   /*GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.1, 0.4 , 1000 );*/
   /*GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0, 0.33, 1000 );
    */
    triangulation.prepare_coarsening_and_refinement();
    SolutionTransfer<dim, BlockVector<double> > solution_transfer(dof_handler);
    solution_transfer.prepare_for_coarsening_and_refinement (solution);
    
    SolutionTransfer<dim, Vector<double> > sat_transfer(sat_dh);
    sat_transfer.prepare_for_coarsening_and_refinement (saturation);

    triangulation.execute_coarsening_and_refinement ();
    dof_handler.distribute_dofs (fe);
    sat_dh.distribute_dofs (sat_fe);
        std::cout << "active cells "<< triangulation.n_active_cells()
	    << " total cells: " << triangulation.n_cells()
            << " dofs: " << dof_handler.n_dofs() << std::endl;

    DoFRenumbering::component_wise (dof_handler,component_to_block);
    linear_system.reinit(dof_handler);
    old_solution.reinit (linear_system.get_rhs());
    old_saturation.reinit(old_solution.block(head_b));

    solution_transfer.interpolate (solution, old_solution);
    solution.reinit(linear_system.get_rhs());
    solution=old_solution;
    linear_system.get_constraints().distribute (old_solution);
    linear_system.get_constraints().distribute (solution);

    sat_transfer.interpolate(saturation, old_saturation);
    saturation.reinit(solution.block(head_b));
    // saturation is discontinuous -> no constrains,
    // also we don't copy saturation, since it will be imediatelly updated from solution
    residual.reinit(linear_system.get_rhs());
}


/**
 * update saturation component form pressure head
 */

template <int dim>
void Richards<dim>::saturation_update() {
  for(unsigned int i=0;i<solution.block(head_b).size();i++) {
      //cout << i <<" : " << (solution.block(head_b))(i) << " " << saturation(i) << endl;
      saturation(i)=fq( (solution.block(head_b))(i) );
  }
}
/**
 *  Assembly linear system of one iteration.
 */

template <int dim>
void Richards<dim>::assemble_system() {
    QGauss<dim> quadrature_formula(order + 1);
    QGauss < dim - 1 > face_quadrature_formula(order + 1);

   FEValues<dim> fe_values(fe, quadrature_formula,
            update_values | update_gradients |
            update_quadrature_points | update_JxW_values);
    FEValues<dim> sat_fe_values(sat_fe, quadrature_formula,
            update_values);
    FEFaceValues<dim> fe_face_values(fe, face_quadrature_formula,
            update_values | update_normal_vectors |
            update_quadrature_points | update_JxW_values);


    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> local_rhs(dofs_per_cell);

    std::vector<unsigned int> local_dof_indices(dofs_per_cell);
    std::vector<unsigned int> face_dofs(fe.dofs_per_face);


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
 
    std::map<unsigned int, double> RT_boundary_values;

    linear_system.set_zero();

    double max_m_crit=0;
    // Main cycle over the cells.
    typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
            endc = dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
        sat_cell = sat_dh.begin_active();
    for (; cell != endc; ++cell, ++sat_cell) {
        fe_values.reinit(cell);
        sat_fe_values.reinit(sat_cell);
        local_matrix = 0;
        local_rhs = 0;

        // functions fixed during assembly
        quadrature_points = fe_values.get_quadrature_points();
        right_hand_side.value_list(quadrature_points, rhs_values);
        k_inverse->value_list(quadrature_points, k_inverse_values);

        //fe_values[velocities].get_function_divergences(old_solution, div_u_values);
        fe_values[pressure].get_function_values(solution, pressure_values);
        fe_values[pressure].get_function_values(old_solution, old_pressure_values);
        sat_fe_values.get_function_values(saturation, saturation_values);
        sat_fe_values.get_function_values(old_saturation,saturation_old_values);
        fe_values[velocities].get_function_values(old_solution, old_u_values);
        fe_values[velocities].get_function_divergences(old_solution, old_div_u_values);
        fe_values[velocities].get_function_values(solution, u_values);

        // M-matrix criterium 6*dt > crit.
        double m_crit=0;
        // cycle over quadrature points
        for (unsigned int q = 0; q < n_q_points; ++q) {

            double sat = saturation_values[q];
            double sat_old = saturation_old_values[q];
            double last_p = pressure_values[q];
            double old_p = old_pressure_values[q];
            Tensor < 1, dim> old_u = old_u_values[q];
            double old_div_u = old_div_u_values[q];
            Tensor < 1, dim> velocity = u_values[q];
            B<double> p, f;
            p = (1-time_theta)*old_p + time_theta*last_p;
            f = inv_fk_diff(p);
            f.diff(0, 1);
            double inv_k = f.x();
            double diff_inv_k = p.d(0);
//            cout <<"h: " << last_p << "ik: "<< inv_k << "u: " <<velocity << endl;

            p = last_p;
            f = fq_diff(p);
            f.diff(0, 1);
            //double sat=f.val();
            double diff_sat = p.d(0);
            //double div_u = div_u_values[q];
            //cout << "p "<<last_p << "ik " << inv_k << "dik "<<diff_inv_k<<"ds "<<diff_sat << endl;
            //inv_k=0.01;
            //diff_sat=0.0;
            //sat=sat_old=0.0;

            m_crit+=diff_sat*inv_k*fe_values.JxW(q); // dx*C/K

            for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                const Tensor < 1, dim> phi_i_u = fe_values[velocities].value(i, q);
                const double div_phi_i_u = fe_values[velocities].divergence(i, q);
                const double phi_i_p = fe_values[pressure].value(i, q);

                for (unsigned int j = 0; j < dofs_per_cell; ++j) {
                    const Tensor < 1, dim> phi_j_u = fe_values[velocities].value(j, q);
                    const double div_phi_j_u = fe_values[velocities].divergence(j, q);
                    const double phi_j_p = fe_values[pressure].value(j, q);


                    // posileni diagonaly touto cesotu vede k prilis male rychlosti

                    local_matrix(i,i) +=  (1-weight) * phi_i_u * k_inverse_values[q] *  phi_j_u * inv_k * fe_values.JxW(q);;
                    local_matrix(i, j) += (
                              time_theta* (
                                  weight * phi_i_u * k_inverse_values[q] *  phi_j_u * inv_k
                                - div_phi_i_u * phi_j_p
                                - phi_i_p * div_phi_j_u
                              )
                            - diff_sat / time.dt() * phi_i_p * phi_j_p
                            ) * fe_values.JxW(q);
                    if (newton) {
                        local_matrix(i, j) += phi_i_u * k_inverse_values[q] * velocity * diff_inv_k * phi_j_p;
                    }

                }

                local_rhs(i) += (
                         - phi_i_p * rhs_values[q]
                         + density * gravity
                           * quadrature_points[q](z_component) * div_phi_i_u
                         -(1-time_theta)*(
                            phi_i_u* k_inverse_values[q] * inv_k * old_u
                            - div_phi_i_u * old_p
                            - phi_i_p * old_div_u
                          )
                         - diff_sat / time.dt() * phi_i_p * last_p
                         + (sat - sat_old) / time.dt() * phi_i_p
                         ) * fe_values.JxW(q);
                if (newton) {
                    local_rhs(i) += 
                            phi_i_u * k_inverse_values[q] * velocity * diff_inv_k * last_p;
                }
            }
        } // q points cycle

        m_crit*=cell->measure();

        if (m_crit> max_m_crit) max_m_crit=m_crit;
       // if (m_crit > 6*time.dt()) {
       //     cout << "CRIT: "<< m_crit/(6*time.dt()) << ", " << m_crit/cell->measure() << endl;
       // }
        std::vector<unsigned int> sat_local_dofs(sat_fe.dofs_per_cell);
        sat_cell->get_dof_indices(sat_local_dofs);
        for(int i=0;i<sat_local_dofs.size();++i) m_mat_crit(sat_local_dofs[i])=m_crit/(6*time.dt());

        if (cell->at_boundary()) {
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no) {
                unsigned int b_ind = cell->face(face_no)->boundary_indicator();
                if (b_ind == 255) continue;
                BoundaryCondition<dim> *face_bc = bc[b_ind];
                fe_face_values.reinit(cell, face_no);

                //cout << "BC:" << b_ind << " "<<face_bc->type() <<endl;
                switch (face_bc->type()) {
                    case BoundaryCondition<dim>::Dirichlet:
                        //! add boundary term
                        //! \f$ \int_{\prtl \Omega} V\cdot N ( h_D +z) \f$
                        quadrature_points = fe_face_values.get_quadrature_points();
                        for (unsigned int q = 0; q < n_face_q_points; ++q) {
                            Point<dim> &q_point=quadrature_points[q];
                            //cout << "val:" << face_bc->value(q_point)<<endl;
                            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                                local_rhs(i) += -(fe_face_values[velocities].value(i, q) *
                                    fe_face_values.normal_vector(q) *
                                    (face_bc->value(q_point) + density * gravity * q_point(z_component)) *
                                    fe_face_values.JxW(q));

                        }


                        break;
                        
                    case BoundaryCondition<dim>::Neuman:
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
                    default:
                         Assert (false, ExcNotImplemented());
                         break;
                }
            }
        }

        cout << "------------\n";
        local_matrix.print_formatted(cout);
        cell->get_dof_indices(local_dof_indices);
        linear_system.matrix_add(local_dof_indices, local_matrix);
        linear_system.rhs_add(local_dof_indices, local_rhs);


    }

    MatrixTools::apply_boundary_values(RT_boundary_values,
            linear_system.get_matrix(),
            solution,
            linear_system.get_rhs());


     cout << "CRIT: "<< max_m_crit/(6*time.dt()) << endl;
}


/**
 *  Output method
 */

template <int dim>
void Richards<dim>::output_results ()
{
  std::vector<std::string> solution_names;

  bc_out.output(time, dof_handler, sat_dh, solution, saturation);

 //if (time.t() < print_time) return;
 cout << "PRINT time (" << print_step << "): " << time.t() << endl;
 print_time=time.t() + print_time_step;
 time.add_target_time(print_time);

  switch (dim)
    {
      case 2:
            solution_names.push_back ("u");
            solution_names.push_back ("v");
            solution_names.push_back ("p");
            //solution_names.push_back ("r");

            break;
            
      case 3:
            solution_names.push_back ("u");
            solution_names.push_back ("v");
            solution_names.push_back ("w");
            solution_names.push_back ("p");
//            solution_names.push_back ("S");

            break;
            
      default:
            Assert (false, ExcNotImplemented());
            break;
    }

  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, solution_names);
  data_out.build_patches (order+1);

  std::ostringstream filename;
  filename << "./output/solution-" << setfill('0')<<setw(3) << print_step << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
  
  print_step++;
}


/**
 *  Output residuum and other numericaly important fields during given interval.
 * 
 */ 

template <int dim>
void Richards<dim>::output_residuum()
{
  static int n_output=0;
 
  std::vector<std::string> solution_names;
  solution_names.push_back ("residuum");
  
  DataOut<dim> data_out;
  data_out.attach_dof_handler (sat_dh);
  data_out.add_data_vector (residual.block(head_b), "residuum");
  data_out.add_data_vector (old_saturation, "old_saturation");
  data_out.add_data_vector (saturation, "saturation");
  //data_out.add_data_vector (solution.block(head_b), "crit");
  data_out.add_data_vector (m_mat_crit, "crit");
  data_out.build_patches (order+1);

  std::ostringstream filename; 
  filename << "./output/residuum-" << setfill('0')<<setw(3) << n_output << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
  
  n_output++;    
}


/**
 *   Computation method.
 *
 */


template <int dim>
void Richards<dim>::run () 
{

//    VectorType &system_rhs=linear_system->get_rhs();
//    MatrixType &system_matrix=linear_system->get_matrix();

//    Vector<double> rhs_vel(solution.block(1).size());

    Vector<double> head_last, sat_last;
    Vector<double> head_new;

//    SparseILU<double> preconditioner;

    double res_norm;
    double last_res_norm;
    double decrease, last_decrease;

  {
    ConstraintMatrix constraints;
    constraints.close();
    
    VectorTools::project (dof_handler,
                          constraints,
                          QGauss<dim>(order+2),
                          InitialValues<dim>(),
                          solution);
  }

  head_last.reinit(old_solution.block(head_b).size());
  head_new.reinit(old_solution.block(head_b).size());
  sat_last.reinit(old_solution.block(head_b).size());
  // update saturation
  //saturation=1.;
  //old_saturation=1.;
  saturation_update();
  print_step=0.01;
  print_time=time.t();
  output_results ();
  time.inc();

  double lambda;
  int iter;
  int cum_iter=0;

  do {  // ---------------------------- Time loop
      old_saturation=saturation;
      old_solution=solution;
      grid_refine();
renew_timestep:
      std::cout << "Timestep " << time.n_step() << " t= "<<time.t()<<", dt= " << time.dt() << endl;

      iter=0;
      last_decrease=decrease=0.5;
      last_res_norm=100;
      newton=false;
      do {
            lambda=1;
 line_search:
            sat_last=saturation;
            saturation_update();
        
            assemble_system (); // this also apply boundary conditions thus change solution
            //std::cout << "assembled, ";


            // we take water inbalance / dt as a residual
            // this is harder to meet for shor timesteps, but
            // but leads to error independent of timestep
            // possibly it is to hard to satify ???
            //linear_system.eliminate_block(solution,vel_b);
            res_norm=linear_system.residual(solution,residual);
            res_norm*=time.dt();
            //output_residuum();

            last_decrease=decrease;
            decrease = (iter == 0? 100 : res_norm/last_res_norm);
            
            {
            // vypada to, ze posun cela by mohla lepe zachytit norma zmeny reseni - integral
            Vector<double> diff(saturation);
            Vector<double> diff_integral(triangulation.n_active_cells());
            
            VectorTools::integrate_difference(sat_dh, residual.block(head_b), ZeroFunction<dim>(1),diff_integral ,QGauss<dim>(order+1),VectorTools::L2_norm);
            double res_diff=diff_integral.l2_norm();
            
            diff=saturation; diff-=sat_last;
            VectorTools::integrate_difference(sat_dh, diff, ZeroFunction<dim>(1),diff_integral ,QGauss<dim>(order+1),VectorTools::L2_norm);
            double sat_diff=diff_integral.l2_norm();

            diff=solution.block(head_b); diff-=head_last;
            VectorTools::integrate_difference(sat_dh, diff, ZeroFunction<dim>(1),diff_integral ,QGauss<dim>(order+1),VectorTools::L2_norm);
            double head_diff=diff_integral.l2_norm();

            //cout<< "res norm "<<res_norm << " decrease: "<< decrease <<endl;
            cout<< iter<<" (r,dic,dr,ds,dh): "
                <<res_norm <<" "
                <<decrease <<" "
                <<res_diff <<" "
                <<sat_diff <<" "
                <<head_diff << " " << endl;

            }
            // this is very simple and ugly linesearch the aim is only
            // demonstrate that the equation have solution even for longer timesteps
            // so the adaptivity should not be used to overcome bed nonlinear solver
/*
            if (iter>0 && decrease < 1.1 0.9 && lambda > 0.00001 ) {
                lambda *= 1.0 / 2.0;
                cout << "lambda: "<< lambda << endl;
                solution.block(0).sadd(0.0, lambda, head_new, (1-lambda), head_last);

                goto line_search;
            }*/
            //if (iter > 4) newton =true;
            // test convergency
            // we require at least one iteration
            // if we only test the residuum, it is usualy very low at the beginning of infiltration
            // especialy with small timesteps, the error is hidden in the velocity
            // maybe we should check both the residuum before velocity recovery and after
            if (iter>0 && res_norm < nonlin_tol) break; // ------------------------------------------------------

            if (iter > max_non_lin_iter /*|| (iter>10 && decrease > 0.9 && last_decrease > 0.9 )*/ ) {
                //cout << "!!! divergence of nonlinear solver" << endl;

                solution=old_solution;
                saturation=old_saturation;
                time.reinc_time(0.7);

                goto renew_timestep;
            }



            linear_system.solve(solution);
            head_last=head_new;
            head_new=solution.block(head_b);

            iter++;
            cum_iter++;
            last_res_norm=res_norm;
      } while (1);

      // time error estimator
      if (theta_diff.size() == 0) {
        theta_diff.reinit(saturation.size(),false);
      }

      double diff;
      double err=0;
      for(int i=0;i<saturation.size();++i) {
          diff=(saturation(i)-old_saturation(i))/time.dt();
          //cout << diff <<endl;
          err=std::max(err, diff - theta_diff(i) );
          theta_diff(i) = diff;
      }
      err*=time.dt()/2.0;

      cout << "time err est: " << err << endl;


      // adapt time step
      double factor=1.0;
      if (iter<5) factor=1.0 / 0.9;
      else if (iter>13) factor=0.9;
      time.scale_time_step(factor);

      std::cout << "nonlin iter: "<< iter << " cum: "<<cum_iter<< endl;

      double final_res_norm;
      linear_system.eliminate_block(solution,vel_b);
      final_res_norm=linear_system.residual(solution,residual);
      cout << "final res: " << final_res_norm << endl;
      residual.block(head_b)*=time.dt();
      saturation-= residual.block(head_b);
      
      output_results ();
      output_residuum();
      
      std::cout << endl;
                
    }
  while (time.inc());
}
