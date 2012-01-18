/* 
 * File:   richards_main.cc
 * Author: jb
 *
 * Created on June 10, 2010, 2:03 PM
 */

#include <grid/grid_generator.h>
#include <stdlib.h>
#include <spatial_functions.hh>
#include <richards_lmh.hh>
//#include <richards_bc.hh>
#include "output.hh"

#define DIM 2

/*
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
*/

/**
 *  Simple test setting for the Richards equation.
 */
class TestProblem{
public:
    TestProblem();
    void declare_params();
    void solve()
        { equation->run(); }
    ~TestProblem() {
        delete equation;
        delete data;
    }
private:
    RichardsData<DIM> *data;
    Triangulation<DIM>  coarse_tria;
    Richards_LMH<DIM> * equation;
    ParameterHandler prm;
};

void TestProblem::declare_params () {

    // space discretization
    prm.declare_entry ("x_size", "1.0",
                        Patterns::Double(),
                        "X size of domain.");
    prm.declare_entry ("z_size", "5.0",
                        Patterns::Double(),
                        "Z size of domain.");
    prm.declare_entry ("hx", "1.0",
                        Patterns::Double(),
                        "X direction space step.");
    prm.declare_entry ("hz", "0.01",
                        Patterns::Double(),
                        "Z direction space step.");

    // time discretization
    prm.declare_entry ("print_time_step", "0.01",
                        Patterns::Double(),
                        "Time step for filed output.");

    prm.declare_entry ("t_init", "0.0",
                        Patterns::Double(),
                        "Initial time.");
    prm.declare_entry ("t_end", "1.0",
                                Patterns::Double(),
                                "End time.");
    prm.declare_entry ("dt_init", "0.0",
                                Patterns::Double(),
                                "Initial dt.");
    prm.declare_entry ("dt_min", "0.0",
                                Patterns::Double(),
                                "Minimal dt.");
    prm.declare_entry ("dt_max", "0.0",
                                Patterns::Double(),
                                "Maximal dt.");
    // physical
    prm.declare_entry("hydro_n","0.0",Patterns::Double(),"");
    prm.declare_entry("hydro_alfa","0.0",Patterns::Double(),"");
    prm.declare_entry("hydro_cut","0.0",Patterns::Double(),"");
    prm.declare_entry("hydro_qr","0.0",Patterns::Double(),"");
    prm.declare_entry("hydro_qs","0.0",Patterns::Double(),"");
    prm.declare_entry("hydro_ks","0.0",Patterns::Double(),"");
    prm.declare_entry("hydro_kk","0.0",Patterns::Double(),"");


    // numerical parameters
    prm.declare_entry("nlin_use_homotopy", "false", Patterns::Bool(),"Use Homotopy method for initial approximation.");
    prm.declare_entry ("nlin_tol", "0.0001",
                                Patterns::Double(),
                                "Tolerance of nonlinear solver.");
    prm.declare_entry ("nlin_tol_step_mult", "100",
                                Patterns::Double(),
                                "Multiple of final tolerance to perform homotopy step.");
    prm.declare_entry ("nlin_tol_max_mult", "10",
                                    Patterns::Double(),
                                    "Multiple of step tolerance to limit length of homotopy step.");

    prm.declare_entry ("nlin_max_lambda", "0.2",
                                Patterns::Double(),
                                "Maximim step in homotopy parameter.");

    prm.declare_entry ("nlin_min_lambda", "0.2",
                                Patterns::Double(),
                                "Minimum (relative) linesearch parameter in Newton method.");

    prm.declare_entry ("nlin_max_it", "20",
                                Patterns::Integer(),
                                "Max nonlin. iterationss.");

    prm.declare_entry ("nlin_alpha", "1.0E-4",
                                Patterns::Double(),
                                "Decrease parameter");

    prm.declare_entry ("nlin_check_jacobian", "1.0E-4",
                                Patterns::Bool(),
                                "Check jacobian by finite diferences.");

    prm.declare_entry ("lin_rtol", "0.01",
                                Patterns::Double(),
                                "Relative tolerance of linear solver.");

    prm.declare_entry ("lin_atol", "1.0E-12",
                                Patterns::Double(),
                                "Absolute tolerance of linear solver.");
    prm.declare_entry ("lin_max_it", "1000",
                                Patterns::Double(),
                                "Max iteration in linear solver.");


    // data
    prm.declare_entry ("hetero_k_n_centers", "0",
                                Patterns::Integer(),
                                "Num of centers of K heterogenity.");


    prm.read_input("input.txt");

}

TestProblem::TestProblem ()

{
  srand(100);

  declare_params();
  /*
  HydrologyParams *hydro_data = new HydrologyParams;
  hydro_data->n = prm.get_double("hydro_n");
  hydro_data->alfa = prm.get_double("hydro_alfa");
  hydro_data->cut_fraction = prm.get_double("hydro_cut");
  hydro_data->Qr = prm.get_double("hydro_qr");
  hydro_data->Qs = prm.get_double("hydro_qs");
  hydro_data->Ks = prm.get_double("hydro_ks");
  hydro_data->Kk = prm.get_double("hydro_kk");
 */
  data = new RichardsData<DIM>(prm);


  // set domain and coarse grid
  double size_x = prm.get_double("x_size"),
         size_z = prm.get_double("z_size");
  Point<DIM> pa(0,0),
             pb(size_x, -size_z);


  double hx = prm.get_double("hx"),
         hz = prm.get_double("hz");
  std::vector<unsigned int> grid_steps;
  grid_steps.push_back(ceil(fabs(size_x)/hx));
  grid_steps.push_back(ceil(fabs(size_z)/hz));

  // colorize boudaries (x- 0,1; y- 2,3; z- 4,5) and materials (according to octants)
  GridGenerator::subdivided_hyper_rectangle (coarse_tria,grid_steps ,pa, pb,true);
  //coarse_tria.refine_global(3);

  equation = new Richards_LMH<DIM>(coarse_tria,prm,0);


  //data.k_inverse = new KInverse<DIM> (0, pa, pb);
  //data.initial_value = new InitialValue<DIM>;
  //data->print_mat_table();
  equation->reinit(data);
}


/**
 * 
 */
int main(int argc, char** argv) {

  //deallog.attach(cout);
  try
    {
//      PetscInitialize(&argc,&argv,0,0);

      deallog.depth_console (0);


      TestProblem problem;
      problem.solve();

 //     PetscFinalize();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

   return (EXIT_SUCCESS);
}

