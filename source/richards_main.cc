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
#include <richards_bc.hh>
#include "output.hh"

#define DIM 2

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


/**
 *  Simple test setting for the Richards equation.
 */
class TestProblem{
public:
    TestProblem();
    void solve()
        { equation->run(); }
    ~TestProblem()
        { delete equation; }
private:
    RichardsData<DIM> data;
    Triangulation<DIM>  coarse_tria;
    SmartPointer< Richards_LMH<DIM> > equation;
};

TestProblem::TestProblem ()
: data(h_params)
{
  srand(100);
  // set coarse triangulation
  Point<DIM> pa(0,0), pb(1,-5);
  std::vector<unsigned int> grid_steps;
  grid_steps.push_back(1);
  grid_steps.push_back(5);

  // colorize boudaries (x- 0,1; y- 2,3; z- 4,5) and materials (according to octants)
  GridGenerator::subdivided_hyper_rectangle (coarse_tria,grid_steps ,pa, pb,true);
  //coarse_tria.refine_global(3);

  ParameterHandler prm;
  prm.read_input("input.txt");


  equation = new Richards_LMH<DIM>(coarse_tria,prm,0);
  data.add_bc(0,new BC2D(BC2D::Neuman,0)); // left
  data.add_bc(1,new BC2D(BC2D::Neuman,0)); // right
  //equation->add_bc(2,new BC2D(BC2D::Dirichlet,-150)); // bottom
  data.add_bc(2,new BC2D(BC2D::Neuman,0)); // bottom
  data.add_bc(3,new BC2D(BC2D::Dirichlet,1)); // top
  //equation->add_bc(3,new BC2D(BC2D::Neuman,0)); // top

  data.k_inverse = new KInverse<DIM> (0, pa, pb);
  data.initial_value = new InitialValue<DIM>;
  equation->set_data(&data);
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

