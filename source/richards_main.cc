/* 
 * File:   richards_main.cc
 * Author: jb
 *
 * Created on June 10, 2010, 2:03 PM
 */

#include <grid/grid_generator.h>
#include <stdlib.h>
#include <richards.hh>

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
    Triangulation<2>  coarse_tria;
    SmartPointer< Richards<2> > equation;
};

TestProblem::TestProblem ()
{
  srand(100);
  // set coarse triangulation
  Point<2> pa(0,0), pb(1,-5);
  std::vector<unsigned int> grid_steps;
  grid_steps.push_back(1);
  grid_steps.push_back(5000);

  // colorize boudaries (x- 0,1; y- 2,3; z- 4,5) and materials (according to octants)
  GridGenerator::subdivided_hyper_rectangle (coarse_tria,grid_steps ,pa, pb,true);
  //coarse_tria.refine_global(3);


  
  equation = new Richards<2>(coarse_tria,0);
  equation->add_bc(0,new BC2D(BC2D::Neuman,0)); // left
  equation->add_bc(1,new BC2D(BC2D::Neuman,0)); // right
  //equation->add_bc(2,new BC2D(BC2D::Dirichlet,-150)); // bottom
  equation->add_bc(2,new BC2D(BC2D::Neuman,0)); // bottom
  equation->add_bc(3,new BC2D(BC2D::Dirichlet,1)); // top
  //equation->add_bc(3,new BC2D(BC2D::Neuman,0)); // top
  equation->add_resistivity(new KInverse<2> (0, pa, pb));
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

