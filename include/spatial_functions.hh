/* 
 * File:   spatial_functions.hh
 * Author: jb
 *
 * Created on May 17, 2010, 1:50 PM
 */

#ifndef _SPATIAL_FUNCTIONS_HH
#define	_SPATIAL_FUNCTIONS_HH

#include <iostream>

#include <base/function.h>
#include <base/tensor_function.h>
#include <base/tensor_base.h>
#include <base/tensor.h>

#include <hydro_functions.hh>
#include <richards_bc.hh>
#include "FADBAD++/fadbad.h"
#include "FADBAD++/badiff.h"
using namespace fadbad;

//*****************************************************************
// FUNCTIONS
/**
 *  Water source function
 */
using namespace dealii;


template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    virtual ~RightHandSide() {}
};

template <int dim>
double RightHandSide<dim>::value (const Point<dim>  &/*p*/,
				  const unsigned int /*component*/) const
{
  return 0;
}

/**
 *  Dirichlet pressure head on the boundary.
 */


template <int dim>
class PressureBoundaryValues : public Function<dim>
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

template <int dim>
double PressureBoundaryValues<dim>::value (const Point<dim>  &p,
					   const unsigned int component) const
{

    //return ( (1-p(0)) * (1+p(0)) * (p(1)+1) ); // just constant saturated
    //return ( (p(1)+1)*(p(1)+1)*11.0 / 2.0 -20 ); // just constant saturated
    return 0.5;
}


/*
template <int dim>
class ExactSolution : public Function<dim>
{
  public:
    ExactSolution () : Function<dim>(dim+1) {}

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &value) const;
};








template <int dim>
void
ExactSolution<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &values) const
{
  Assert (values.size() == dim+1,
	  ExcDimensionMismatch (values.size(), dim+1));

  const double alpha = 0.3;
  const double beta = 1;

  values(0) = alpha*p[1]*p[1]/2 + beta - alpha*p[0]*p[0]/2;
  values(1) = alpha*p[0]*p[1];
  values(2) = -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
}
*/

//*************************************************************************
/**
 * KInverse - class that provides possibly space variable conductivity tensor.
 */


template <int dim>
class KInverse : public TensorFunction<2,dim>
{
  public:
    KInverse (const unsigned int n_centers, Point<dim> pa, Point<dim> pb);

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<2,dim> >    &values) const;

  private:
    std::vector<Point<dim> > centers;
    double radius;
};


template <int dim>
KInverse<dim>::KInverse (const unsigned int n_centers, Point<dim> pa, Point<dim> pb)
{
  centers.resize (n_centers);
  for (unsigned int i=0; i<n_centers; ++i)
    for (unsigned int d=0; d<dim; ++d) {
      centers[i][d] = pa[d] + ( ((double)rand()) * (pb[d]-pa[d]) ) / RAND_MAX;
      std::cout << centers[i][d] << " ";
    }
  radius = 1.0E100;
  for (unsigned int d=0; d<dim; ++d)
    radius =   std::min(radius, std::abs( pb[d] - pa[d] ) );

  std::cout << std::endl << "rad: "<< radius<< std::endl;
}


template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
          ExcDimensionMismatch (points.size(), values.size()));

  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();

      double permeability = 0;
      if (centers.size() > 0) {
        for (unsigned int i=0; i<centers.size(); ++i)
          permeability += std::exp(-(points[p]-centers[i]).square() / radius / radius ) + 0.5; // minimal relative permeability

        permeability= std::max(permeability, 0.005);
      } else {
          permeability= 1.0;
      }
      for (unsigned int d=0; d<dim; ++d)
        values[p][d][d] = 1./permeability;
    }
}

/**
 *  initial pressure head profile
 *
 */

template <int dim>
class InitialValue : public Function<dim>
{

  public:
    InitialValue () : Function<dim>(1) {}

    virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const
    {
        //return -100 + p[dim-1]; // piezometric head
        // analytical solution
        double s = p[dim-1];
        return -tan( (exp(s)-1) / (exp(s) +1 ));

    }


};

#define MAX_NUM_OF_DEALII_BOUNDARIES 255

template <int dim>
struct RichardsData {
public:

    RichardsData(const HydrologyParams & h_params):
    fq_diff(h_params),
    //fc(h_params),
    fq(h_params),
    inv_fk(h_params),
    inv_fk_diff(h_params),
    density(1.0),
    gravity(0.0),
    bc(MAX_NUM_OF_DEALII_BOUNDARIES, NULL)
    {}

    ~RichardsData() {
        if (k_inverse != NULL) delete k_inverse;
        if (initial_value != NULL) delete initial_value;
        for(unsigned int i=0; i<MAX_NUM_OF_DEALII_BOUNDARIES;i++)
            if (bc[i] !=NULL) delete bc[i];
    }

    void add_bc(const unsigned int boundary_index,  BoundaryCondition<dim> *one_bc)
        { Assert( boundary_index < MAX_NUM_OF_DEALII_BOUNDARIES, ExcMessage("invalid index.") );
        bc[boundary_index]=one_bc;
        }

    inline double pressure(double p_head, Point<dim> point)
    {
        return (p_head - point[dim-1])*density*gravity;
    }

    void print_mat_table() {
        cout << "MATERIAL TABLE" <<
        cout << "(h, sat, cond)" << endl;
        for(double h = -50; h < 0.1; h+=0.1) {
            cout << h << " " << fq(h) << " " << 1.0/inv_fk(h) << endl;
        }
    }

    KInverse<dim> *k_inverse;
    InitialValue<dim> *initial_value;

    //! there should be whole BC descriptor object which returns
    //! BC objects for given index with checking, possibly returning None type of BC
    std::vector< BoundaryCondition<dim> *> bc;

    FQ_analytical< B<double> > fq_diff;
    //FC<double> fc;
    FQ_analytical<double> fq;
    INV_FK_analytical<double> inv_fk;
    INV_FK_analytical< B<double> > inv_fk_diff;

    double density;
    double gravity;
};



#endif	/* _SPATIAL_FUNCTIONS_HH */

