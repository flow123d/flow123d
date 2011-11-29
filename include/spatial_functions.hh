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
//#include <richards_bc.hh>
#include "FADBAD++/fadbad.h"
#include "FADBAD++/badiff.h"
#include "base/parameter_handler.h"
using namespace fadbad;
using namespace dealii;

//*****************************************************************
// FUNCTIONS
//*****************************************************************


/**
 *  Water source function
 */
template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>(1) {}
    virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const
    {
        return 0;
    }
    virtual ~RightHandSide() {}
};

/**
 *  BOUNDARY FUNCTIONS
 */

template <int dim>
class ASol_ATan : public Function<dim> {

    virtual double value (const Point<dim>   &p, const unsigned int component = 0) const
    {
        double s = p[dim-1] - this->get_time();
        if (s <=0) return -s/2;
        else return -tan( (exp(s)-1) / (exp(s) +1 ));
    }
    virtual ~ASol_ATan() {}
};

template <int dim>
class AFlux_ATan : public Function<dim> {

    virtual double value (const Point<dim>   &p, const unsigned int component = 0) const
    {
        double s = p[dim-1] - this->get_time();
        if (s <=0) return 0.5*2.0;
        else {
            double tmp=cos( (exp(s)-1)/(exp(s)+1) ) * ( exp(s)+1 );
            double h=-tan( (exp(s)-1) / (exp(s) +1 ));
            return 2*exp(s)/tmp/tmp *2.0 / (1+h*h);
        }
    }
    virtual ~AFlux_ATan() {}
};

template <int dim>
class ASol_lin : public Function<dim> {

    virtual double value (const Point<dim>   &p, const unsigned int component = 0) const
    {
        double s = p[dim-1] - this->get_time();
        return  exp(-s);
    }
    virtual ~ASol_lin() {}
};

template <int dim>
class ASol_lin_steady : public Function<dim> {

    virtual double value (const Point<dim>   &p, const unsigned int component = 0) const
    {
        double s = p[dim-1];
        return exp(s/2.0) + exp(-s/2.0);
    }
    virtual ~ASol_lin_steady() {}
};


//*************************************************************************
/**
 * KInverse - class that provides possibly space variable conductivity tensor.
 */


template <int dim>
class KInverse : public TensorFunction<2,dim>
{
  public:
    KInverse (ParameterHandler &prm);

    virtual void value_list (const std::vector<Point<dim> > &points,
                             std::vector<Tensor<2,dim> >    &values) const;

  private:
    std::vector<Point<dim> > centers;
    double radius;
};


template <int dim>
KInverse<dim>::KInverse (ParameterHandler &prm)
{
  unsigned int n_centers = prm.get_integer("hetero_k_n_centers");
  double size[dim];
  size[0] = prm.get_double("x_size");
  size[1] = prm.get_double("z_size");

  centers.resize (n_centers);
  for (unsigned int i=0; i<n_centers; ++i)
    for (unsigned int d=0; d<dim; ++d) {
      centers[i][d] = ( ((double)rand()) * size[d] ) / RAND_MAX;
      std::cout << centers[i][d] << " ";
    }
  radius = 1.0E100;
  for (unsigned int d=0; d<dim; ++d)
    radius =   std::min(radius, std::abs( size[d] ) );

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


#define MAX_NUM_OF_DEALII_BOUNDARIES 255

/**
 * Collect all equation data.
 */
template <int dim>
class RichardsData {
public:
    typedef enum { None, Dirichlet, Neuman } BCType;
    KInverse<dim> k_inverse;
    Function<dim> *initial;
    Function<dim> *anal_sol;
    Function<dim> *anal_flux;

private:
    //! there should be whole BC descriptor object which returns
    //! BC objects for given index with checking, possibly returning None type of BC
    std::vector< std::pair<BCType, Function<dim> *> > bc_segments;

    //FQ_lin< B<double> > fq_diff;
    //FK_lin< B<double> > fk_diff;

    FQ_analytical< B<double> > fq_diff;
    FK_analytical< B<double> > fk_diff;

    double density;
    double gravity;

public:

    RichardsData(const HydrologyParams & h_params, ParameterHandler &prm):
    k_inverse(prm),
    initial(NULL),
    bc_segments(MAX_NUM_OF_DEALII_BOUNDARIES, std::pair<BCType, Function<dim> * >(None, (Function<dim> *)(NULL) ) ),
    fq_diff(h_params),
    //fc(h_params),
    fk_diff(h_params),
    density(1.0),
    gravity(0.0)
    {
        add_bc(0, Neuman, new ZeroFunction<dim>() ); // left
        add_bc(1, Neuman, new ZeroFunction<dim>() ); // right
        add_bc(2, Dirichlet, new ASol_ATan<dim>() ); // bottom
        add_bc(3, Dirichlet, new ASol_ATan<dim>()); // top
        //add_bc(3,Dirichlet, new ConstantFunction<dim>(1.0); // top
        //add_bc(3,Dirichlet, new ConstantFunction<dim>(1.0); // top

        initial = new ASol_ATan<dim>();
        anal_sol = new ASol_ATan<dim>();
        anal_flux = new AFlux_ATan<dim>();

        //anal_sol = new ASol_lin<dim>();
        //initial = new ASol_lin<dim>();


    }

    ~RichardsData() {
        for(unsigned int i=0; i<MAX_NUM_OF_DEALII_BOUNDARIES;i++)
            if (bc_segments[i].second != NULL) delete bc_segments[i].second;
    }

    void add_bc(const unsigned int boundary_index, const BCType bt, Function<dim> * f)
    {
        Assert( boundary_index < MAX_NUM_OF_DEALII_BOUNDARIES, ExcMessage("invalid BC index.") );
        Assert( f->n_components == 1, ExcMessage("Function for BC should have 1 component.") );
        bc_segments[boundary_index].first= bt;
        bc_segments[boundary_index].second= f;
    }

    double fq(double h) {return fq(h,h);}

    double fq(const double h, double &dfdx) {
      B<double> x(h), f(fq_diff(x));
      f.diff(0,1);
      dfdx=x.d(0);
      return f.val();
    }

    double fk(double h) {return fk(h,h);}

    double fk(const double h, double &dfdx) {
      B<double> x(h), f(fk_diff(x));
      f.diff(0,1);
      dfdx=x.d(0);
      return f.val();
    }

    void set_time(double t) {
        k_inverse.set_time(t);
        initial->set_time(t);
        anal_sol->set_time(t);
        anal_flux->set_time(t);
        for(unsigned int i=0; i<MAX_NUM_OF_DEALII_BOUNDARIES; i++)
            if (bc_segments[i].second != NULL) bc_segments[i].second->set_time(t);

    }

    inline double pressure(double p_head, Point<dim> point)
    {
        return p_head - point[dim-1]*density*gravity;
    }

    BCType bc_type(unsigned int idx)
    {
        Assert( idx < MAX_NUM_OF_DEALII_BOUNDARIES, ExcMessage("invalid BC index.") );
        return bc_segments[idx].first;
    }

    Function<dim> *bc_func(unsigned int idx)
    {
        Assert( idx < MAX_NUM_OF_DEALII_BOUNDARIES, ExcMessage("invalid BC index.") );
        Assert( bc_segments[idx].second != NULL, ExcMessage("BC function not initialized.") );
        return bc_segments[idx].second;
    }

    void print_mat_table() {
        cout << "MATERIAL TABLE" <<
        cout << "(h, sat, cond)" << endl;
        for(double h = -50; h < 0.1; h+=0.1) {
            cout << h << " " << fq(h) << " " << fk(h) << endl;
        }
    }



};





#endif	/* _SPATIAL_FUNCTIONS_HH */

