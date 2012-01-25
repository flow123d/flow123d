/* 
 * File:   spatial_functions.hh
 * Author: jb
 *
 * Created on May 17, 2010, 1:50 PM
 */

#ifndef _SPATIAL_FUNCTIONS_HH
#define	_SPATIAL_FUNCTIONS_HH

#include <iostream>
#include <math.h>

#include <base/function.h>
#include <base/tensor_function.h>
#include <base/tensor_base.h>
#include <base/tensor.h>

#include <hydro_functions.hh>
//#include <richards_bc.hh>
#include "FADBAD++/fadbad.h"
#include "FADBAD++/badiff.h"
#include "FADBAD++/fadiff.h"
#include "base/parameter_handler.h"
using namespace fadbad;
using namespace dealii;
using namespace std;

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


template <int dim>
class ZLinear : public Function<dim>
{
  public:
    ZLinear (ParameterHandler &prm, double bottom_val, double top_val) : Function<dim>(1) {

        double x0 = - prm.get_double("z_size");
        double x1 =0;
        a1 = (bottom_val - top_val)/(x0-x1);
        a0 = bottom_val - a1*x0;
    }
    virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const
    {
        return (a1+1)*p[dim-1]+a0 ;
    }
    virtual ~ZLinear() {}

  private:
    double a0, a1;
};


/**
 *  BOUNDARY FUNCTIONS
 */

template <int dim>
class TLinear : public Function<dim>
{
  public:
    TLinear (ParameterHandler &prm) : Function<dim>(1) {

        double t_limit = 100;
        a0 = -2.0;
        limit= 0.2;
        a1 = (limit-a0)/t_limit;
    }
    virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const
    {
        //std::cout << "bc t " << this->get_time() << std::endl;
        return min(limit,  this->get_time() * a1 + a0) + p[dim-1];
    }
    virtual ~TLinear() {}

  private:
    double a0, a1, limit;
};


/**
 * Sliding analytical solution with ATan saturation. (pressure)
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

/**
 * Sliding analytical solution with ATan saturation. (flux = -K * grad piezo_head), however piezo_head == p_head, zero gravity
 */
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

    double cap_max_;
    double cap_arg_max_;
    double lambda_cap_max_half;
private:
    //! there should be whole BC descriptor object which returns
    //! BC objects for given index with checking, possibly returning None type of BC
    std::vector< std::pair<BCType, Function<dim> *> > bc_segments;

    HydrologyModel hydro_model;

    double density;
    double gravity;


    bool no_exact_solution;

public:

    RichardsData( ParameterHandler &prm):
    k_inverse(prm),
    initial(NULL),
    bc_segments(MAX_NUM_OF_DEALII_BOUNDARIES, std::pair<BCType, Function<dim> * >(None, (Function<dim> *)(NULL) ) ),
    hydro_model(prm),
    //fq_diff(h_params),
    //fc(h_params),
    //fk_diff(h_params),
    density(1.0),
    gravity(1.0)
    {
        cap_arg_max_=hydro_model.cap_arg_max();
        fq(cap_arg_max_/1.5, lambda_cap_max_half);
        lambda_cap_max_half/=fk(cap_arg_max_/1.5);
        cout << "lh "<<lambda_cap_max_half <<endl;
        fq(cap_arg_max_, cap_max_);


        add_bc(0, Neuman, new ZeroFunction<dim>() ); // left
        add_bc(1, Neuman, new ZeroFunction<dim>() ); // right
        add_bc(3, Dirichlet, new TLinear<dim>(prm) ); // top
        add_bc(2, Dirichlet, new ConstantFunction<dim>(-2) ); // bottom

        initial = new ZLinear<dim>(prm, -1, -2);
        anal_sol = new ZeroFunction<dim>();
        anal_flux = new ZeroFunction<dim>();
        no_exact_solution =true;

        /*
        // Analytical solution setting

        add_bc(0, Neuman, new ZeroFunction<dim>() ); // left
        add_bc(1, Neuman, new ZeroFunction<dim>() ); // right
        add_bc(2, Dirichlet, new ASol_ATan<dim>() ); // bottom
        add_bc(3, Dirichlet, new ASol_ATan<dim>()); // top
        //add_bc(3,Dirichlet, new ConstantFunction<dim>(1.0); // top
        //add_bc(3,Dirichlet, new ConstantFunction<dim>(1.0); // top

        initial = new ASol_ATan<dim>();
        anal_sol = new ASol_ATan<dim>();
        anal_flux = new AFlux_ATan<dim>();
*/
        //anal_sol = new ASol_lin<dim>();
        //initial = new ASol_lin<dim>();

        print_mat_table();
    }

    ~RichardsData() {
        for(unsigned int i=0; i<MAX_NUM_OF_DEALII_BOUNDARIES;i++)
            if (bc_segments[i].second != NULL) delete bc_segments[i].second;
    }

    //inline double lambda(const double cap) const {return hydro_model.lambda(cap / cap_max_);}

    bool has_exact_solution() { return ! no_exact_solution; }

    void add_bc(const unsigned int boundary_index, const BCType bt, Function<dim> * f)
    {
        Assert( boundary_index < MAX_NUM_OF_DEALII_BOUNDARIES, ExcMessage("invalid BC index.") );
        Assert( f->n_components == 1, ExcMessage("Function for BC should have 1 component.") );
        bc_segments[boundary_index].first= bt;
        bc_segments[boundary_index].second= f;
    }

    double fq(double h) {return fq(h,h);}

    double fq(const double h, double &dfdx) {
      B<double> x(h);
      B<double> f( hydro_model.FQ(x) );
      f.diff(0,1);
      dfdx=x.d(0);
      return f.val();
    }

    // move this into hydro model ?
    double fqq(const double h) {
        B< F<double> > x(h);
        x.x().diff(0,2);
        B< F<double> > f( hydro_model.FQ(x) );
        f.diff(0,1);
        return x.d(0).d(0);
    }

    double fk(double h) {return fk(h,h);}

    double fk(const double h, double &dfdx) {
      B<double> x(h);
      B<double> f( hydro_model.FK(x) );
      f.diff(0,1);
      dfdx=x.d(0);
      return f.val();
    }

    void set_time(double t, double dt) {
        k_inverse.set_time(t);
        initial->set_time(t);
        anal_sol->set_time(t);
        anal_flux->set_time(t-0.5*dt);
        //anal_flux->set_time(t);
        for(unsigned int i=0; i<MAX_NUM_OF_DEALII_BOUNDARIES; i++)
            if (bc_segments[i].second != NULL) bc_segments[i].second->set_time(t);

    }

    inline double pressure(const double p_head,const  Point<dim> &point)
    {
        return p_head - point[dim-1]*density*gravity;
    }

    inline double pressure(const double p_head,const double z_coord)
    {
        return p_head - z_coord*density*gravity;
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
        std::cout << "MATERIAL TABLE" <<
        std::cout << "(h, sat, cap, con, con_diff, cap_diff)" << std::endl;
        double cap=0,k_diff=0;
        for(double e = -5; e < 2; e+=0.1) {
            double h=-5*exp(e);
            std::cout << h << " " << fq(h,cap) << " " << cap << " "<< fk(h,k_diff) << " "<<k_diff <<" "<< fqq(h) <<std::endl;
        }
    }



};





#endif	/* _SPATIAL_FUNCTIONS_HH */

