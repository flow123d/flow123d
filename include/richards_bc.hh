/* 
 * File:   richards_bc.hh
 * Author: jb
 *
 * Created on June 10, 2010, 1:25 PM
 */

#ifndef _RICHARDS_BC_HH
#define	_RICHARDS_BC_HH

#include <base/subscriptor.h>

/**
 * Boundary condition class provides type of the boundary condition and
 * its value. For further use, there can be also supprot for Newton (returns a vector)
 * Seepage, and time dependent, to this end ther should be some notification capability,
 * to update the condition values (e.g according to time)
 *
 * The assumption about the cocept is that every BC can be resolved as Diriclet, Neuman or Newton
 * this may not be the case fi we would like advanced resolve for seepage or atmo.
 */

template <int dim>
class BoundaryCondition : public Function<dim> {
public:
    typedef enum { Dirichlet, Neuman } BCType;
    BoundaryCondition(BCType type, double value)
      : Function<dim>(1),bc_type(type), bc_value(value) {}

    BCType type() const
        {return bc_type;}

    virtual double value (const Point<dim>   &p, const unsigned int component = 0) const
        { return bc_value; }

private:
    BCType bc_type;
    double bc_value;
};

typedef BoundaryCondition<2> BC2D;
typedef BoundaryCondition<3> BC3D;


#endif	/* _RICHARDS_BC_HH */

