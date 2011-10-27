/**
   This example illustrates the 2D elliptic problem corresponding to a
   circular domain with an inward edge and peridodic Dirichlet
   boundary conditions.
 */

#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include <gmpxx.h>
#include "hdnum.hh"

#include "examples/laplace.hh"

//! Pi constant
const double pi = 3.141592653;

//! Helper function to calculate the angle of vector x with the
//! positive x-axes.
template<class N>
N getPhi(const Vector<N> & x){
  N phi = 0;
  if(x[0] < 0 && x[1] >= 0){
    phi += pi / 2.;
    phi += abs(std::atan(x[0]/x[1]));
  }
  else if(x[0] < 0 && x[1] < 0){
    phi += pi;
    phi += abs(std::atan(x[1]/x[0]));
  }
  else if(x[0] >= 0 && x[1] < 0){
    phi += pi * 1.5;
    phi += abs(std::atan(x[0]/x[1]));
  }
  else{
    phi += abs(std::atan(x[1]/x[0]));
  }
  return phi;
}

//! The domain function is passed to the model class. It defines the
//! computational domain Omega (it cuts Omega out of the positive
//! quadrant defined by the vector "extent", which is defined further
//! below).
template<class N>
class InwardEdgeDomainFunction
{
public:
  const N Phi;

  InwardEdgeDomainFunction(N Phi_)
    : Phi(Phi_) {}

  typedef N number_type;
  bool evaluate(Vector<number_type> & c) const
  {
    Vector<number_type> x = copy(c);
    x[0] -= 1.;
    x[1] -= 1.;

    if(x[0]*x[0] + x[1]*x[1] > 1)
      return false;

    N phi = getPhi(x);

    if(phi < Phi)
      return true;
    else
      return false;
  }
};

//! The boundary functor determines for each grid boundary node,
//! whether it is a Dirichlet boundary node (otherwise it is assumed
//! to be a Neumann boundary node). It also provides the Dirichlet or
//! Neumann values.
template<class N>
class InwardEdgeBoundaryFunction
{
public:

  const N Phi;

  InwardEdgeBoundaryFunction(N Phi_)
    : Phi(Phi_) {}

  typedef N number_type;
  bool isDirichlet(const Vector<number_type> & ) const
  {
    return true;
  }

  number_type getDirichletValue(const number_type , 
                                const Vector<number_type> & c) const
  {
    Vector<number_type> x = copy(c);
    x[0] -= 1.;
    x[1] -= 1.;

    const N phi = getPhi(x);

    const N av = norm(x);
    
    if(av < 0.8)
      return 0;
    
    N v = std::sin(pi /Phi * phi);
    return v;
  }

  Vector<number_type> getNeumannValue(const number_type , 
                                      const Vector<number_type> & c) const
  {
    Vector<number_type> sl(c.size());
    sl = 1;
    return copy(sl);
  }
};

using namespace hdnum; 
int main ()
{
  // The number type.
  typedef double Number; 

  // The spatial dimension.
  const int dim = 2;

  // The extent of the grid (meaning the grid domain spans the
  // quadrant [0,0] X [extent[0],extent[1]]. Notice, that the actual
  // computational domain may be only a subset of this space,
  // depending on the domain function.
  Vector<Number> extent(dim);
  extent[0]=2.; extent[1]=2.;

  // The grid resolution. Here it is 31x31 nodes i.e. 30x30 cells.
  Vector<size_t> size(dim);
  size[0]=31; size[1]=31;

  std::cout << "Inward Edge Problem" << std::endl;
  {
    typedef InwardEdgeBoundaryFunction<Number> BoundaryFunction;
    typedef InwardEdgeDomainFunction<Number> DomainFunction;
    Number Phi = pi * 1.25;
    BoundaryFunction bf(Phi);
    DomainFunction df(Phi);
    typedef LaplaceCentralDifferences<Number,Number,DomainFunction,BoundaryFunction,dim>
      Model;

    std::cout << "Setting up the model" << std::endl;
    Model model(extent,size,df,bf);

    typedef StationarySolver<Model> Solver;
    std::cout << "Setting up the solver" << std::endl;
    Solver solver(model);

    std::cout << "Solving" << std::endl;
    solver.solve();
    Vector<Number> x = solver.get_state();

    std::cout << "Output" << std::endl;
    pde_gnuplot2d("laplace.gp",x,model.getGrid()); // output model result
  }

  return 0;
}
