// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_NEWTON_HH
#define HDNUM_NEWTON_HH

#include "lr.hh"

/** @file
 *  @brief Newton's method with line search
 */

namespace hdnum {

  /** @brief Example class for a nonlinear model F(x) = 0;

      This example solves F(x) = x*x - a = 0

      \tparam N a type representing x and F components
  */
  template<class N>
  class SquareRootProblem
  {
  public:
    /** \brief export size_type */
    typedef std::size_t size_type;

    /** \brief export number_type */
    typedef N number_type;

    //! constructor stores parameter lambda
    SquareRootProblem (number_type a_)
      : a(a_)
    {}

    //! return number of componentes for the model
    std::size_t size () const
    {
      return 1;
    }

    //! model evaluation
    void F (const Vector<N>& x, Vector<N>& result) const
    {
      result[0] = x[0]*x[0] - a;
    }

    //! jacobian evaluation needed for implicit solvers
    void F_x (const Vector<N>& x, DenseMatrix<N>& result) const
    {
      result[0][0] = number_type(2.0)*x[0];
    }

  private:
    number_type a;
  };


  /** @brief Solve nonlinear problem using a damped Newton method

      The Newton solver is parametrized by a model. The model also
      exports all relevant types for types.

  */
  class Newton
  {
    typedef std::size_t size_type;

  public:
    //! constructor stores reference to the model
    Newton ()
      : maxit(25), linesearchsteps(10), verbosity(0), 
	reduction(1e-14), abslimit(1e-30), converged(false)
    {}

    //! maximum number of iterations before giving up
    void set_maxit (size_type n)
    {
      maxit = n;
    }

    //! maximum number of steps in linesearch before giving up
    void set_linesearchsteps (size_type n)
    {
      linesearchsteps = n;
    }

    //! control output given 0=nothing, 1=summary, 2=every step, 3=include line search
    void set_verbosity (size_type n)
    {
      verbosity = n;
    }

    //! basolute limit for defect
    void set_abslimit (double l)
    {
      abslimit = l;
    }

    //! reduction factor
    void set_reduction (double l)
    {
      reduction = l;
    }

    //! do one step
    template<class M>
    void solve (const M& model, Vector<typename M::number_type> & x) const
    {
      typedef typename M::number_type N;
      Vector<N> r(model.size());              // residual
      DenseMatrix<N> A(model.size(),model.size()); // Jacobian matrix
      Vector<N> y(model.size());              // temporary solution in line search
      Vector<N> z(model.size());              // solution of linear system
      Vector<N> s(model.size());              // scaling factors
      Array<size_type> p(model.size());                 // row permutations
      Array<size_type> q(model.size());                 // column permutations

      model.F(x,r);                                     // compute nonlinear residual
      N R0(norm(r));                          // norm of initial residual
      N R(R0);                                // current residual norm
      if (verbosity>=1)
	{
	  std::cout << "Newton " 
		    << "   norm=" << std::scientific << std::showpoint 
		    << std::setprecision(4) << R0
		    << std::endl;
	}

      converged = false;
      for (size_type i=1; i<=maxit; i++)                // do Newton iterations
	{
	  // check absolute size of residual
	  if (R<=abslimit)
	    {
	      converged = true;
	      return;
	    } 

	  // solve Jacobian system for update
	  model.F_x(x,A);                               // compute Jacobian matrix
	  row_equilibrate(A,s);                         // equilibrate rows
	  lr_fullpivot(A,p,q);                          // LR decomposition of A
	  z = N(0.0);                                   // clear solution
	  apply_equilibrate(s,r);                       // equilibration of right hand side
	  permute_forward(p,r);                         // permutation of right hand side
	  solveL(A,r,r);                                // forward substitution
	  solveR(A,z,r);                                // backward substitution
	  permute_backward(q,z);                        // backward permutation

	  // line search
	  N lambda(1.0);                      // start with lambda=1
	  for (size_type k=0; k<linesearchsteps; k++)
	    {
	      y = x;                                    
	      y.update(-lambda,z);                       // y = x+lambda*z
	      model.F(y,r);                             // r = F(y)
	      N newR(norm(r));                // compute norm
	      if (verbosity>=3)
		{
		  std::cout << "    line search "  << std::setw(2) << k 
			    << " lambda=" << std::scientific << std::showpoint 
			    << std::setprecision(4) << lambda
			    << " norm=" << std::scientific << std::showpoint 
			    << std::setprecision(4) << newR
			    << " red=" << std::scientific << std::showpoint 
			    << std::setprecision(4) << newR/R
			    << std::endl;
		}
	      if (newR<(1.0-0.25*lambda)*R)            // check convergence
		{
		  if (verbosity>=2)
		    {
		      std::cout << "  step"  << std::setw(3) << i 
				<< " norm=" << std::scientific << std::showpoint 
				<< std::setprecision(4) << newR
				<< " red=" << std::scientific << std::showpoint 
				<< std::setprecision(4) << newR/R
				<< std::endl;
		    }
		  x = y;
		  R = newR;
		  break;                                // continue with Newton loop
		}
	      else lambda *= 0.5;                       // reduce damping factor
	      if (k==linesearchsteps-1)
		{
		  if (verbosity>=3)
		    std::cout << "    line search not converged within " << linesearchsteps << " steps" << std::endl;
		  return;
		}
	    }

	  // check convergence
	  if (R<=reduction*R0)
	    {
	      if (verbosity>=1)
		{
		  std::cout << "Newton converged in "  << i << " steps"
			    << " reduction=" << std::scientific << std::showpoint 
			    << std::setprecision(4) << R/R0
			    << std::endl;
		}
	      converged = true;
	      return;
	    }
	  if (i==maxit)
	    {
	      if (verbosity>=2)
		std::cout << "Newton not converged within " << maxit << " iterations" << std::endl;
	    }
	}
    }
    
    bool has_converged () const
    {
      return converged;
    }

  private:
    size_type maxit;
    size_type linesearchsteps;
    size_type verbosity;
    double reduction;
    double abslimit;
    mutable bool converged;
  };




  /** @brief Solve nonlinear problem using a fixed point iteration

      solve F(x) = 0.

      x = x - \sigma*F(x)

  */
  class Banach
  {
    typedef std::size_t size_type;

  public:
    //! constructor stores reference to the model
    Banach ()
      : maxit(25), linesearchsteps(10), verbosity(0), 
        reduction(1e-14), abslimit(1e-30),  sigma(1.0), converged(false)
    {}

    //! maximum number of iterations before giving up
    void set_maxit (size_type n)
    {
      maxit = n;
    }

    //! damping parameter
    void set_sigma (double sigma_)
    {
      sigma = sigma_;
    }

    //! maximum number of steps in linesearch before giving up
    void set_linesearchsteps (size_type n)
    {
      linesearchsteps = n;
    }

    //! control output given 0=nothing, 1=summary, 2=every step, 3=include line search
    void set_verbosity (size_type n)
    {
      verbosity = n;
    }

    //! basolute limit for defect
    void set_abslimit (double l)
    {
      abslimit = l;
    }

    //! reduction factor
    void set_reduction (double l)
    {
      reduction = l;
    }

    //! do one step
    template<class M>
    void solve (const M& model, Vector<typename M::number_type> x) const
    {
      typedef typename M::number_type N;
      Vector<N> r(model.size());              // residual
      Vector<N> y(model.size());              // temporary solution in line search

      model.F(x,r);                           // compute nonlinear residual
      N R0(norm(r));                          // norm of initial residual
      N R(R0);                                // current residual norm
      if (verbosity>=1)
	{
	  std::cout << "Banach " 
		    << " norm=" << std::scientific << std::showpoint 
		    << std::setprecision(4) << R0
		    << std::endl;
	}

      converged = false;
      for (size_type i=1; i<=maxit; i++)                // do iterations
	{
	  // check absolute size of residual
	  if (R<=abslimit)
	    {
	      converged = true;
	      return;
	    } 

	  // next iterate
	  y = x;                                    
	  y.update(-sigma,r);                       // y = x+lambda*z
	  model.F(y,r);                             // r = F(y)
	  N newR(norm(r));                // compute norm
	  if (verbosity>=2)
	    {
	      std::cout << "    "  << std::setw(3) << i 
			<< " norm=" << std::scientific << std::showpoint 
			<< std::setprecision(4) << newR
			<< " red=" << std::scientific << std::showpoint 
			<< std::setprecision(4) << newR/R
			<< std::endl;
	    }
	  x = y;                                // accept new iterate
	  R = newR;                             // remember new norm

	  // check convergence
	  if (R<=reduction*R0 || R<=abslimit)
	    {
	      if (verbosity>=1)
		{
		  std::cout << "Banach converged in "  << i << " steps"
			    << " reduction=" << std::scientific << std::showpoint 
			    << std::setprecision(4) << R/R0
			    << std::endl;
		}
	      converged = true;
	      return;
	    }
	}
    }
    
    bool has_converged () const
    {
      return converged;
    }

  private:
    size_type maxit;
    size_type linesearchsteps;
    size_type verbosity;
    double reduction;
    double abslimit;
    double sigma;
    mutable bool converged;
  };

} // namespace hdnum

#endif
