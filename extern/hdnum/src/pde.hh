// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_PDE_HH
#define HDNUM_PDE_HH

#include<vector>
#include "newton.hh"

/** @file
 *  @brief solvers for partial differential equations
 */

namespace hdnum {

  /** @brief Stationary problem solver. E.g. for elliptic problmes

      The PDE solver is parametrized by a model. The model also
      exports all relevant types for the solution.
      The PDE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class StationarySolver
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    StationarySolver (const M& model_)
      : model(model_), x(model.size())
    {
    }

    //! do one step
    void solve ()
    {
      const size_t n_dofs = model.size();
  
      DenseMatrix<number_type> A(n_dofs,n_dofs,0.);
      Vector<number_type> b(n_dofs,0.);

      Vector<number_type> s(n_dofs);              // scaling factors
      Array<size_t> p(n_dofs);                 // row permutations
      Array<size_t> q(n_dofs);                 // column permutations

      number_type t = 0.;

      x = 0.;

      model.f_x(t, x, A);
      model.f(t, x, b);

      b*=-1.;

      row_equilibrate(A,s);                         // equilibrate rows
      lr_fullpivot(A,p,q);                          // LR decomposition of A
      apply_equilibrate(s,b);                       // equilibration of right hand side
      permute_forward(p,b);                         // permutation of right hand side
      solveL(A,b,b);                                // forward substitution
      solveR(A,x,b);                                // backward substitution
      permute_backward(q,x);                        // backward permutation
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return x;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 2;
    }
    
  private:
    const M& model;
    Vector<number_type> x;
  };


  //! gnuplot output for stationary state
  template<class N, class G>
  inline void pde_gnuplot2d (const std::string& fname, const Vector<N> solution, 
                             const G & grid)
  {

    const std::vector<Vector<N> > coords = grid.getNodeCoordinates();
    Vector<typename G::size_type> gsize = grid.getGridSize();

    std::fstream f(fname.c_str(),std::ios::out);
    // f << "set dgrid3d ";

    // f << gsize[0] << "," << gsize[1] << std::endl;

    // f << "set hidden3d" << std::endl;
    f << "set ticslevel 0" << std::endl;
    f << "splot \"-\" using 1:2:3 with points" << std::endl;
    f << "#" << std::endl;
    for (typename Vector<N>::size_type n=0; n<solution.size(); n++)
      {
        for (typename Vector<N>::size_type d=0; d<coords[n].size(); d++){
          f << std::scientific << std::showpoint 
            << std::setprecision(16) << coords[n][d] << " ";
        }

        f << std::scientific << std::showpoint 
          << std::setprecision(solution.precision()) << solution[n];

        f << std::endl;
      }
    f << "end" << std::endl;
    f << "pause -1" << std::endl;
    f.close();
  }


}
#endif
