// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_LR_HH
#define HDNUM_LR_HH

#include "vector.hh"
#include "densematrix.hh"

/** @file
 *  @brief This file implements a generic and dynamic vector class
 */

namespace hdnum {

  //! compute lr decomposition of A with first nonzero pivoting
  template<class T>
  void lr (DenseMatrix<T>& A, Array<std::size_t>& p)
  {
    if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
      HDNUM_ERROR("need square and nonempty matrix");
    if (A.rowsize()!=p.size())
      HDNUM_ERROR("permutation vector incompatible with matrix");
	
    // transformation to upper triangular
    for (std::size_t k=0; k<A.rowsize()-1; ++k) 
      {
        // find pivot element and exchange rows
        for (std::size_t r=k; r<A.rowsize(); ++r)
          if (A[r][k]!=0) 
            {
              p[k] = r; // store permutation in step k
              if (r>k) // exchange complete row if r!=k
                for (std::size_t j=0; j<A.colsize(); ++j)
                  {
                    T temp(A[k][j]);
                    A[k][j] = A[r][j];
                    A[r][j] = temp;
                  }
              break;
            }
        if (A[k][k]==0) HDNUM_ERROR("matrix is singular");

        // modification
        for (std::size_t i=k+1; i<A.rowsize(); ++i) 
          {
            T qik(A[i][k]/A[k][k]);
            A[i][k] = qik;
            for (std::size_t j=k+1; j<A.colsize(); ++j)
              A[i][j] -= qik * A[k][j];
          }
      }
  }

  //! our own abs class that works also for multiprecision types
  template<class T>
  T abs (const T& t)
  {
    if (t<0.0)
      return -t;
    else
      return t;
  }

  //! lr decomposition of A with column pivoting
  template<class T>
  void lr_partialpivot (DenseMatrix<T>& A, Array<std::size_t>& p)
  {
    if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
      HDNUM_ERROR("need square and nonempty matrix");
    if (A.rowsize()!=p.size())
      HDNUM_ERROR("permutation vector incompatible with matrix");
	
    // initialize permutation
    for (std::size_t k=0; k<A.rowsize(); ++k)
      p[k] = k;

    // transformation to upper triangular
    for (std::size_t k=0; k<A.rowsize()-1; ++k) 
      {
        // find pivot element
        for (std::size_t r=k+1; r<A.rowsize(); ++r)
          if (abs(A[r][k])>abs(A[k][k])) 
            p[k] = r; // store permutation in step k

        if (p[k]>k) // exchange complete row if r!=k
          for (std::size_t j=0; j<A.colsize(); ++j)
            {
              T temp(A[k][j]);
              A[k][j] = A[p[k]][j];
              A[p[k]][j] = temp;
            }

        if (A[k][k]==0) HDNUM_ERROR("matrix is singular");

        // modification
        for (std::size_t i=k+1; i<A.rowsize(); ++i) 
          {
            T qik(A[i][k]/A[k][k]);
            A[i][k] = qik;
            for (std::size_t j=k+1; j<A.colsize(); ++j)
              A[i][j] -= qik * A[k][j];
          }
      }
  }

  //! lr decomposition of A with full pivoting
  template<class T>
  void lr_fullpivot (DenseMatrix<T>& A, Array<std::size_t>& p, Array<std::size_t>& q)
  {
    if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
      HDNUM_ERROR("need square and nonempty matrix");
    if (A.rowsize()!=p.size())
      HDNUM_ERROR("permutation vector incompatible with matrix");

    // initialize permutation
    for (std::size_t k=0; k<A.rowsize(); ++k)
      p[k] = q[k] = k;

    // transformation to upper triangular
    for (std::size_t k=0; k<A.rowsize()-1; ++k) 
      {
        // find pivot element
        for (std::size_t r=k; r<A.rowsize(); ++r)
          for (std::size_t s=k; s<A.colsize(); ++s)
            if (abs(A[r][s])>abs(A[k][k]))
              { 
                p[k] = r; // store permutation in step k
                q[k] = s;
              }

        if (p[k]>k) // exchange complete row if r!=k
          for (std::size_t j=0; j<A.colsize(); ++j)
            {
              T temp(A[k][j]);
              A[k][j] = A[p[k]][j];
              A[p[k]][j] = temp;
            }
        if (q[k]>k) // exchange complete column if s!=k
          for (std::size_t i=0; i<A.rowsize(); ++i)
            {
              T temp(A[i][k]);
              A[i][k] = A[i][q[k]];
              A[i][q[k]] = temp;
            }

        if (A[k][k]==0) HDNUM_ERROR("matrix is singular");

        // modification
        for (std::size_t i=k+1; i<A.rowsize(); ++i) 
          {
            T qik(A[i][k]/A[k][k]);
            A[i][k] = qik;
            for (std::size_t j=k+1; j<A.colsize(); ++j)
              A[i][j] -= qik * A[k][j];
          }
      }
  }

  //! apply permutations to a right hand side vector
  template<class T>
  void permute_forward (const Array<std::size_t>& p, Vector<T>& b)
  {
    if (b.size()!=p.size())
      HDNUM_ERROR("permutation vector incompatible with rhs");
  
    for (std::size_t k=0; k<b.size()-1; ++k)
      if (p[k]!=k)
        {
          T temp(b[k]);
          b[k] = b[p[k]];
          b[p[k]] = temp;
        }
  }

  //! apply permutations to a solution vector
  template<class T>
  void permute_backward (const Array<std::size_t>& q, Vector<T>& z)
  {
    if (z.size()!=q.size())
      HDNUM_ERROR("permutation vector incompatible with z");
  
    for (int k=z.size()-2; k>=0; --k)
      if (q[k]!=std::size_t(k))
        {
          T temp(z[k]);
          z[k] = z[q[k]];
          z[q[k]] = temp;
        }
  }

  //! perform a row equilibration of a matrix; return scaling for later use
  template<class T>
  void row_equilibrate (DenseMatrix<T>& A, Vector<T>& s)
  {
    if (A.rowsize()*A.colsize()==0) 
      HDNUM_ERROR("need nonempty matrix");
    if (A.rowsize()!=s.size())
      HDNUM_ERROR("scaling vector incompatible with matrix");
  
    // equilibrate row sums
    for (std::size_t k=0; k<A.rowsize(); ++k) 
      {
        s[k] = T(0.0);
        for (std::size_t j=0; j<A.colsize(); ++j)
          s[k] += abs(A[k][j]);
        if (s[k]==0) HDNUM_ERROR("row sum is zero");
        for (std::size_t j=0; j<A.colsize(); ++j)
          A[k][j] /= s[k];
      }
  }

  //! apply row equilibration to right hand side vector
  template<class T>
  void apply_equilibrate (Vector<T>& s, Vector<T>& b)
  {
    if (s.size()!=b.size())
      HDNUM_ERROR("s and b incompatible");
  
    // equilibrate row sums
    for (std::size_t k=0; k<b.size(); ++k) 
      b[k] /= s[k];
  }

  //! Assume L = lower triangle of A with l_ii=1, solve L x = b
  template<class T>
  void solveL (const DenseMatrix<T>& A, Vector<T>& x, const Vector<T>& b)
  {
    if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
      HDNUM_ERROR("need square and nonempty matrix");
    if (A.rowsize()!=b.size())
      HDNUM_ERROR("right hand side incompatible with matrix");

    for (std::size_t i=0; i<A.rowsize(); ++i)
      {
        T rhs(b[i]);
        for (std::size_t j=0; j<i; j++)
          rhs -= A[i][j] * x[j];
        x[i] = rhs;
      }
  }

  //! Assume R = upper triangle of A and solve R x = b
  template<class T>
  void solveR (const DenseMatrix<T>& A, Vector<T>& x, const Vector<T>& b)
  {
    if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
      HDNUM_ERROR("need square and nonempty matrix");
    if (A.rowsize()!=b.size())
      HDNUM_ERROR("right hand side incompatible with matrix");

    for (int i=A.rowsize()-1; i>=0; --i)
      {
        T rhs(b[i]);
        for (std::size_t j=i+1; j<A.colsize(); j++)
          rhs -= A[i][j] * x[j];
        x[i] = rhs/A[i][i];
      }
  }

}
#endif

