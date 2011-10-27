// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_EXAMPLES_NBODY_HH
#define HDNUM_EXAMPLES_NBODY_HH

#include<math.h>

/** @brief generic astrophysical n-body model

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
    \tparam n number of bodies
    \tparam d space dimension

    Uses the following data layout:

    (x0_0,...,x0_{d-1},v0_0,...,v0_{d-1}, x01_0,...,x1_{d-1},v1_0,...,v1_{d-1}, ...)

    The page http://ssd.jpl.nasa.gov/?planets#elem might be interesting
    for someone who wants to build in a complete solar system.
*/
template<class T, class N, int d>
class NBody
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 2*n*d;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    for (size_type i=0; i<n; i++)
      {
        size_type I=2*d*i;
        for (size_type k=0; k<d; k++) 
          {
            result[I+k] = x[I+d+k];
            result[I+d+k] = 0.0;
          }
        for (size_type j=0; j<i; j++)
          {
            size_type J=2*d*j;
            number_type r2(0.0);
            for (size_type k=0; k<d; k++) r2 += (x[J+k]-x[I+k])*(x[J+k]-x[I+k]);
            number_type r3(r2*sqrt(r2));
            for (size_type k=0; k<d; k++)
              {
                result[I+d+k] += m[j]*(x[J+k]-x[I+k])/r3;
                result[J+d+k] += m[i]*(x[I+k]-x[J+k])/r3;
              }
          }
      }
    for (size_type i=0; i<n; i++)
      for (size_type k=0; k<d; k++) result[2*d*i+d+k] *= G;
    ctr++;
  }

  // set center of mass and its velocity to zero
  void normalize (hdnum::Vector<N>& x) const
  {
    number_type s[d],t[d];
    for (size_type k=0; k<d; k++) s[k]=t[k]=0.0;
    
    for (size_type i=0; i<n; i++)
      for (size_type k=0; k<d; k++) 
        {
          s[k] += m[i]*x[2*d*i+k];   // center of mass
          t[k] += m[i]*x[2*d*i+d+k]; // first moment of velocity
        }
    for (size_type i=0; i<n; i++)
      for (size_type k=0; k<d; k++) 
        {
          x[2*d*i+k] -= s[k];
          x[2*d*i+d+k] -= t[k];
        }
  }

  // set center of mass and its velocity to zero
  number_type energy (const hdnum::Vector<N>& x) const
  {
    number_type e(0.0);
    for (size_type i=0; i<n; i++)
      {
        size_type I=2*d*i;
        number_type ei(0.0);
        for (size_type k=0; k<d; k++) 
          ei += 0.5*m[i]*x[I+d+k]*x[I+d+k];
        e += ei;
        ei = 0.0;
        for (size_type j=0; j<i; j++)
          {
            size_type J=2*d*j;
            number_type r2(0.0);
            for (size_type k=0; k<d; k++) r2 += (x[J+k]-x[I+k])*(x[J+k]-x[I+k]);
            number_type r(sqrt(r2));
            ei += G*m[i]*m[j]/r;
          }
        e -= ei;
      }
    return e;
  }

  size_type get_count () const
  {
    return ctr;
  }

protected:
  // may only be called by derived class
  NBody (size_type n_) : n(n_), ctr(0), m(2*n_*d)
  {}

  size_type n; 
  mutable size_type ctr;
  number_type G;                // gravitational constant
  hdnum::Vector<number_type> m; // masses
};

#endif
