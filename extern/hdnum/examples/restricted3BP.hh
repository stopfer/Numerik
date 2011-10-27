#include "nbody.hh"

/** @brief restricted 3 body example

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class Restricted3Body : public NBody<T,N,2>
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  // make a the two body model 
  Restricted3Body () : NBody<T,N,2>(3)
  {
    this->m[0] = 1; 
    this->m[1] = 1; 
    this->m[2] = 0.01;
    this->G = 1.0;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = 0;

    number_type v(0.5);

    x0[0] = 1;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = v;

    x0[4] = -1;
    x0[5] = 0;
    x0[6] = 0;
    x0[7] = -v;

    x0[8] = -0.08;
    //x0[8] = 0;  
    x0[9] = 0;
    x0[10] = -0.3575;
    //x0[10] = -0.3540;
    x0[11] = 0;

    normalize(x0);
  }

};
