#include "nbody.hh"

/** @brief three bodies moving on a figure eight

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class FigureEight : public NBody<T,N,2>
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  // make a the two body model 
  FigureEight () : NBody<T,N,2>(3)
  {
    this->m[0] = 1; 
    this->m[1] = 1; 
    this->m[2] = 1;
    this->G = 1.0;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = 0;

    x0[0] = 0.9700436;
    x0[1] = -0.24308753;
    x0[2] = 0.466203685;
    x0[3] = 0.43236573;

    x0[4] = -0.9700436;
    x0[5] = 0.24308753;
    x0[6] = 0.466203685;
    x0[7] = 0.43236573;

    x0[8] = 0;
    x0[9] = 0;
    x0[10] = -0.93240737;
    x0[11] = -0.86473146;

    normalize(x0);
  }

};
