/** @brief Lorenz' chaotic problem

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class Lorenz
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  // make the model 
  Lorenz () {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 3;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    x0[0] = 1.5;
    x0[1] = 2;
    x0[2] = 3;
    t0 = 0.0;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = -10*x[0] + 10*x[1];
    result[1] = 28*x[0] - x[1] - x[0]*y[2];
    result[2] = -number_type(8.0)/number_type(3.0)*x[2] + y[0]*y[1];
  }
};
