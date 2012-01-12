#ifndef RENTIERE_MODEL_HH_
#define RENTIERE_MODEL_HH_

/** @brief differential equation model of reindeer population

    The model is

    u'(t) = 0.015u(t) − 3 ∗ 10^(−5) u(t)2 + 0.3, u(t_0) = u_0.

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class ReindeerModel
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  ReindeerModel (const N& u_0_, const N& t_0_)
    : u_0(u_0_), t_0(t_0_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 1;
  }

  //! set initial state including time value
  void initialize (T& t0, hdnum::Vector<N>& x0) const
  {
    t0 = t_0;//0;//t_0;
    x0[0] = u_0;//100;//u_0;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = 0.015 * x[0] - 0.00003 * x[0] * x[0] + 0.3;
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
    throw("Jacobian evaluation not implemented");
  }

private:
  N u_0;
  N t_0;
};


#endif /* RENTIERE_MODEL_HH_ */
