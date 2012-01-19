/** @brief Example class for a differential equation model 

    The model is

    u'(t) = - lambda * t * u(t)^2, t>=t_0, u(t_0) = u_0.

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class ConvergenceModel
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  ConvergenceModel (const N& lambda_, const N& u_0_, const N& t_0_)
    : lambda(lambda_), u_0(u_0_), t_0(t_0_), counter(0)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 1;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = t_0;
    x0[0] = u_0;
  }

  //! model evaluation
  void f (const T& t, const Vector<N>& x, Vector<N>& result) const
  {
    result[0] = - t * lambda * x[0] * x[0];
    counter=counter+1;
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const Vector<N>& x, DenseMatrix<N>& result) const
  {
    N f_value; f(t,x,f_value);
    result[0] = - ( lambda * x[0] * x[0] + t * 2 * lambda * x[0] * f_value);
  }

private:
  N lambda;
  N u_0;
  N t_0;
public:
  mutable int counter;
};
