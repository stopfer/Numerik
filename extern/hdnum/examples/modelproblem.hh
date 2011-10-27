/** @brief Example class for a differential equation model 

    The model is

    u'(t) = lambda*u(t), t>=t_0, u(t_0) = u_0.

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class ModelProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  ModelProblem (const N& lambda_)
    : lambda(lambda_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 1;
  }

  //! set initial state including time value
  void initialize (T& t0, hdnum::Vector<N>& x0) const
  {
    t0 = 0;
    x0[0] = 1.0;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = lambda*x[0];
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
    result[0] = lambda;
  }

private:
  N lambda;
};
