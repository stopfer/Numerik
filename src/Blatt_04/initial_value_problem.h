#ifndef INITIAL_VALUE_PROBLEM_H_
#define INITIAL_VALUE_PROBLEM_H_

/** @brief Example class for a differential equation model

    The model is

    u'(t) = lambda*t*u^(-1)(t)

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class InitialValueProblem {
private:
  const N lambda;
  const N u0;
  const N t0;

public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  InitialValueProblem (const N& lambda_, const N& u0_, const N& t0_)
    : lambda(lambda_), u0(u0_), t0(t0_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 1;
  }

  //! set initial state including time value
  void initialize (T& t0, hdnum::Vector<N>& x0) const
  {
    t0 = this->t0;
    x0[0] = u0;
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    result[0] = (lambda*t)/x[0];
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const hdnum::Vector<N>& x, hdnum::DenseMatrix<N>& result) const
  {
    throw std::string("Jacobian evaluation not implemented!");
  }
};



#endif /* INITIAL_VALUE_PROBLEM_H_ */
