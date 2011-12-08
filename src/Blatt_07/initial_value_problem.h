#ifndef INITIAL_VALUE_PROBLEM_H_
#define INITIAL_VALUE_PROBLEM_H_

/** @brief Example class for a differential equation model

    The model is

    u'(t) = A*u(t)

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class N, class T = N>
class InitialValueProblem {
private:
  const hdnum::Vector<N> u0;
  const T t0;

public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  InitialValueProblem ( const hdnum::Vector<N>& u0_, const T& t0_)
    : u0(u0_), t0(t0_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return u0.size();
  }

  //! set initial state including time value
  void initialize (T& t0, hdnum::Vector<N>& x0) const
  {
    t0 = this->t0;
    x0[0] = u0[0];
    x0[1] = u0[1];
  }

  //! model evaluation
  void f (const T& t, const hdnum::Vector<N>& x, hdnum::Vector<N>& result) const
  {
    // u1' = sin(u1)*sin(u2)
    // u2' = sin(u1) *sin(u2)
   result.resize(2);
    result[0]=std::sin(x[0])*sin(x[1]);
    result[1]=std::sin(x[0])*sin(x[1]);
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const hdnum::Vector<N>& x,hdnum::DenseMatrix<N>& result) const
  {
     result= hdnum::DenseMatrix<N>(2,2);
     result(0,0)=std::cos(x[0])*sin(x[1]);
     result(0,1)=std::sin(x[0])*cos(x[1]);
     result(1,0)=std::cos(x[0])*sin(x[1]);
     result(1,1)=std::sin(x[0])*cos(x[1]);
  }
};



#endif /* INITIAL_VALUE_PROBLEM_H_ */
