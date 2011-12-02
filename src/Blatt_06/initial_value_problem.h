#ifndef INITIAL_VALUE_PROBLEM_H_
#define INITIAL_VALUE_PROBLEM_H_

/** @brief Example class for a differential equation model

    The model is

    u'(t) = A*u(t)

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<template <typename> class MATRIX, template <typename> class VECTOR, class N, class T = N>
class InitialValueProblem {
private:
  const MATRIX<N> A;
  const VECTOR<N> u0;
  const T t0;

public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  InitialValueProblem (const MATRIX<N>& A_, const VECTOR<N>& u0_, const T& t0_)
    : A(A_), u0(u0_), t0(t0_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return u0.size();
  }

  //! set initial state including time value
  void initialize (T& t0, VECTOR<N>& x0) const
  {
    t0 = this->t0;
    x0[0] = u0;
  }

  //! model evaluation
  void f (const T& t, const VECTOR<N>& x, VECTOR<N>& result) const
  {
    A.mv(result,x);
  }

  //! jacobian evaluation needed for implicit solvers
  void f_x (const T& t, const VECTOR<N>& x, MATRIX<N>& result) const
  {
    result = A;
  }
};



#endif /* INITIAL_VALUE_PROBLEM_H_ */
