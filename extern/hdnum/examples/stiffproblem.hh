template<class T, class N=T>
class StiffProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  StiffProblem () : ctr(0)
  {
  }

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 2;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = 0;
    x0[0] = 1.0;
    x0[1] = 2.0;
  }

  //! model evaluation
  void f (const T& t, const Vector<N>& x, Vector<N>& result) const
  {
    result[0] = 998.0*x[0] + 1998.0*x[1];
    result[1] = -999.0*x[0] + -1999.0*x[1];
    ctr++;
  }

  //! model evaluation
  void f_x (const T& t, const Vector<N>& x, Matrix<N>& result) const
  {
    result[0][0] = 998.0;      result[0][1] = 1998.0;
    result[1][0] = -999.0;  result[1][1] = -1999.0;
  }

  size_type get_count () const
  {
    return ctr;
  }

private:
  mutable size_type ctr;
};

