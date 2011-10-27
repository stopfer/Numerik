// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/*
 * File:   vector.hh
 * Author: ngo
 *
 * Created on April 14th, 2011
 */

#ifndef _VECTOR_HH
#define	_VECTOR_HH

#include <vector>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>


namespace hdnum {

  /*! \brief Class with mathematical vector operations
   */

  template<typename REAL>
  class Vector : public std::vector<REAL>  // inherit from the STL vector
  {
  public:
    /** \brief Type used for array indices */
    typedef std::size_t size_type;

  private:
    static bool bScientific;
    static std::size_t nIndexWidth;
    static std::size_t nValueWidth;
    static std::size_t nValuePrecision;

  public:

    // default constructor, also inherited from the STL vector default constructor
    Vector() : std::vector<REAL>()
    {
    }

    // another constructor, with arguments, setting the default value for all entries of the vector of given size
    Vector( const size_t size,                 // user must specify the size
            const REAL defaultvalue_ = 0    // if not specified, the value 0 will take effect
            )
      : std::vector<REAL>( size, defaultvalue_ )
    {
    }


    // Methods:

    /*!
      \brief Assign all values of the Vector from one scalar value: x = value

      \b Example:
      \code
      hdnum::Vector<double> x(4);
      x = 1.23;
      std::cout << "x=" << x << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]  1.2340000e+00
      [ 1]  1.2340000e+00
      [ 2]  1.2340000e+00
      [ 3]  1.2340000e+00
      \endverbatim
    */
    Vector& operator=( const REAL value )
    {
      const size_t s = this->size();
      Vector & self = *this;
      for(size_t i=0; i<s; ++i)
        self[i] = value;
      return *this;
    }

    /*!
      \brief Subvector extraction

      Returns a new vector that is a subset of the components
      of the given vector.

      \param[in] i first index of the new vector
      \param[in] m size of the new vector, i.e. it has components [i,i+m-1]
    */
    Vector sub (size_type i, size_type m)
    {
      Vector v(m);
      Vector &self = *this;
      size_type k=0;
      for (size_type j=i; j<i+m; j++){
        v[k]=self[j];
        k++;
      }
      return v;
    }



    /*!
      \brief Assigning a vector from a given vector: x = y

      \b Example:
      \code
      hdnum::Vector<double> x(4);
      hdnum::Vector<double> y(4);
      x[0] = 1.23;
      x[1] = 2.31;
      x[2] = 4.54;
      x[3] = 9.98;
      std::cout << "x=" << x << std::endl;
      y = x;
      std::cout << "y=" << y << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]  1.2300000e+00
      [ 1]  2.3100000e+00
      [ 2]  4.5400000e+00
      [ 3]  9.9800000e+00

      y=
      [ 0]  1.2300000e+00
      [ 1]  2.3100000e+00
      [ 2]  4.5400000e+00
      [ 3]  9.9800000e+00
      \endverbatim
    */
#ifdef DOXYGEN
    Vector& operator=( const Vector& y )
    {
      // It is already implemented in the STL vector class itself!
    }
#endif



    // Multiplication by a scalar value: x *= value
    Vector& operator*=( const REAL value )
    {
      Vector &self = *this;
      for (size_t i = 0; i < this->size(); ++i)
        self[i] *= value;
      return *this;
    }


    // Division by a scalar value: x /= value
    Vector& operator/=( const REAL value )
    {
      Vector &self = *this;
      for (size_t i = 0; i < this->size(); ++i)
        self[i] /= value;
      return *this;
    }


    // Add another vector: x += y
    Vector& operator+=( const Vector & y )
    {
      assert( this->size() == y.size());
      Vector &self = *this;
      for (size_t i = 0; i < this->size(); ++i)
        self[i] += y[i];
      return *this;
    }


    // Subtract another vector: x -= y
    Vector& operator-=( const Vector & y )
    {
      assert( this->size() == y.size());
      Vector &self = *this;
      for (size_t i = 0; i < this->size(); ++i)
        self[i] -= y[i];
      return *this;
    }


    // Update vector by addition of a scaled vector:
    // x += alpha * y
    //
    Vector & update(const REAL alpha, const Vector & y)
    {
      assert( this->size() == y.size());
      Vector &self = *this;
      for (size_t i = 0; i < this->size(); ++i)
        self[i] += alpha * y[i];
      return *this;
    }


    /*!
      \brief Inner product with another vector

      \b Example:
      \code
      hdnum::Vector<double> x(2);
      x.scientific(false); // set fixed point display mode
      x[0] = 12.0;
      x[1] = 3.0;
      std::cout << "x=" << x << std::endl;
      hdnum::Vector<double> y(2);
      y[0] = 4.0;
      y[1] = -1.0;
      std::cout << "y=" << y << std::endl;
      double s = x*y;
      std::cout << "s = x*y = " << s << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]     12.0000000
      [ 1]      3.0000000

      y=
      [ 0]      4.0000000
      [ 1]     -1.0000000

      s = x*y = 45.0000000
      \endverbatim
    */
    REAL operator*(Vector & x) const
    {
      assert( x.size() == this->size() );   // checks if the dimensions of the two vectors are equal
      REAL sum( 0 );
      const Vector & self = *this;
      for( size_t i = 0; i < this->size(); ++i )
        sum += self[i] * x[i];
      return sum;
    }




    /*!
      \brief Adding two vectors x+y

      \b Example:
      \code
      hdnum::Vector<double> x(2);
      x.scientific(false); // set fixed point display mode
      x[0] = 12.0;
      x[1] = 3.0;
      std::cout << "x=" << x << std::endl;
      hdnum::Vector<double> y(2);
      y[0] = 4.0;
      y[1] = -1.0;
      std::cout << "y=" << y << std::endl;
      std::cout << "x+y = " << x+y << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]     12.0000000
      [ 1]      3.0000000

      y=
      [ 0]      4.0000000
      [ 1]     -1.0000000

      x+y =
      [ 0]     16.0000000
      [ 1]      2.0000000
      \endverbatim
    */
    Vector operator+(Vector & x) const
    {
      assert( x.size() == this->size() );   // checks if the dimensions of the two vectors are equal
      Vector sum( *this );
      sum += x;
      return sum;
    }



    /*!
      \brief vector subtraction x-y

      \b Example:
      \code
      hdnum::Vector<double> x(2);
      x.scientific(false); // set fixed point display mode
      x[0] = 12.0;
      x[1] = 3.0;
      std::cout << "x=" << x << std::endl;
      hdnum::Vector<double> y(2);
      y[0] = 4.0;
      y[1] = -1.0;
      std::cout << "y=" << y << std::endl;
      std::cout << "x-y = " << x-y << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]     12.0000000
      [ 1]      3.0000000

      y=
      [ 0]      4.0000000
      [ 1]     -1.0000000

      x-y =
      [ 0]      8.0000000
      [ 1]      4.0000000
      \endverbatim
    */
    Vector operator-(Vector & x) const
    {
      assert( x.size() == this->size() );   // checks if the dimensions of the two vectors are equal
      Vector sum( *this );
      sum -= x;
      return sum;
    }



    //! Square of the Euclidean norm
    REAL two_norm_2() const
    {
      REAL sum( 0 );
      const Vector & self = *this;
      for (size_t i = 0; i < (size_t) this->size(); ++i)
        sum += self[i] * self[i];
      return sum;
    }

    /*!
      \brief Euclidean norm of a vector

      \b Example:
      \code
      hdnum::Vector<double> x(3);
      x.scientific(false); // set fixed point display mode
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 1.0;
      std::cout << "x=" << x << std::endl;
      std::cout << "euclidean norm of x = " << x.two_norm() << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]      2.0000000
      [ 1]      2.0000000
      [ 2]      1.0000000

      euclidean norm of x = 3.0000000
      \endverbatim
    */
    REAL two_norm() const
    {
      return sqrt(two_norm_2());
    }

    //! pretty-print output property: true = scientific, false = fixed point representation
    bool scientific() const
    {
      return bScientific;
    }

    /*!
      \brief scientific(true) is the default, scientific(false) switches to the fixed point representation

      \b Example:
      \code
      hdnum::Vector<double> x(3);
      x[0] = 2.0;
      x[1] = 2.0;
      x[2] = 1.0;
      std::cout << "x=" << x << std::endl;
      x.scientific(false); // set fixed point display mode
      std::cout << "x=" << x << std::endl;
      \endcode

      \b Output:
      \verbatim
      x=
      [ 0]  2.0000000e+00
      [ 1]  2.0000000e+00
      [ 2]  1.0000000e+00

      x=
      [ 0]      2.0000000
      [ 1]      2.0000000
      [ 2]      1.0000000
      \endverbatim
    */
    void scientific(bool b) const
    {
      bScientific=b;
    }

    std::size_t width() const
    {
      return nValueWidth;
    }

    void width (std::size_t i) const
    {
      nValueWidth=i;
    }

    std::size_t iwidth() const
    {
      return nIndexWidth;
    }

    void iwidth (std::size_t i) const
    {
      nIndexWidth=i;
    }

    std::size_t precision() const
    {
      return nValuePrecision;
    }

    void precision (std::size_t i) const
    {
      nValuePrecision=i;
    }

  };



  template<typename REAL>
  bool Vector<REAL>::bScientific = true;

  template<typename REAL>
  std::size_t Vector<REAL>::nIndexWidth = 2;

  template<typename REAL>
  std::size_t Vector<REAL>::nValueWidth = 15;

  template<typename REAL>
  std::size_t Vector<REAL>::nValuePrecision = 7;


  /*!
    \relates Vector
    \brief Output operator for Vector

    \b Example:
    \code
    hdnum::Vector<double> x(3);
    x[0] = 2.0;
    x[1] = 2.0;
    x[2] = 1.0;
    std::cout << "x=" << x << std::endl;
    \endcode

    \b Output:
    \verbatim
    x=
    [ 0]  2.0000000e+00
    [ 1]  2.0000000e+00
    [ 2]  1.0000000e+00
    \endverbatim
  */
  template <typename REAL>
  inline std::ostream & operator <<(std::ostream & os, const Vector<REAL> & x)
  {
    os << std::endl;

    for (size_t r = 0; r < x.size(); ++r)
      {
        if( x.scientific() )
          {
            os << "["
               << std::setw(x.iwidth())
               << r
               << "]"
               << std::scientific
               << std::showpoint
               << std::setw( x.width() )
               << std::setprecision( x.precision() )
               << x[r]
               << std::endl;
          }
        else
          {
            os << "["
               << std::setw(x.iwidth())
               << r
               << "]"
               << std::fixed
               << std::showpoint
               << std::setw( x.width() )
               << std::setprecision( x.precision() )
               << x[r]
               << std::endl;
          }
      }
    return os;
  }



  /*!
    \relates Vector
    \brief Output contents of a Vector x to a text file named fname

    \b Example:
    \code
    hdnum::Vector<double> x(5);
    unitvector(x,3);
    x.scientific(false); // set fixed point display mode
    gnuplot("test.dat",x);
    \endcode

    \b Output:
    \verbatim
    Contents of 'test.dat':
    0      0.0000000
    1      0.0000000
    2      0.0000000
    3      1.0000000
    4      0.0000000
    \endverbatim
  */
  template<typename REAL>
  inline void gnuplot(
                      const std::string& fname,
                      const Vector<REAL> x
                      )
  {
    std::fstream f(fname.c_str(),std::ios::out);
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      {
        if( x.scientific() )
          {
            f << std::setw(x.width())
              << i
              << std::scientific
              << std::showpoint
              << std::setw( x.width() )
              << std::setprecision( x.precision() )
              << x[i]
              << std::endl;
          }
        else
          {
            f << std::setw(x.width())
              << i
              << std::fixed
              << std::showpoint
              << std::setw( x.width() )
              << std::setprecision( x.precision() )
              << x[i]
              << std::endl;
          }
      }
    f.close();
  }

  //! gnuplot output for two Vectors
  template<typename REAL>
  inline void gnuplot(
                      const std::string& fname,
                      const Vector<REAL> x,
                      const Vector<REAL> y
                      )
  {
    std::fstream f(fname.c_str(),std::ios::out);
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      {
        if( x.scientific() )
          {
            f << std::setw(x.width())
              << i
              << std::scientific
              << std::showpoint
              << std::setw( x.width() )
              << std::setprecision( x.precision() )
              << x[i]
              << " "
              << std::setw( x.width() )
              << std::setprecision( x.precision() )
              << y[i]
              << std::endl;
          }
        else
          {
            f << std::setw(x.width())
              << i
              << std::fixed
              << std::showpoint
              << std::setw( x.width() )
              << std::setprecision( x.precision() )
              << x[i]
              << " "
              << std::setw( x.width() )
              << std::setprecision( x.precision() )
              << y[i]
              << std::endl;
          }
      }

    f.close();
  }



  /*!
    \relates Vector
    \brief Read vector from a text file

    \param[in] filename name of the text file
    \param[in,out] x reference to a Vector

    \b Example:
    \code
    hdnum::Vector<number> x;
    readVectorFromFile("x.dat", x );
    std::cout << "x=" << x << std::endl;
    \endcode

    \b Output:
    \verbatim
    Contents of "x.dat":
    1.0
    2.0
    3.0

    would give:
    x=
    [ 0]  1.0000000e+00
    [ 1]  2.0000000e+00
    [ 2]  3.0000000e+00
    \endverbatim
  */
  template<typename REAL>
  inline void readVectorFromFile (const std::string& filename, Vector<REAL> &vector)
  {
    std::string buffer;
    std::ifstream fin( filename.c_str() );
    if( fin.is_open() ){
      while( fin ){
        std::string sub;
        fin >> sub;
        //std::cout << " sub = " << sub.c_str() << ": ";
        if( sub.length()>0 ){
          REAL a = atof(sub.c_str());
          //std::cout << std::fixed << std::setw(10) << std::setprecision(5) << a;
          vector.push_back(a);
        }
      }
      fin.close();
    }
    else{
      HDNUM_ERROR("Could not open file!");
    }
  }


  //! annulize vector
  template<class REAL>
  inline void zero (Vector<REAL>& x)
  {
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      x[i] = REAL(0);
  }

  //! norm of a vector
  template<class REAL>
  inline REAL norm (Vector<REAL> x)
  {
    REAL sum(0.0);
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      sum += x[i]*x[i];
    return sqrt(sum);
  }

  //! fill vector, all with the same entry
  template<class REAL>
  inline void fill (Vector<REAL>& x, const REAL t)
  {
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      x[i] = t;
  }

  /*!
    \relates Vector
    \brief Fill vector, with entries starting at t, consecutively shifted by dt

    \b Example:
    \code
    hdnum::Vector<double> x(5);
    fill(x,2.01,0.1);
    x.scientific(false); // set fixed point display mode
    std::cout << "x=" << x << std::endl;
    \endcode

    \b Output:
    \verbatim
    x=
    [ 0]      2.0100000
    [ 1]      2.1100000
    [ 2]      2.2100000
    [ 3]      2.3100000
    [ 4]      2.4100000
    \endverbatim
  */
  template<class REAL>
  inline void fill (Vector<REAL>& x, const REAL& t, const REAL& dt)
  {
    REAL myt(t);
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      {
        x[i] = myt;
        myt += dt;
      }
  }


  /*!
    \relates Vector
    \brief Defines j-th unitvector (j=0,...,n-1) where n = length of the vector

    \b Example:
    \code
    hdnum::Vector<double> x(5);
    unitvector(x,3);
    x.scientific(false); // set fixed point display mode
    std::cout << "x=" << x << std::endl;
    \endcode

    \b Output:
    \verbatim
    x=
    [ 0]      0.0000000
    [ 1]      0.0000000
    [ 2]      0.0000000
    [ 3]      1.0000000
    [ 4]      0.0000000
    \endverbatim
  */
  template<class REAL>
  inline void unitvector (Vector<REAL> & x, std::size_t j)
  {
    for (typename Vector<REAL>::size_type i=0; i<x.size(); i++)
      if (i==j)
        x[i] = REAL(1);
      else
        x[i] = REAL(0);
  }


} // end of namespace hdnum

#endif	/* _VECTOR_HH */
