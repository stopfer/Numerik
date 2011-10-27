// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
/* 
 * File:   densematrix.hh
 * Author: ngo
 *
 * Created on April 15, 2011
 */

#ifndef DENSEMATRIX_HH
#define	DENSEMATRIX_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include "countablearray.hh"
#include "vector.hh"


namespace hdnum {

  /*! \brief Class with mathematical matrix operations
   */

  template<typename REAL>
  class DenseMatrix
  {
  public:
	/** \brief Type used for array indices */
	typedef std::size_t size_type;
	typedef typename std::vector<REAL> VType;
	typedef typename VType::const_iterator ConstVectorIterator;
	typedef typename VType::iterator VectorIterator;

  private:
	VType m_data;  // Matrix data is stored in an STL vector!
	std::size_t m_rows;          // Number of Matrix rows
	std::size_t m_cols;          // Number of Matrix columns

	static bool bScientific;
	static std::size_t nIndexWidth;
	static std::size_t nValueWidth;
	static std::size_t nValuePrecision;
	

	// function to calculate the modulus of a value
    REAL myabs (REAL x) const
    {
      if (x>=REAL(0)) 
		return x; 
	  else 
		return -x;
    }
	

	// get matrix element for write access:
	inline REAL & at(const std::size_t row, const std::size_t col)
	{
	  return m_data[row * m_cols + col];
	}
	
	// get matrix element for read-only access:
	inline const REAL & at(const std::size_t row, const std::size_t col) const
	{
	  return m_data[row * m_cols + col];
	}
	
  public:
	
	// default constructor (empty Matrix)
	DenseMatrix()
	  : m_data( 0, 0 )
	  , m_rows( 0 )
	  , m_cols( 0 )
	{
	}
	
	// constructor
	DenseMatrix( const std::size_t _rows, 
				 const std::size_t _cols, 
				 const REAL def_val=0 
				 )
	  : m_data( _rows*_cols, def_val )
	  , m_rows( _rows )
	  , m_cols( _cols )
	{
	}

    void addNewRow( const hdnum::Vector<REAL> & rowvector ){
      m_rows++;
      m_cols = rowvector.size();
      for(std::size_t i=0; i<m_cols; i++ )
        m_data.push_back( rowvector[i] );
    }
    /*
	// copy constructor (not needed, since it inherits from the STL vector)
	DenseMatrix( const DenseMatrix& A )
    {
    this->m_data = A.m_data;
    m_rows = A.m_rows;
    m_cols = A.m_cols;
    }
    */
	/*! 
      \brief get number of rows of the matrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(4,5);
      size_t nRows = A.rowsize();
      std::cout << "Matrix A has " << nRows << " rows." << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      Matrix A has 4 rows.
	  \endverbatim
    */
	size_t rowsize () const
	{
	  return m_rows;
	}

	/*! 
      \brief get number of columns of the matrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(4,5);
      size_t nColumns = A.colsize();
      std::cout << "Matrix A has " << nColumns << " columns." << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      Matrix A has 5 columns.
	  \endverbatim
    */
	size_t colsize () const
	{
	  return m_cols;
	}


	// pretty-print output properties
	bool scientific() const 
	{
	  return bScientific;
	}
	
	/*! 
      \brief   Switch between floating point (default=true) and fixed point (false) display

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(4,4);
      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(3);
      identity(A);  // Defines the identity matrix of the same dimension
      std::cout << "A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2        3 
      0     1.000    0.000    0.000    0.000 
      1     0.000    1.000    0.000    0.000 
      2     0.000    0.000    1.000    0.000 
      3     0.000    0.000    0.000    1.000 
	  \endverbatim
    */
    void scientific(bool b) const 
    {
      bScientific=b;
    }

    //! get index field width for pretty-printing
    std::size_t iwidth () const
    {
      return nIndexWidth;
    }

    //! get data field width for pretty-printing
    std::size_t width () const
    {
      return nValueWidth;
    }

    //! get data precision for pretty-printing
    std::size_t precision () const
    {
      return nValuePrecision;
    }

    //! set index field width for pretty-printing
    void iwidth (std::size_t i) const
    {
      nIndexWidth=i;
    }

    //! set data field width for pretty-printing
    void width (std::size_t i) const
    {
      nValueWidth=i;
    }

    //! set data precision for pretty-printing
    void precision (std::size_t i) const
    {
      nValuePrecision=i;
    }


	/*! 
      \brief (i,j)-operator for accessing entries of a (m x n)-matrix directly

      \param[in] i row index (0...m-1)
      \param[in] j column index (0...n-1)

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(4,4);
      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(3);

      identity(A);  // Defines the identity matrix of the same dimension
      std::cout << "A=" << A << std::endl;

      std::cout << "reading A(0,0)=" << A(0,0) << std::endl;

      std::cout << "resetting A(0,0) and A(2,3)..." << std::endl;
      A(0,0) = 1.234;
      A(2,3) = 432.1;
  
      std::cout << "A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2        3 
      0     1.000    0.000    0.000    0.000 
      1     0.000    1.000    0.000    0.000 
      2     0.000    0.000    1.000    0.000 
      3     0.000    0.000    0.000    1.000 

      reading A(0,0)=1.000
      resetting A(0,0) and A(2,3)...
      A=
      0        1        2        3 
      0     1.234    0.000    0.000    0.000 
      1     0.000    1.000    0.000    0.000 
      2     0.000    0.000    1.000  432.100 
      3     0.000    0.000    0.000    1.000 
	  \endverbatim
    */
	// overloaded element access operators
	// write access on matrix element A_ij using A(i,j)
	inline REAL & operator()(const std::size_t row, const std::size_t col)
	{
	  assert(row < m_rows|| col < m_cols);
	  return at(row,col);
	}
	
	// read-access on matrix element A_ij using A(i,j)
	inline const REAL & operator()(const std::size_t row, const std::size_t col) const
	{
	  assert(row < m_rows|| col < m_cols);
	  return at(row,col);
	}


	// read-access on matrix element A_ij using A[i][j]
	const ConstVectorIterator operator[](const std::size_t row) const
	{
	  assert(row < m_rows);
	  return m_data.begin() + row * m_cols;
	}

	/*! 
      \brief [i][j]-operator for accessing entries of a (m x n)-matrix directly

      \param[in] i row index (0...m-1)
      \param[in] j column index (0...n-1)

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(4,4);
      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(3);

      identity(A);  // Defines the identity matrix of the same dimension
      std::cout << "A=" << A << std::endl;

      std::cout << "reading A[0][0]=" << A[0][0] << std::endl;
      std::cout << "resetting A[0][0] and A[2][3]..." << std::endl;
      A[0][0] = 1.234;
      A[2][3] = 432.1;

      std::cout << "A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2        3 
      0     1.000    0.000    0.000    0.000 
      1     0.000    1.000    0.000    0.000 
      2     0.000    0.000    1.000    0.000 
      3     0.000    0.000    0.000    1.000 

      reading A[0][0]=1.000
      resetting A[0][0] and A[2][3]...
      A=
      0        1        2        3 
      0     1.234    0.000    0.000    0.000 
      1     0.000    1.000    0.000    0.000 
      2     0.000    0.000    1.000  432.100 
      3     0.000    0.000    0.000    1.000 
	  \endverbatim
    */

	// write-access on matrix element A_ij using A[i][j]
	VectorIterator operator[](const std::size_t row)
	{
	  assert(row < m_rows);
	  return m_data.begin() + row * m_cols;
	}



	/*! 
      \brief assignment operator

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(4,4);
      spd(A);
      hdnum::DenseMatrix<double> B(4,4);
      B = A;
      std::cout << "B=" << B << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      B=
      0          1          2          3 
      0   4.000e+00 -1.000e+00 -2.500e-01 -1.111e-01 
      1  -1.000e+00  4.000e+00 -1.000e+00 -2.500e-01 
      2  -2.500e-01 -1.000e+00  4.000e+00 -1.000e+00 
      3  -1.111e-01 -2.500e-01 -1.000e+00  4.000e+00 
	  \endverbatim
    */
    DenseMatrix& operator= (const DenseMatrix& A)
    {
	  m_data = A.m_data;
	  m_rows = A.m_rows;
	  m_cols = A.m_cols;
	  return *this;
    }



	/*! 
      \brief assignment from a scalar value

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(2,3);
      A = 5.432;
      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(3);
      std::cout << "A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2 
      0     5.432    5.432    5.432
      1     5.432    5.432    5.432
	  \endverbatim
    */
	DenseMatrix& operator= (const REAL value)
    {
      for (std::size_t i=0; i<rowsize(); i++)
        for (std::size_t j=0; j<colsize(); j++)
          (*this)(i,j) = value;
      return *this;
    }

    /*! 
      \brief Submatrix extraction

      Returns a new matrix that is a subset of the components
      of the given matrix.

      \param[in] i first row index of the new matrix
      \param[in] j first column index of the new matrix
      \param[in] rows row size of the new matrix, i.e. it has components [i,i+rows-1]
      \param[in] cols column size of the new matrix, i.e. it has components [j,j+m-1]
    */
    DenseMatrix sub (size_type i, size_type j, size_type rows, size_type cols)
    {
      DenseMatrix A(rows,cols);
      DenseMatrix &self = *this;
      size_type k1=0, k2=0;
      for (size_type i_=i; i_ < i+rows; i_++){
        for (size_type j_=j; j_ < j+cols; j_++){
          A[k1][k2] = self[i_][j_];
          k2++;
        }
        k1++;
      }
      return A;
    }



	// Basic Matrix Operations

	/*! 
      \brief Addition assignment 
	  
      Implements A += B matrix  addition
	  
      \param[in] B another Matrix
    */
    DenseMatrix& operator+= (const DenseMatrix& B)
    {
      for (std::size_t i=0; i<rowsize(); ++i) 
        for (std::size_t j=0; j<colsize(); ++j) 
		  (*this)(i,j) += B(i,j);
      return *this;
    }

	/*! 
      \brief Subtraction assignment
	  
      Implements A -= B matrix subtraction
	  
      \param[in] B another matrix
    */
    DenseMatrix& operator-= (const DenseMatrix& B)
    {
      for (std::size_t i=0; i<rowsize(); ++i) 
        for (std::size_t j=0; j<colsize(); ++j) 
		  (*this)(i,j) -= B(i,j);
      return *this;
    }


	/*! 
      \brief Scalar multiplication assignment
	  
      Implements A *= s where s is a scalar

      \param[in] s scalar value to multiply with

	  \b Example:
	  \code
	  double s = 0.5;
	  hdnum::DenseMatrix<double> A(2,3,1.0);
	  std::cout << "A=" << A << std::endl;
	  A *= s;
	  std::cout << "A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0          1          2 
      0   1.000e+00  1.000e+00  1.000e+00 
      1   1.000e+00  1.000e+00  1.000e+00 

      0.5*A =
      0          1          2 
      0   5.000e-01  5.000e-01  5.000e-01 
      1   5.000e-01  5.000e-01  5.000e-01 
	  \endverbatim
    */
    DenseMatrix& operator*= (const REAL s)
    {
      for (std::size_t i=0; i<rowsize(); ++i) 
        for (std::size_t j=0; j<colsize(); ++j) 
		  (*this)(i,j) *= s;
      return *this;
    }


	/*! 
      \brief Scalar division assignment
	  
      Implements A /= s where s is a scalar
	  
      \param[in] s scalar value to multiply with

	  \b Example:
	  \code
	  double s = 0.5;
	  hdnum::DenseMatrix<double> A(2,3,1.0);
	  std::cout << "A=" << A << std::endl;
	  A /= s;
	  std::cout << "A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0          1          2 
      0   1.000e+00  1.000e+00  1.000e+00 
      1   1.000e+00  1.000e+00  1.000e+00 

      A/0.5 =
      0          1          2 
      0   2.000e+00  2.000e+00  2.000e+00 
      1   2.000e+00  2.000e+00  2.000e+00 
	  \endverbatim

    */
    DenseMatrix& operator/= (const REAL s)
    {
      for (std::size_t i=0; i<rowsize(); ++i) 
        for (std::size_t j=0; j<colsize(); ++j) 
		  (*this)(i,j) /= s;
      return *this;
    }


    /*! 
      \brief Scaled update of a Matrix
	  
      Implements A += s*B where s is a scalar and B a matrix
	  
      \param[in] s scalar value to multiply with
      \param[in] B another matrix

	  \b Example:
	  \code
	  double s = 0.5;
	  hdnum::DenseMatrix<double> A(2,3,1.0);
	  hdnum::DenseMatrix<double> B(2,3,2.0);
	  A.update(s,B);
	  std::cout << "A + s*B =" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A + s*B =
      0          1          2 
      0       1.500      1.500      1.500 
      1       1.500      1.500      1.500 
	  \endverbatim

    */
	void update (const REAL s, const DenseMatrix& B)
    {
      for (std::size_t i=0; i<rowsize(); ++i) 
        for (std::size_t j=0; j<colsize(); ++j) 
		  (*this)(i,j) += s*B(i,j);
    }
	


    /*! 
      \brief matrix vector product y = A*x
	  
      Implements y = A*x where x and y are a vectors and A is a matrix
	  
      \param[in] y reference to the resulting Vector
      \param[in] x constant reference to a Vector

	  \b Example:
	  \code
      hdnum::Vector<double> x(3,10.0);
      hdnum::Vector<double> y(2);
      hdnum::DenseMatrix<double> A(2,3,1.0);

      x.scientific(false); // fixed point representation for all Vector objects
      A.scientific(false); // fixed point representation for all DenseMatrix objects

      std::cout << "A =" << A << std::endl;
      std::cout << "x =" << x << std::endl;
      A.mv(y,x);
      std::cout << "y = A*x =" << y << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A =
      0          1          2 
      0       1.000      1.000      1.000 
      1       1.000      1.000      1.000 

      x =
      [ 0]     10.0000000
      [ 1]     10.0000000
      [ 2]     10.0000000

      y = A*x =
      [ 0]     30.0000000
      [ 1]     30.0000000
	  \endverbatim

    */
    template<class V>
    void mv (Vector<V>& y, const Vector<V>& x) const
    {
      if (this->rowsize()!=y.size()) 
        HDNUM_ERROR("mv: size of A and y do not match");
      if (this->colsize()!=x.size()) 
        HDNUM_ERROR("mv: size of A and x do not match");
      for (std::size_t i=0; i<rowsize(); ++i) 
        {
          y[i] = 0;
		  for (std::size_t j=0; j<colsize(); ++j) 
            y[i] += (*this)(i,j)*x[j];
        }
    }


    /*! 
      \brief update matrix vector product y += A*x
	  
      Implements y += A*x where x and y are a vectors and A is a matrix
	  
      \param[in] y reference to the resulting Vector
      \param[in] x constant reference to a Vector

	  \b Example:
	  \code
      hdnum::Vector<double> x(3,10.0);
      hdnum::Vector<double> y(2,5.0);
      hdnum::DenseMatrix<double> A(2,3,1.0);

      x.scientific(false); // fixed point representation for all Vector objects
      A.scientific(false); // fixed point representation for all DenseMatrix objects

      std::cout << "y =" << y << std::endl;
      std::cout << "A =" << A << std::endl;
      std::cout << "x =" << x << std::endl;
      A.umv(y,x);
      std::cout << "y = A*x =" << y << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      y =
      [ 0]      5.0000000
      [ 1]      5.0000000

      A =
      0          1          2 
      0       1.000      1.000      1.000 
      1       1.000      1.000      1.000 

      x =
      [ 0]     10.0000000
      [ 1]     10.0000000
      [ 2]     10.0000000

      y + A*x =
      [ 0]     35.0000000
      [ 1]     35.0000000
	  \endverbatim

    */
    template<class V>
    void umv (Vector<V>& y, const Vector<V>& x) const
    {
      if (this->rowsize()!=y.size()) 
        HDNUM_ERROR("mv: size of A and y do not match");
      if (this->colsize()!=x.size()) 
        HDNUM_ERROR("mv: size of A and x do not match");
      for (std::size_t i=0; i<rowsize(); ++i) 
        {
		  for (std::size_t j=0; j<colsize(); ++j) 
            y[i] += (*this)(i,j)*x[j];
        }
    }


    /*! 
      \brief update matrix vector product y += sA*x
	  
      Implements y += sA*x where s is a scalar value, x and y are a vectors and A is a matrix
	  
      \param[in] y reference to the resulting Vector
      \param[in] s constant reference to a number type
      \param[in] x constant reference to a Vector

	  \b Example:
	  \code
      double s=0.5;
      hdnum::Vector<double> x(3,10.0);
      hdnum::Vector<double> y(2,5.0);
      hdnum::DenseMatrix<double> A(2,3,1.0);

      x.scientific(false); // fixed point representation for all Vector objects
      A.scientific(false); // fixed point representation for all DenseMatrix objects

      std::cout << "y =" << y << std::endl;
      std::cout << "A =" << A << std::endl;
      std::cout << "x =" << x << std::endl;
      A.umv(y,s,x);
      std::cout << "y = s*A*x =" << y << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      y =
      [ 0]      5.0000000
      [ 1]      5.0000000

      A =
      0          1          2 
      0       1.000      1.000      1.000 
      1       1.000      1.000      1.000 

      x =
      [ 0]     10.0000000
      [ 1]     10.0000000
      [ 2]     10.0000000

      y = s*A*x =
      [ 0]     20.0000000
      [ 1]     20.0000000
	  \endverbatim

    */
    template<class V>
	void umv (Vector<V>& y, const V& s, const Vector<V>& x) const
    {
      if (this->rowsize()!=y.size()) 
        HDNUM_ERROR("mv: size of A and y do not match");
      if (this->colsize()!=x.size()) 
        HDNUM_ERROR("mv: size of A and x do not match");
      for (std::size_t i=0; i<rowsize(); ++i) 
        {
		  for (std::size_t j=0; j<colsize(); ++j) 
            y[i] += s*(*this)(i,j)*x[j];
        }
    }



    /*! 
      \brief assign to matrix product C = A*B to matrix C
	  
      Implements C = A*B where A and B are matrices
	  
      \param[in] x constant reference to a DenseMatrix
      \param[in] x constant reference to a DenseMatrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(2,6,1.0);
      hdnum::DenseMatrix<double> B(6,3,-1.0);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(6);          // use at least 6 columns for displaying matrix entries
      A.precision(3);      // display 3 digits behind the point 

      std::cout << "A =" << A << std::endl;
      std::cout << "B =" << B << std::endl;

      hdnum::DenseMatrix<double> C(2,3);
      C.mm(A,B);
      std::cout << "C = A*B =" << C << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A =
      0      1      2      3      4      5 
      0   1.000  1.000  1.000  1.000  1.000  1.000 
      1   1.000  1.000  1.000  1.000  1.000  1.000 

      B =
      0      1      2 
      0  -1.000 -1.000 -1.000 
      1  -1.000 -1.000 -1.000 
      2  -1.000 -1.000 -1.000 
      3  -1.000 -1.000 -1.000 
      4  -1.000 -1.000 -1.000 
      5  -1.000 -1.000 -1.000 

      C = A*B =
      0      1      2 
      0  -6.000 -6.000 -6.000 
      1  -6.000 -6.000 -6.000 
	  \endverbatim

    */
    void mm (const DenseMatrix<REAL>& A, const DenseMatrix<REAL>& B)
    {
      if (this->rowsize()!=A.rowsize()) 
        HDNUM_ERROR("mm: size incompatible");
      if (this->colsize()!=B.colsize()) 
        HDNUM_ERROR("mm: size incompatible");
      if (A.colsize()!=B.rowsize()) 
        HDNUM_ERROR("mm: size incompatible");

      for (std::size_t i=0; i<rowsize(); i++)
        for (std::size_t j=0; j<colsize(); j++)
          {
            (*this)(i,j) = 0;
            for (std::size_t k=0; k<A.colsize(); k++)
              (*this)(i,j) += A(i,k)*B(k,j);
          }
    }


    /*! 
      \brief add matrix product A*B to matrix C
	  
      Implements C += A*B where A, B and C are matrices
	  
      \param[in] x constant reference to a DenseMatrix
      \param[in] x constant reference to a DenseMatrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(2,6,1.0);
      hdnum::DenseMatrix<double> B(6,3,-1.0);
      hdnum::DenseMatrix<double> C(2,3,0.5);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(6);
      A.precision(3);

      std::cout << "C =" << C << std::endl;
      std::cout << "A =" << A << std::endl;
      std::cout << "B =" << B << std::endl;

      C.umm(A,B);
      std::cout << "C + A*B =" << C << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      C =
      0      1      2 
      0   0.500  0.500  0.500 
      1   0.500  0.500  0.500 

      A =
      0      1      2      3      4      5 
      0   1.000  1.000  1.000  1.000  1.000  1.000 
      1   1.000  1.000  1.000  1.000  1.000  1.000 

      B =
      0      1      2 
      0  -1.000 -1.000 -1.000 
      1  -1.000 -1.000 -1.000 
      2  -1.000 -1.000 -1.000 
      3  -1.000 -1.000 -1.000 
      4  -1.000 -1.000 -1.000 
      5  -1.000 -1.000 -1.000 

      C + A*B =
      0      1      2 
      0  -5.500 -5.500 -5.500 
      1  -5.500 -5.500 -5.500 
	  \endverbatim

    */
    void umm (const DenseMatrix<REAL>& A, const DenseMatrix<REAL>& B)
    {
      if (this->rowsize()!=A.rowsize()) 
        HDNUM_ERROR("mm: size incompatible");
      if (this->colsize()!=B.colsize()) 
        HDNUM_ERROR("mm: size incompatible");
      if (A.colsize()!=B.rowsize()) 
        HDNUM_ERROR("mm: size incompatible");

      for (std::size_t i=0; i<rowsize(); i++)
        for (std::size_t j=0; j<colsize(); j++)
		  for (std::size_t k=0; k<A.colsize(); k++)
            (*this)(i,j) += A(i,k)*B(k,j);
    }



    /*! 
      \brief set column: make x the k'th column of A
	  	  
      \param[in] x constant reference to a Vector
      \param[in] k number of the column of A to be set

	  \b Example:
	  \code
      hdnum::Vector<double> x(2,434.0);
      hdnum::DenseMatrix<double> A(2,6);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(1);
  
      std::cout << "original A=" << A << std::endl;
      A.sc(x,3);   // redefine fourth column of the matrix
      std::cout << "modified A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      original A=
      0        1        2        3        4        5 
      0       0.0      0.0      0.0      0.0      0.0      0.0 
      1       0.0      0.0      0.0      0.0      0.0      0.0 

      modified A=
      0        1        2        3        4        5 
      0       0.0      0.0      0.0    434.0      0.0      0.0 
      1       0.0      0.0      0.0    434.0      0.0      0.0 
	  \endverbatim

    */
    void sc (const Vector<REAL>& x, std::size_t k)
    {
      if (this->rowsize()!=x.size()) 
        HDNUM_ERROR("cc: size incompatible");

      for (std::size_t i=0; i<rowsize(); i++)
        (*this)(i,k) = x[i];
    }

    /*! 
      \brief set row: make x the k'th row of A
	  	  
      \param[in] x constant reference to a Vector
      \param[in] k number of the row of A to be set

	  \b Example:
	  \code
      hdnum::Vector<double> x(3,434.0);
      hdnum::DenseMatrix<double> A(3,3);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(1);
  
      std::cout << "original A=" << A << std::endl;
      A.sr(x,1);   // redefine second row of the matrix
      std::cout << "modified A=" << A << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      original A=
      0        1        2 
      0       0.0      0.0      0.0 
      1       0.0      0.0      0.0 
      2       0.0      0.0      0.0 

      modified A=
      0        1        2 
      0       0.0      0.0      0.0 
      1     434.0    434.0    434.0 
      2       0.0      0.0      0.0 
	  \endverbatim

    */
    void sr (const Vector<REAL>& x, std::size_t k)
    {
      if (this->colsize()!=x.size()) 
        HDNUM_ERROR("cc: size incompatible");

      for (std::size_t i=0; i<colsize(); i++)
        (*this)(k,i) = x[i];
    }


    //! compute row sum norm
    REAL norm_infty () const
    {
      REAL norm(0.0);
      for (std::size_t i=0; i<rowsize(); i++)
        {
          REAL sum(0.0);
          for (std::size_t j=0; j<colsize(); j++)
            sum += myabs((*this)(i,j));
          if (sum>norm) norm = sum;
        }
      return norm;
    }

    //! compute column sum norm
    REAL norm_1 () const
    {
      REAL norm(0.0);
      for (std::size_t j=0; j<colsize(); j++)
        {
          REAL sum(0.0);
          for (std::size_t i=0; i<rowsize(); i++)
            sum += myabs((*this)(i,j));
          if (sum>norm) norm = sum;
        }
      return norm;
    }
	


    /*! 
      \brief vector = matrix * vector
	  	  
      \param[in] x constant reference to a Vector

	  \b Example:
	  \code
      hdnum::Vector<double> x(3,4.0);
      hdnum::DenseMatrix<double> A(3,3,2.0);
      hdnum::Vector<double> y(3);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(1);
  
      x.scientific(false); // fixed point representation for all Vector objects
      x.width(8);
      x.precision(1);


      std::cout << "A=" << A << std::endl;
      std::cout << "x=" << x << std::endl;
      y=A*x;
      std::cout << "y=A*x" << y << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2 
      0       2.0      2.0      2.0 
      1       2.0      2.0      2.0 
      2       2.0      2.0      2.0 

      x=
      [ 0]     4.0
      [ 1]     4.0
      [ 2]     4.0

      y=A*x
      [ 0]    24.0
      [ 1]    24.0
      [ 2]    24.0
	  \endverbatim

    */
	Vector<REAL> operator* (const Vector<REAL> & x)
	{
	  assert( x.size() == rowsize() );
	  
	  Vector<REAL> y( rowsize() );
	  for(std::size_t r=0; r<rowsize(); ++r){
		for(std::size_t c=0; c<colsize(); ++c){
		  y[r]+= at(r,c) * x[c];
		}
	  }
	  return y;
	}
	
	

	
    /*! 
      \brief matrix = matrix * matrix
	  	  
      \param[in] x constant reference to a DenseMatrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(3,3,2.0);
      hdnum::DenseMatrix<double> B(3,3,4.0);
      hdnum::DenseMatrix<double> C(3,3);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(1);

      std::cout << "A=" << A << std::endl;
      std::cout << "B=" << B << std::endl;
      C=A*B;
      std::cout << "C=A*B=" << C << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2 
      0       2.0      2.0      2.0 
      1       2.0      2.0      2.0 
      2       2.0      2.0      2.0 

      B=
      0        1        2 
      0       4.0      4.0      4.0 
      1       4.0      4.0      4.0 
      2       4.0      4.0      4.0 

      C=A*B=
      0        1        2 
      0      24.0     24.0     24.0 
      1      24.0     24.0     24.0 
      2      24.0     24.0     24.0 
	  \endverbatim

    */
	DenseMatrix operator* (const DenseMatrix & x) const
	{
	  assert(colsize() == x.rowsize());
	  
	  const std::size_t out_rows = rowsize();
	  const std::size_t out_cols = x.colsize();
	  DenseMatrix y(out_rows, out_cols,0.0);
	  for(std::size_t r=0; r<out_rows; ++r)
		for(std::size_t c=0; c<out_cols; ++c)
		  for(std::size_t i=0; i<colsize(); ++i)
			y(r,c) += at(r,i) * x(i,c);
	  
	  return y;
	}



    /*! 
      \brief matrix = matrix + matrix
	  	  
      \param[in] x constant reference to a DenseMatrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(3,3,2.0);
      hdnum::DenseMatrix<double> B(3,3,4.0);
      hdnum::DenseMatrix<double> C(3,3);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(1);

      std::cout << "A=" << A << std::endl;
      std::cout << "B=" << B << std::endl;
      C=A+B;
      std::cout << "C=A+B=" << C << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2 
      0       2.0      2.0      2.0 
      1       2.0      2.0      2.0 
      2       2.0      2.0      2.0 

      B=
      0        1        2 
      0       4.0      4.0      4.0 
      1       4.0      4.0      4.0 
      2       4.0      4.0      4.0 

      C=A+B=
      0        1        2 
      0       6.0      6.0      6.0 
      1       6.0      6.0      6.0 
      2       6.0      6.0      6.0 
	  \endverbatim

    */
	DenseMatrix operator+ (const DenseMatrix & x) const
	{
	  assert(colsize() == x.rowsize());
	  
	  const std::size_t out_rows = rowsize();
	  const std::size_t out_cols = x.colsize();
	  DenseMatrix y(out_rows, out_cols,0.0);
      y = *this;
      y+=x;
	  return y;
	}



    /*! 
      \brief matrix = matrix - matrix
	  	  
      \param[in] x constant reference to a DenseMatrix

	  \b Example:
	  \code
      hdnum::DenseMatrix<double> A(3,3,2.0);
      hdnum::DenseMatrix<double> B(3,3,4.0);
      hdnum::DenseMatrix<double> C(3,3);

      A.scientific(false); // fixed point representation for all DenseMatrix objects
      A.width(8);
      A.precision(1);

      std::cout << "A=" << A << std::endl;
      std::cout << "B=" << B << std::endl;
      C=A-B;
      std::cout << "C=A-B=" << C << std::endl;
	  \endcode

	  \b Output:
	  \verbatim
      A=
      0        1        2 
      0       2.0      2.0      2.0 
      1       2.0      2.0      2.0 
      2       2.0      2.0      2.0 

      B=
      0        1        2 
      0       4.0      4.0      4.0 
      1       4.0      4.0      4.0 
      2       4.0      4.0      4.0 

      C=A-B=
      0        1        2 
      0      -2.0     -2.0     -2.0 
      1      -2.0     -2.0     -2.0 
      2      -2.0     -2.0     -2.0 
	  \endverbatim

    */
	DenseMatrix operator- (const DenseMatrix & x) const
	{
	  assert(colsize() == x.rowsize());
	  
	  const std::size_t out_rows = rowsize();
	  const std::size_t out_cols = x.colsize();
	  DenseMatrix y(out_rows, out_cols,0.0);
      y = *this;
      y-=x;
	  return y;
	}


  };



  template<typename REAL>
  bool DenseMatrix<REAL>::bScientific = true;
  template<typename REAL>
  std::size_t DenseMatrix<REAL>::nIndexWidth = 10;
  template<typename REAL>
  std::size_t DenseMatrix<REAL>::nValueWidth = 10;
  template<typename REAL>
  std::size_t DenseMatrix<REAL>::nValuePrecision = 3;

  
  /*! 
    \brief   Output operator for Matrix

    \b Example:
    \code
    hdnum::DenseMatrix<double> A(4,4);
    A.scientific(false); // fixed point representation for all DenseMatrix objects
    A.width(8);
    A.precision(3);
    identity(A);  // Defines the identity matrix of the same dimension
    std::cout << "A=" << A << std::endl;
    \endcode

    \b Output:
    \verbatim
    A=
    0        1        2        3 
    0     1.000    0.000    0.000    0.000 
    1     0.000    1.000    0.000    0.000 
    2     0.000    0.000    1.000    0.000 
    3     0.000    0.000    0.000    1.000 
    \endverbatim
  */
  template <typename REAL>
  inline std::ostream& operator<< (std::ostream& s, const DenseMatrix<REAL>& A)
  {
	s << std::endl;
	s << " " << std::setw(A.iwidth()) << " " << "  ";
	for (typename DenseMatrix<REAL>::size_type j=0; j<A.colsize(); ++j)
	  s << std::setw(A.width()) << j << " ";
	s << std::endl;

	for (typename DenseMatrix<REAL>::size_type i=0; i<A.rowsize(); ++i)
	  {
		s << " " << std::setw(A.iwidth()) << i << "  ";
		for (typename DenseMatrix<REAL>::size_type j=0; j<A.colsize(); ++j)
		  {
			if( A.scientific() )
			  {
				s << std::setw(A.width()) 
				  << std::scientific 
				  << std::showpoint 
				  << std::setprecision(A.precision()) 
				  << A[i][j] << " ";
			  }
			else
			  {
				s << std::setw(A.width()) 
				  << std::fixed
				  << std::showpoint 
				  << std::setprecision(A.precision()) 
				  << A[i][j] << " ";
			  }
		  }
		s << std::endl;
	  }
    return s;
  }
  
  /*!
    make a matrix filled with one scalar value
    
    \param[in] A reference to a DenseMatrix that shall be filled with entries
    \param[in] t scalar value
  */
  template<typename REAL>
  inline void fill (DenseMatrix<REAL> A, const REAL& t)
  {
    for (typename DenseMatrix<REAL>::size_type i=0; i<A.rowsize(); ++i)
      for (typename DenseMatrix<REAL>::size_type j=0; j<A.colsize(); ++j)
        A[i][j] = t;
  }

  //! make a zero matrix
  template<typename REAL>
  inline void zero (DenseMatrix<REAL> &A)
  {
    for (std::size_t i=0; i<A.rowsize(); ++i)
	  for (std::size_t j=0; j<A.colsize(); ++j)
        A(i,j) = REAL(0);
  }

  /*! 
    \relates DenseMatrix
	\n	
	\b Function: make identity matrix
	\code
    template<class T>
    inline void identity (DenseMatrix<T> &A)
	\endcode
	\param[in] A reference to a DenseMatrix that shall be filled with entries
	
	\b Example:
	\code
    hdnum::DenseMatrix<double> A(4,4);
    identity(A);

    A.scientific(false); // fixed point representation for all DenseMatrix objects
    A.width(10);
    A.precision(5);

    std::cout << "A=" << A << std::endl;
	\endcode
	
	\b Output:
	\verbatim
    A=
    0          1          2          3 
    0     1.00000    0.00000    0.00000    0.00000 
    1     0.00000    1.00000    0.00000    0.00000 
    2     0.00000    0.00000    1.00000    0.00000 
    3     0.00000    0.00000    0.00000    1.00000 
	\endverbatim
	
  */
  template<class T>
  inline void identity (DenseMatrix<T> &A)
  {
    for (typename DenseMatrix<T>::size_type i=0; i<A.rowsize(); ++i)
      for (typename DenseMatrix<T>::size_type j=0; j<A.colsize(); ++j)
        if (i==j)
          A[i][i] = T(1);
        else
          A[i][j] = T(0);
  }

  /*! 
    \relates DenseMatrix
	
	\n
	\b Function: make a symmetric and positive definite matrix
	\code
    template<typename REAL>
    inline void spd (DenseMatrix<REAL> &A)
	\endcode
	
	\param[in] A reference to a DenseMatrix that shall be filled with entries
	
	\b Example:
	\code
    hdnum::DenseMatrix<double> A(4,4);
    spd(A);

    A.scientific(false); // fixed point representation for all DenseMatrix objects
    A.width(10);
    A.precision(5);

    std::cout << "A=" << A << std::endl;
	\endcode
	
	\b Output:
	\verbatim
    A=
    0          1          2          3 
    0     4.00000   -1.00000   -0.25000   -0.11111 
    1    -1.00000    4.00000   -1.00000   -0.25000 
    2    -0.25000   -1.00000    4.00000   -1.00000 
    3    -0.11111   -0.25000   -1.00000    4.00000 
	\endverbatim
	
  */
  template<typename REAL>
  inline void spd (DenseMatrix<REAL> &A)
  {
	if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
	  HDNUM_ERROR("need square and nonempty matrix");
    for (std::size_t i=0; i<A.rowsize(); ++i)
	  for (std::size_t j=0; j<A.colsize(); ++j)
		if (i==j)
		  A(i,i) = REAL(4.0);
		else
		  A(i,j) = - REAL(1.0)/((i-j)*(i-j));
  }

  /*!
    \relates DenseMatrix
	\n
	\b Function:  make a vandermonde matrix
	\code
    template<typename REAL>
    inline void vandermonde (DenseMatrix<REAL> &A, const Vector<REAL> x)
	\endcode
	
	\param[in] A reference to a DenseMatrix that shall be filled with entries
	\param[in] x constant reference to a Vector
	
	\b Example:
	\code
    hdnum::Vector<double> x(4);
    fill(x,2.0,1.0);
    hdnum::DenseMatrix<double> A(4,4);
    vandermonde(A,x);

    A.scientific(false); // fixed point representation for all DenseMatrix objects
    A.width(10);
    A.precision(5);

    x.scientific(false); // fixed point representation for all Vector objects
    x.width(10);
    x.precision(5);

    std::cout << "x=" << x << std::endl;
    std::cout << "A=" << A << std::endl;
	\endcode
	
	\b Output:
	\verbatim
    x=
    [ 0]   2.00000
    [ 1]   3.00000
    [ 2]   4.00000
    [ 3]   5.00000

    A=
    0          1          2          3 
    0     1.00000    2.00000    4.00000    8.00000 
    1     1.00000    3.00000    9.00000   27.00000 
    2     1.00000    4.00000   16.00000   64.00000 
    3     1.00000    5.00000   25.00000  125.00000 

	\endverbatim
	
  */
  template<typename REAL>
  inline void vandermonde (DenseMatrix<REAL> &A, const Vector<REAL> x)
  {
	if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
	  HDNUM_ERROR("need square and nonempty matrix");
	if (A.rowsize()!=x.size()) 
	  HDNUM_ERROR("need A and x of same size");
    for (typename DenseMatrix<REAL>::size_type i=0; i<A.rowsize(); ++i)
      {
        REAL p(1.0);
        for (typename DenseMatrix<REAL>::size_type j=0; j<A.colsize(); ++j)
          {
            A[i][j] = p;
            p *= x[i];
          }
      }
  }
  
  //! gnuplot output for matrix
  template<typename REAL>
  inline void gnuplot (const std::string& fname, const DenseMatrix<REAL> &A)
  {
    std::fstream f(fname.c_str(),std::ios::out); 
    for (typename DenseMatrix<REAL>::size_type i=0; i<A.rowsize(); ++i)
      {
        for (typename DenseMatrix<REAL>::size_type j=0; j<A.colsize(); ++j)
		  {
			if( A.scientific() )
			  {
				f << std::setw(A.width()) 
				  << std::scientific 
				  << std::showpoint 
				  << std::setprecision(A.precision()) << A[i][j];
			  }
			else
			  {
				f << std::setw(A.width()) 
				  << std::fixed
				  << std::showpoint 
				  << std::setprecision(A.precision()) << A[i][j];
			  }
		  }
		f << std::endl;
      }
    f.close();
  }

  

  /*!
    \relates DenseMatrix
    \brief Read matrix from a text file	

    \param[in] filename name of the text file
    \param[in,out] A reference to a DenseMatrix

	\b Example:
	\code
    hdnum::DenseMatrix<number> L;
    readMatrixFromFile("matrixL.dat", L );
    std::cout << "L=" << L << std::endl;
	\endcode
	
	\b Output:
	\verbatim
    Contents of "matrixL.dat":
    1.000e+00  0.000e+00  0.000e+00	 
    2.000e+00  1.000e+00  0.000e+00	 
    3.000e+00  2.000e+00  1.000e+00	 

    would give:
    L=
    0          1          2
    0   1.000e+00  0.000e+00  0.000e+00
    1   2.000e+00  1.000e+00  0.000e+00
    2   3.000e+00  2.000e+00  1.000e+00
	\endverbatim
  */
  template<typename REAL>
  inline void readMatrixFromFile (const std::string& filename, DenseMatrix<REAL> &A)
  {
    std::string buffer;
    std::ifstream fin( filename.c_str() );
    std::size_t i=0;
    std::size_t j=0;
    if( fin.is_open() ){
      while( std::getline( fin, buffer ) ){
        std::istringstream iss(buffer);
        hdnum::Vector<REAL> rowvector;
        while( iss ){
          std::string sub;
          iss >> sub;
          //std::cout << " sub = " << sub.c_str() << ": ";
          if( sub.length()>0 ){
            REAL a = atof(sub.c_str());
            //std::cout << std::fixed << std::setw(10) << std::setprecision(5) << a;
            rowvector.push_back(a);
          }
          j++;
        }
        if( rowvector.size()>0 ){
          A.addNewRow( rowvector );
          i++;
          //std::cout << std::endl;
        }
      }
      fin.close();
    }
    else{
      HDNUM_ERROR("Could not open file!");
    }
  }

}

#endif	// DENSEMATRIX_HH

