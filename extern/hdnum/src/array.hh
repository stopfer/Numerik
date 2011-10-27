// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_ARRAY_HH
#define HDNUM_ARRAY_HH

/** @file
 *  @brief This file implements a basic dynamic array class
 */

namespace hdnum {

  /** \brief A basic dynamic array class

      Provides a dyamically allocated array with access operator, resizing
      and size method.
   */
  template<class T>
  class Array {
  public:

	//! Remember the storage type 
	typedef T value_type;

	/** \brief Reference to an object */
	typedef value_type& reference;

	/** \brief Const reference to an object */
	typedef const value_type& const_reference;

	/** \brief Type used for array indices */
	typedef std::size_t size_type;

	/** \brief Difference type */
	typedef std::ptrdiff_t difference_type;

	//! make empty array
	Array ()
	{
	  n = 0;
	  p = 0;
	}

	//! make array with _n uninitialized components
	Array (size_type _n) : n(_n)
	{
	  if (n>0) 
		p = new T[n];
	  else
		{
		  n = 0;
		  p = 0;
		}
	}

	//! make array with _n initialized components
	Array (size_type _n, const T& _t) : n(_n)
	{
	  if (n>0) 
		p = new T[n];
	  else
		{
		  n = 0;
		  p = 0;
		}
	  for (size_type i=0; i<n; i++) p[i] = _t;
	}

	//! copy constructor
	Array (const Array& a) 
	{
	  // allocate memory with same size as a
	  n = a.n;
	  if (n>0) 
		p = new T[n];
	  else
		{
		  n = 0;
		  p = 0;
		}

	  // and copy elements
	  for (size_type i=0; i<n; i++) p[i]=a.p[i];
	}

	//! destructor, free dynamic memory
	~Array () 
	{ 
	  if (n>0) delete [] p; 
	}

	//! reallocate array to given size, any data is lost
	void resize (size_type _n)
	{
	  if (n==_n) return;
	  if (n>0) delete [] p; 
	  n = _n;
	  if (n>0) 
		p = new T[n];
	  else
		{
		  n = 0;
		  p = 0;
		}
	}

	//! reallocate array to given size, any data is lost
	void resize (size_type _n, const T& _t)
	{
	  if (n==_n) return;
	  if (n>0) delete [] p; 
	  n = _n;
	  if (n>0) 
		p = new T[n];
	  else
		{
		  n = 0;
		  p = 0;
		}
	  for (size_type i=0; i<n; i++) p[i] = _t;
	}

	//! assignment
	Array& operator= (const Array& a)
	{
	  if (&a!=this) // check if this and a are different objects
		{
		  // adjust size of array
		  if (n!=a.n) // check if size is different
			{
			  if (n>0) delete [] p; // delete old memory
			  n = a.n;
			  if (n>0) 
				p = new T[n];
			  else
				{
				  n = 0;
				  p = 0;
				}
			}
		  // copy data
		  for (size_type i=0; i<this->n; i++) p[i]=a.p[i];
		}
	  return *this;
	}

	//! Component access
	reference operator[] (size_type i)
	{
	  return p[i];
	}

	//! Const component access
	const_reference operator[] (size_type i) const
	{
	  return p[i];
	}

	//! get array size
	size_type size () const
	{
	  return n;
	}

  private:
	size_type n;
	T* p;
  };

} // namespace hdnum

#endif
