// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_COUNTABLEARRAY_HH
#define HDNUM_COUNTABLEARRAY_HH

#include "countingptr.hh"
#include "array.hh"

/** @file
 *  @brief This file implements a basic dynamic array class
 */

namespace hdnum {

  /** \brief Dynamic array that can be used with the reference counting pointer

      Provides a dyamically allocated array with access operator, resizing
      and size method.
   */
  template<class T>
  class CountableArray : public Countable, public Array<T>
  {
  public:
	/** \brief Type used for array indices */
	typedef std::size_t size_type;

	//! make empty array
	CountableArray () : Countable(), Array<T>()
	{}

	//! make array with _n uninitialized components
	CountableArray (size_type _n) : Countable(), Array<T>(_n)
	{}

	//! make array with _n initialized components
	CountableArray (size_type _n, const T& _t) : Countable(), Array<T>(_n,_t)
	{}
  };

} // namespace hdnum

#endif
