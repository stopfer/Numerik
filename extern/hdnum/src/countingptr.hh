// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_COUNTINGPTR_HH
#define HDNUM_COUNTINGPTR_HH

#include <iostream>

/** @file
 *  @brief This file implements a counting pointer with configurable memory
 *         management policy
 *  Adapted from dune-pdelab
 */

namespace hdnum {

	/** @brief Don't delete target if reference count reaches zero
     *
     *  If this class is given to CP as the memory management policy, the CP
     *  objects won't delete the pointed to object if the reference count
     *  reaches zero.
     */
    class NondeletingMemoryManagementPolicy
	{
	public:
	  template<typename T>
	  static void delete_action (T* p)
	  {}
	};

	/** @brief Delete target if reference count reaches zero
     *
     *  If this class is given to CP as the memory management policy, the CP
     *  objects will delete the pointed to object if the reference count
     *  reaches zero.
     */
	class DeletingMemoryManagementPolicy
	{
	public:
	  template<typename T>
	  static void delete_action (T* p)
	  {
		if (p->reference_counter_zero())
		  delete p;
	  }
	};

	template<typename T, typename P=DeletingMemoryManagementPolicy>
    /** @brief Pointer with a reference count in the pointed-to object
     *
     *  @tparam T The type of the pointed-to object.  Must be derived from
     *            Countable.
     *  @tparam P What to do when the reference count reaches 0.  Two
     *            predefined policy classes are available:
     *            NondeletingMemoryManagementPolicy and
     *            DeletingMemoryManagementPolicy (the default).
     *
     *  An object cp of class CP points to another object of a class derived
     *  from Countable, or to 0.  If it does not point to 0, it will keep
     *  track of how many CP objects point to the same object.  If cp stops
     *  pointing to the target object, it will decrement its reference count,
     *  and if the reference count reaches zero may or may not delete the
     *  target object, depending on what the memory managment policy dictates.
     *
     *  cp may be set via assingment from an apropriate C pointer or another
     *  CP of the same type.  To access the pointed to object, the expressions
     *  *cp and cp->member may be used, where member is a member of the
     *  pointed to object.  Finally, CP objects may be compared using == and
     *  != to find out whether they point to the same object.
     */
	class CP
	{
	  T* p;

	public:
      //! Construct a CP object which points to 0
	  CP ()
	  {
		p = 0;
	  }

      //! Construct a CP object which points to p_ (which may be 0)
	  explicit CP (T* p_)
	  {
		p = p_;
		if (p!=0)
		  p->reference_counter_increment();
	  }

      //! Copy constructor
	  CP (const CP<T>& cp)
	  {
		p = cp.p;
		if (p!=0)
		  p->reference_counter_increment();
	  }

      //! Destructor
	  ~CP ()
	  {
		if (p!=0)
		  {
			p->reference_counter_decrement();
			P::delete_action(p);
		  }
	  }

      //! assignment from a C pointer
	  CP<T>& operator= (T* p_)
	  {
		if (p!=p_)
		  {
			if (p!=0)
			  p->reference_counter_decrement();
			p = p_;
			if (p!=0)
			  p->reference_counter_increment();
		  }
		return *this;
	  }

      //! copy operator
	  CP<T>& operator= (const CP<T>& cp)
	  {
		if (p!=cp.p)
		  {
			if (p!=0)
			  p->reference_counter_decrement();
			p = cp.p;
			if (p!=0)
			  p->reference_counter_increment();
		  }
		return *this;		
	  }

      //! target element access
	  T* operator-> () const
	  {
		return p;
	  }

      //! dereference operator
	  T& operator* () const
	  {
		return *p;
	  }

      //! check whether both point to same target
	  bool operator== (const CP<T>& cp) const
	  {
		return p==cp.p;
	  }

      //! check whether target have different adress
	  bool operator!= (const CP<T>& cp) const
	  {
		return p!=cp.p;
	  }

	};

	class CountableException
	{
	  int counter;
	public:
	  CountableException (int i) : counter(i) {}
	  int get_counter () const
	  {
		return counter;
	  }
	};

    /** @brief Base class for object pointed to by CP
     *
     *  This provides the necessary functionality in the target object for the
     *  CP template class to work.
     */
	class Countable 
	{
	  mutable int counter;

	public:

      //! Default constructor
	  Countable () : counter(0) 
	  {
	  }

	  //! copy constructor: new object, no pointer exists
	  Countable (const Countable& )
	  {
		counter = 0;
	  } 

	  //! number of pointers does not change
	  Countable& operator= (const Countable& )
	  {
        return *this;
	  } 

      //! increment reference counter
	  void reference_counter_increment () const
	  {
		counter++;
	  }

      //! decrement reference counter
	  void reference_counter_decrement () const
	  {
		counter--;
	  }

      //! check wether the reference counter is zero
	  bool reference_counter_zero () const
	  {
		return counter==0;
	  }

      //! get value of reference counter
	  int get_reference_counter () const
	  {
		return counter;
	  }

      /** @brief Destructor
       *
       *  Warn if any CP is still pointing to us.
       */
	  ~Countable ()
	  {
		if (counter!=0)
		  {
			std::cout << counter << " counting pointer(s) point to object at "
					  << this << " while it is deleted" << std::endl;
            throw CountableException(counter);
		  }
	  }
	};

} // namespace hdnum

#endif
