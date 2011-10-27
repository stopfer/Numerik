#ifndef HDNUM_EXCEPTIONS_HH
#define HDNUM_EXCEPTIONS_HH

#include <string>
#include <sstream>

namespace hdnum {

  /*! \file
	\brief A few common exception classes
  
	This file defines a common framework for generating exception
	subclasses and to throw them in a simple manner.
	Taken from the DUNE project www.dune-project.org

  */


  /*! \brief Base class for Exceptions
  
	all HDNUM exceptions are derived from this class via trivial subclassing:

	\code
	class MyException : public Dune::Exception {};
	\endcode

	You should not \c throw a Dune::Exception directly but use the macro
	DUNE_THROW() instead which fills the message-buffer of the exception
	in a standard way and features a way to pass the result in the
	operator<<-style

	\see HDNUM_THROW, IOError, MathError

  */
  class Exception {
  public:
	void message(const std::string &message); //!< store string in internal message buffer
	const std::string& what() const;          //!< output internal message buffer
  private:
	std::string _message;
  };

  inline void Exception::message(const std::string &message)
  {
	_message = message;
  }

  inline const std::string& Exception::what() const
  {
	return _message;
  }

  inline std::ostream& operator<<(std::ostream &stream, const Exception &e)
  {
	return stream << e.what();
  }

  // the "format" the exception-type gets printed.  __FILE__ and
  // __LINE__ are standard C-defines, the GNU cpp-infofile claims that
  // C99 defines __func__ as well. __FUNCTION__ is a GNU-extension
#ifdef HDNUM_DEVEL_MODE
# define THROWSPEC(E) #E << " [" << __func__ << ":" << __FILE__ << ":" << __LINE__ << "]: "
#else
# define THROWSPEC(E) #E << ": "
#endif

  /*! Macro to throw an exception

	\param E exception class derived from Dune::Exception
	\param m reason for this exception in ostream-notation

	Example:

	\code
	if (filehandle == 0)
    DUNE_THROW(FileError, "Could not open " << filename << " for reading!")
	\endcode

	DUNE_THROW automatically adds information about the exception thrown
	to the text. If DUNE_DEVEL_MODE is defined more detail about the
	function where the exception happened is included. This mode can be
	activated via the \c --enable-dunedevel switch of \c ./configure
  */
  // this is the magic: use the usual do { ... } while (0) trick, create
  // the full message via a string stream and throw the created object
#define HDNUM_THROW(E, m) do { E th__ex; std::ostringstream th__out;		\
	th__out << THROWSPEC(E) << m; th__ex.message(th__out.str()); throw th__ex; \
  } while (0)

  /*! \brief Default exception class for I/O errors

	This is a superclass for any errors dealing with file/socket I/O problems
	like

	- file not found
	- could not write file
	- could not connect to remote socket
  */
  class IOError : public Exception {};

  /*! \brief Default exception class for mathematical errors

	This is the superclass for all errors which are caused by
	mathematical problems like

	- matrix not invertible
	- not convergent
  */
  class MathError : public Exception {};

  /*! \brief Default exception class for range errors

	This is the superclass for all errors which are caused because
	the user tries to access data that was not allocated before.
	These can be problems like

	- accessing array entries behind the last entry
	- adding the fourth non zero entry in a sparse matrix
	with only three non zero entries per row
  
  */
  class RangeError : public Exception {};

  /*! \brief Default exception for dummy implementations 

	This exception can be used for functions/methods

	- that have to be implemented but should never be called
	- that are missing
  */
  class NotImplemented : public Exception {};

  /*! \brief Default exception class for OS errors

	This class is thrown when a system-call is used and returns an
	error.
  
  */
  class SystemError : public Exception {};

  /*! \brief Default exception if memory allocation fails

   */
  class OutOfMemoryError : public SystemError {};

  /*! \brief Default exception if a function was called while
   * the object is not in a valid state for that function.
   */
  class InvalidStateException : public Exception {};

  /*! \brief General Error
   */
  class ErrorException : public Exception {};
  
  // throw ErrorException with message
#define HDNUM_ERROR(m) do { hdnum::ErrorException th__ex; std::ostringstream th__out;		\
	th__out << THROWSPEC(hdnum::ErrorException) << m; \
	th__ex.message(th__out.str()); \
	std::cout << th__ex.what() << std::endl; \
	throw th__ex;				   \
  } while (0)

} // end namespace

#endif
