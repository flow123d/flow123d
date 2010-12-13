/* base/include/base/config.h.  Generated from config.h.in by configure.  */
/* base/include/base/config.h.in.  Generated from configure.in by autoheader.  */


//----------------------------  config.h  ---------------------------
//    $Id: config.h.in 18724 2009-04-23 23:06:37Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  config.h  ---------------------------
#ifndef __deal2__config_h
#define __deal2__config_h

// Note: you should not usually have to change this file, as it is
// automatically generated from the ./configure file in the top level
// directory. If there are problems with the contents of this file,
// rather than changing it, try to modify the mechanisms in
// configure.in that generated this output. The reason is that you
// would have to make these changes each time you compile a new
// version of the library, or on a different computer. Furthermore, it
// is important not to build different parts of the library with
// different versions of this file.
//
// In case of problems in autodetection of features of your build
// environment, contact the authors of the library.


/**
 * Two macro names that we put at the top and bottom of all deal.II files
 * and that will be expanded to "namespace dealii {" and "}".
 */
#define DEAL_II_NAMESPACE_OPEN namespace dealii {
#define DEAL_II_NAMESPACE_CLOSE }



/* Defined if the prototype of abort() has a no-throw exception specification.
   */
#define DEAL_II_ABORT_NOTHROW_EXCEPTION 1

/* Flag indicating whether there is a bug in the compiler that leads to bogus
   warnings for inline class members in anonymous namespaces */
/* #undef DEAL_II_ANON_NAMESPACE_BOGUS_WARNING */

/* Defined if the compiler needs to see the static keyword even for functions
   in anonymous namespaces, to avoid duplicate symbol errors when linking. For
   the details, look at aclocal.m4 in the top-level directory. */
#define DEAL_II_ANON_NAMESPACE_BUG 1

/* Another test if the compiler needs to see the static keyword even for
   functions in anonymous namespaces, to avoid duplicate symbol errors when
   linking. For the details, look at aclocal.m4 in the top-level directory. */
#define DEAL_II_ANON_NAMESPACE_LINKAGE_BUG 1

/* Defined if the compiler has a problem with assigning arrays in conditionals
   */
/* #undef DEAL_II_ARRAY_CONDITIONAL_DECAY_BUG */

/* Define if the use of socket functionality leads to strange results with
   floating point computations on cygwin systems. */
/* #undef DEAL_II_BROKEN_SOCKETS */

/* Defined if the compiler we use supports the upcoming C++1x standard. */
/* #undef DEAL_II_CAN_USE_CXX1X */

/* Backward compatibility support for functions and classes that do not take
   an explicit mapping variable, but rather use a default Q1 mapping instead
   */
#define DEAL_II_COMPAT_MAPPING true

/* Defined if the compiler supports including <mpi.h> */
#define DEAL_II_COMPILER_SUPPORTS_MPI 1

/* disable the function parser in contrib */
/* #undef DEAL_II_DISABLE_PARSER */

/* Defined if the compiler does not honor the explicit keyword on template
   constructors. */
/* #undef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG */

/* Defined if the compiler needs a workaround for certain problems with taking
   the address of template template functions. For the details, look at
   aclocal.m4 in the top-level directory. */
/* #undef DEAL_II_FUNPTR_TEMPLATE_TEMPLATE_BUG */

/* Defined if std::isfinite is available */
#define DEAL_II_HAVE_ISFINITE 1

/* Flag indicating whether the library shall be compiled to use the Tecplot
   interface */
/* #undef DEAL_II_HAVE_TECPLOT */

/* Defined if the compiler refuses to compile the definition of a function
   that was previously declared abstract. */
#define DEAL_II_IMPLEMENTED_PURE_FUNCTION_BUG 1

/* Define if we have to work around a bug in Sun's Forte compiler. See the
   aclocal.m4 file in the top-level directory for a description of this bug.
   */
/* #undef DEAL_II_LOCAL_TYPEDEF_COMP_WORKAROUND */

/* Defined if the compiler gets an internal compiler upon some code involving
   long doubles, and with optimization. For the details, look at aclocal.m4 in
   the top-level directory. */
/* #undef DEAL_II_LONG_DOUBLE_LOOP_BUG */

/* Major version number of deal.II */
#define DEAL_II_MAJOR 6

/* Defined if the compiler refuses to allow the explicit specialization of
   static member arrays. For the exact failure mode, look at aclocal.m4 in the
   top-level directory. */
/* #undef DEAL_II_MEMBER_ARRAY_SPECIALIZATION_BUG */

/* Define if we have to work around a bug in gcc with explicitly instantiating
   template member operators. See the aclocal.m4 file in the top-level
   directory for a description of this bug. */
#define DEAL_II_MEMBER_OP_TEMPLATE_INST 

/* Defined if the compiler refuses to specialize an outer class template while
   keeping a member as a template. For the exact failure mode, look at
   aclocal.m4 in the top-level directory. */
/* #undef DEAL_II_MEMBER_TEMPLATE_SPECIALIZATION_BUG */

/* Defined if the compiler refuses to allow the explicit specialization of
   static member variables. For the exact failure mode, look at aclocal.m4 in
   the top-level directory. */
/* #undef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG */

/* Minor version number of deal.II */
#define DEAL_II_MINOR 2

/* Set to the minimal number of elements a std::vector<bool> can always hold,
   i.e. its minimal capacity. */
#define DEAL_II_MIN_BOOL_VECTOR_CAPACITY 64

/* Set to the minimal number of elements a std::vector<T> can always hold,
   i.e. its minimal capacity. */
#define DEAL_II_MIN_VECTOR_CAPACITY 1

/* Define if we have to work around a bug in gcc with marking all instances of
   a template class as friends to this class if the class is inside a
   namespace. See the aclocal.m4 file in the top-level directory for a
   description of this bug. */
/* #undef DEAL_II_NAMESP_TEMPL_FRIEND_BUG */

/* Define if we have to work around another bug in gcc with marking all
   instances of a template class as friends to this class if the class is
   inside a namespace. See the aclocal.m4 file in the top-level directory for
   a description of this bug. */
/* #undef DEAL_II_NAMESP_TEMPL_FRIEND_BUG2 */

/* Defined if the compiler does not properly implement the resolution of
   defect report #45 to the C++ standard, which makes nested types implicit
   friends of the enclosing class. */
/* #undef DEAL_II_NESTED_CLASS_FRIEND_BUG */

/* Defined if the compiler does not understand friend declarations for nested
   member classes when giving a full class specification. */
/* #undef DEAL_II_NESTED_CLASS_TEMPL_FRIEND_BUG */

/* Path to the deal.II directory */
#define DEAL_II_PATH "/home/jb/local/deal.II-6.2.1"

/* Defined if the compiler does not support the
   substitution-failure-is-not-an-error paradigm. For the details, look at
   aclocal.m4 in the top-level directory. */
/* #undef DEAL_II_SFINAE_BUG */

/* Define if we have to work around a bug in Sun's Forte compiler. See the
   aclocal.m4 file in the top-level directory for a description of this bug.
   */
/* #undef DEAL_II_TEMPLATE_SPEC_ACCESS_WORKAROUND */

/* Defined if the compiler refuses to allow a typedef to a template template
   class template parameter. For the exact failure mode, look at aclocal.m4 in
   the top-level directory. */
/* #undef DEAL_II_TEMPLATE_TEMPLATE_TYPEDEF_BUG */

/* Defined if the compiler requires the use of the template keyword for
   disambiguation keyword in certain contexts in which it is not supposed to
   do so. For the exact failure mode, look at aclocal.m4 in the top-level
   directory. */
/* #undef DEAL_II_TEMPL_OP_DISAMBIGUATION_BUG */

/* Define if we have to work around a bug with some compilers that will not
   allow us to specify a fully specialized class of a template as a friend.
   See the aclocal.m4 file in the top-level directory for a description of
   this bug. */
/* #undef DEAL_II_TEMPL_SPEC_FRIEND_BUG */

/* Define if the compiler provides a <errno.g> header file which does not
   define all error codes such as EINTR. In that case, use the system include
   file at /usr/include instead. There is probably a better way to do this,
   but it is not apparent by looking at the C/C++ compatibility header
   provided by the compiler. */
/* #undef DEAL_II_USE_DIRECT_ERRNO_H */

/* Defined if a Metis installation was found and is going to be used */
/* #undef DEAL_II_USE_METIS */

/* Flag indicating whether the library shall be compiled for multithreaded
   applications. If so, then it is set to one, otherwise to zero. */
#define DEAL_II_USE_MT 0

/* Defined if multi-threading is to be achieved by using the POSIX functions
   */
/* #undef DEAL_II_USE_MT_POSIX */

/* Defined if POSIX is supported but not the newer POSIX barrier functions.
   Barriers will then not work in the library, but the other threading
   functionality is available. */
/* #undef DEAL_II_USE_MT_POSIX_NO_BARRIERS */

/* Defined if a PETSc installation was found and is going to be used */
#define DEAL_II_USE_PETSC 1

/* Defined if a SLEPc installation was found and is going to be used */
/* #undef DEAL_II_USE_SLEPC */

/* Defined if a Trilinos installation was found and is going to be used */
/* #undef DEAL_II_USE_TRILINOS */

/* This error appears in the Apple edition of the gcc 3.3, which ships with
   Darwin7.9.0 and probably previous version. It leads to problems during
   linking. For the details, look at aclocal.m4 in the top-level directory. */
/* #undef DEAL_II_WEAK_LINKAGE_BUG */

/* Define to 1 if you have the <Amesos.h> header file. */
/* #undef HAVE_AMESOS_H */

/* Define to 1 if you have the <AztecOO.h> header file. */
/* #undef HAVE_AZTECOO_H */

/* Define to 1 if you have the <AztecOO_Operator.h> header file. */
/* #undef HAVE_AZTECOO_OPERATOR_H */

/* Define if the compiler provides __builtin_expect */
#define HAVE_BUILTIN_EXPECT 1

/* Define to 1 if you have the `daxpy_' function. */
#define HAVE_DAXPY_ 1

/* Define to 1 if you have the `dgeevx_' function. */
#define HAVE_DGEEVX_ 1

/* Define to 1 if you have the `dgeev_' function. */
#define HAVE_DGEEV_ 1

/* Define to 1 if you have the `dgemm_' function. */
#define HAVE_DGEMM_ 1

/* Define to 1 if you have the `dgemv_' function. */
#define HAVE_DGEMV_ 1

/* Define to 1 if you have the `dgesvd_' function. */
#define HAVE_DGESVD_ 1

/* Define to 1 if you have the `dgetrf_' function. */
#define HAVE_DGETRF_ 1

/* Define to 1 if you have the `dgetri_' function. */
#define HAVE_DGETRI_ 1

/* Define to 1 if you have the `dgetrs_' function. */
#define HAVE_DGETRS_ 1

/* Define to 1 if you have the `dstev_' function. */
#define HAVE_DSTEV_ 1

/* Define to 1 if you have the <Epetra_CrsGraph.h> header file. */
/* #undef HAVE_EPETRA_CRSGRAPH_H */

/* Define to 1 if you have the <Epetra_CrsMatrix.h> header file. */
/* #undef HAVE_EPETRA_CRSMATRIX_H */

/* Define to 1 if you have the <Epetra_Import.h> header file. */
/* #undef HAVE_EPETRA_IMPORT_H */

/* Define to 1 if you have the <Epetra_LinearProblem.h> header file. */
/* #undef HAVE_EPETRA_LINEARPROBLEM_H */

/* Define to 1 if you have the <Epetra_Map.h> header file. */
/* #undef HAVE_EPETRA_MAP_H */

/* Define to 1 if you have the <Epetra_MultiVector.h> header file. */
/* #undef HAVE_EPETRA_MULTIVECTOR_H */

/* Define to 1 if you have the <Epetra_Operator.h> header file. */
/* #undef HAVE_EPETRA_OPERATOR_H */

/* Define to 1 if you have the <Epetra_SerialComm.h> header file. */
/* #undef HAVE_EPETRA_SERIALCOMM_H */

/* Define to 1 if you have the <Epetra_Vector.h> header file. */
/* #undef HAVE_EPETRA_VECTOR_H */

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to 1 if you have the `getpid' function. */
#define HAVE_GETPID 1

/* Define if deal.II is linked against a libc that provides stacktrace debug
   information that can be printed out in the exception class */
#define HAVE_GLIBC_STACKTRACE 1

/* Availability of the MA27 algorithm from HSL */
/* #undef HAVE_HSL_MA27 */

/* Availability of the MA47 algorithm from HSL */
/* #undef HAVE_HSL_MA47 */

/* Define to 1 if you have the <Ifpack.h> header file. */
/* #undef HAVE_IFPACK_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Defined if deal.II was configured with BLAS support */
#define HAVE_LIBBLAS 1

/* Defined if deal.II was configured with LAPACK support */
#define HAVE_LIBLAPACK 1

/* Define to 1 if you have the `NetCDF' library (-lnetcdf). */
/* #undef HAVE_LIBNETCDF */

/* Define if the std c++ library provides a demangler conforming to the GCC
   libstdc++ interface. */
#define HAVE_LIBSTDCXX_DEMANGLER 1

/* UMFPACK is */
#define HAVE_LIBUMFPACK 1

/* Define to 1 if you have the `z' library (-lz). */
/* #undef HAVE_LIBZ */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <ml_MultiLevelPreconditioner.h> header file. */
/* #undef HAVE_ML_MULTILEVELPRECONDITIONER_H */

/* Define if you have the rand_r function */
/* #undef HAVE_RAND_R */

/* Define to 1 if you have the <Sacado.hpp> header file. */
/* #undef HAVE_SACADO_HPP */

/* Define to 1 if you have the `saxpy_' function. */
#define HAVE_SAXPY_ 1

/* Define to 1 if you have the `sgeevx_' function. */
#define HAVE_SGEEVX_ 1

/* Define to 1 if you have the `sgeev_' function. */
#define HAVE_SGEEV_ 1

/* Define to 1 if you have the `sgemm_' function. */
#define HAVE_SGEMM_ 1

/* Define to 1 if you have the `sgemv_' function. */
#define HAVE_SGEMV_ 1

/* Define to 1 if you have the `sgesvd_' function. */
#define HAVE_SGESVD_ 1

/* Define to 1 if you have the `sgetrf_' function. */
#define HAVE_SGETRF_ 1

/* Define to 1 if you have the `sgetri_' function. */
#define HAVE_SGETRI_ 1

/* Define to 1 if you have the `sgetrs_' function. */
#define HAVE_SGETRS_ 1

/* Define to 1 if you have the `sstev_' function. */
#define HAVE_SSTEV_ 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define if the compiler provides an <iosfwd> header file */
#define HAVE_STD_IOSFWD_HEADER 1

/* Define if the compiler's library in use provides a std::iterator class
   (early gcc versions did not) */
#define HAVE_STD_ITERATOR_CLASS 1

/* Define if the compiler's library in use provides std::numeric_limits
   classes in the appropriate header file */
#define HAVE_STD_NUMERIC_LIMITS 1

/* Define if the compiler provides an <ostream> header file */
#define HAVE_STD_OSTREAM_HEADER 1

/* Define if the compiler's library in use provides std::i/ostringstream
   classes (early gcc versions did not) */
#define HAVE_STD_STRINGSTREAM 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/syscall.h> header file. */
#define HAVE_SYS_SYSCALL_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <Teuchos_ParameterList.hpp> header file. */
/* #undef HAVE_TEUCHOS_PARAMETERLIST_HPP */

/* Define to 1 if you have the <Teuchos_RCP.hpp> header file. */
/* #undef HAVE_TEUCHOS_RCP_HPP */

/* Define to 1 if you have the <Teuchos_RefCountPtr.hpp> header file. */
/* #undef HAVE_TEUCHOS_REFCOUNTPTR_HPP */

/* Define to 1 if you have the <Thyra_AztecOOLinearOpWithSolveFactory.hpp>
   header file. */
/* #undef HAVE_THYRA_AZTECOOLINEAROPWITHSOLVEFACTORY_HPP */

/* Define to 1 if you have the <Thyra_DefaultBlockedLinearOpDecl.hpp> header
   file. */
/* #undef HAVE_THYRA_DEFAULTBLOCKEDLINEAROPDECL_HPP */

/* Define to 1 if you have the <Thyra_DefaultBlockedLinearOp.hpp> header file.
   */
/* #undef HAVE_THYRA_DEFAULTBLOCKEDLINEAROP_HPP */

/* Define to 1 if you have the <Thyra_DefaultInverseLinearOp.hpp> header file.
   */
/* #undef HAVE_THYRA_DEFAULTINVERSELINEAROP_HPP */

/* Define to 1 if you have the <Thyra_EpetraLinearOp.hpp> header file. */
/* #undef HAVE_THYRA_EPETRALINEAROP_HPP */

/* Define to 1 if you have the <Thyra_EpetraThyraWrappers.hpp> header file. */
/* #undef HAVE_THYRA_EPETRATHYRAWRAPPERS_HPP */

/* Define to 1 if you have the <Thyra_InverseLinearOperator.hpp> header file.
   */
/* #undef HAVE_THYRA_INVERSELINEAROPERATOR_HPP */

/* Define to 1 if you have the <Thyra_LinearOperatorDecl.hpp> header file. */
/* #undef HAVE_THYRA_LINEAROPERATORDECL_HPP */

/* Define to 1 if you have the <Thyra_LinearOperatorImpl.hpp> header file. */
/* #undef HAVE_THYRA_LINEAROPERATORIMPL_HPP */

/* Define to 1 if you have the <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
   header file. */
/* #undef HAVE_THYRA_LINEAROPWITHSOLVEFACTORYHELPERS_HPP */

/* Define to 1 if you have the <Thyra_MultiVectorBase.hpp> header file. */
/* #undef HAVE_THYRA_MULTIVECTORBASE_HPP */

/* Define to 1 if you have the <Thyra_MultiVectorDefaultBase.hpp> header file.
   */
/* #undef HAVE_THYRA_MULTIVECTORDEFAULTBASE_HPP */

/* Define to 1 if you have the <Thyra_VectorDecl.hpp> header file. */
/* #undef HAVE_THYRA_VECTORDECL_HPP */

/* Define to 1 if you have the <Thyra_VectorImpl.hpp> header file. */
/* #undef HAVE_THYRA_VECTORIMPL_HPP */

/* Define to 1 if you have the <Thyra_VectorSpaceImpl.hpp> header file. */
/* #undef HAVE_THYRA_VECTORSPACEIMPL_HPP */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if the compiler provides __verbose_terminate_handler */
#define HAVE_VERBOSE_TERMINATE 1

/* On SunOS 4.x, the getrusage() function exists, but is not declared in the
   respective header file <resource.h>, as one would think when reading the
   man pages. Then we have to declare this function ourselves in those files
   that use this function. The question whether we have to do so is controlled
   by the preprocessor variable. */
/* #undef NO_HAVE_GETRUSAGE */

/* Define to the address where bug reports for this package should be sent. */
#define DEAL_II_PACKAGE_BUGREPORT "dealii@dealii.org"

/* Define to the full name of this package. */
#define DEAL_II_PACKAGE_NAME "deal.II"

/* Define to the full name and version of this package. */
#define DEAL_II_PACKAGE_STRING "deal.II 6.2.1"

/* Define to the one symbol short name of this package. */
#define DEAL_II_PACKAGE_TARNAME "deal.II"

/* Define to the version of this package. */
#define DEAL_II_PACKAGE_VERSION "6.2.1"

/* Make sure PETSc doesn't re-define the underscore through the preprocessor,
   since this interferes with boost. PETSc redefines the underscore to be
   "__gterr =", but then forgets to undef this thing. Boost simply wants to
   concatenate the underscore with another string to form a class name, which
   then of course isn't valid any more. See mails in early Feb 2006. */
#define PETSC_SKIP_UNDERSCORE_CHKERR 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* If already available, do not define at all. Otherwise, define to __func__
   if that is available. In all other cases, indicate that no information
   about the present function is available for this compiler. */
/* #undef __PRETTY_FUNCTION__ */


/**
 * Depending on the use of threads, we will have to make some variables
 * volatile. We do this here in a very old-fashioned C-style, but still
 * convenient way.
 */
#if DEAL_II_USE_MT != 0
#  define DEAL_VOLATILE volatile
#else
#  define DEAL_VOLATILE
#endif

#include <numbers.h>

/**
 * If the compiler supports the upcoming C++1x standard, allow us to refer
 * to things in namespace std through namespace std_cxx1x (the namespace
 * into which we import BOOST components if we don't have C++1x support).
 */
#ifdef DEAL_II_CAN_USE_CXX1X
namespace std_cxx1x = std;
#endif

#endif

