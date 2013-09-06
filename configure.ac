#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([overly_ambitious_isqo], [0.1], [traviscj@traviscj.com])
AM_INIT_AUTOMAKE(foreign -Wall -Werror)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])


# Checks for programs.
AC_LANG([C++])
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AC_PROG_MAKE_SET



# Checks for libraries.
AC_CHECK_LIB([blas], [dgemm_],[],AC_MSG_ERROR([Could not find blas!!]))
AC_CHECK_LIB([gfortran], [_gfortran_st_write],[],AC_MSG_ERROR([Could not find -lgfortran!!]))
AC_CHECK_LIB([lapack], [dtrtrs_],[],AC_MSG_ERROR([Could not find -llapack!!]))

# adding the lib to the files to link
LIBS="-lqpOASES_aw $LIBS"
# introduce the optional configure parameter for a non-standard install prefix of qpOASES_aw
AC_ARG_WITH([qpOASES_aw],
    [AS_HELP_STRING([--with-qpOASES_aw=prefix],
        [try this for a non-standard install prefix of the qpOASES_aw library])],
    [qpOASES_awPATHSET=1],
    [qpOASES_awPATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $qpOASES_awPATHSET = 1 ; then

    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_qpOASES_aw/include_aw"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_qpOASES_aw/bin"
fi

##########################################################################
# distribute the changed variables among the Makefiles
# AC_SUBST([LIBS])
# AC_SUBST([AM_CXXFLAGS])
# AC_SUBST([AM_LDFLAGS])

# to check qpOASES, need to pass the library flags.
AX_CXX_CHECK_LIB([qpOASES_aw], [qpOASES::SparseMatrix], [], AC_MSG_ERROR([Could not find qpOASES_aw!!]), $AM_LDFLAGS)


##########################################################################
# adding the amplsolver library (e.g. with static name 'amplsolver.a')
##########################################################################

AC_CHECK_HEADER([asl.h],
[],
[AC_MSG_ERROR([required header file "asl.h" not found. (specify with --with-amplsolver= ?)])])


# introduce the optional configure parameter for a non-standard install prefix of amplsolver
AC_ARG_WITH([amplsolver],
    [AS_HELP_STRING([--with-amplsolver=prefix],
        [try this for a non-standard install prefix of the qpOASES_aw library])],
    [amplsolverPATHSET=1],
    [amplsolverPATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $amplsolverPATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_amplsolver/"
    ASL_LOCATION="$with_amplsolver"
fi

AC_CHECK_FILE([$ASL_LOCATION/amplsolver.a], [], AC_MSG_ERROR([Could not find amplsolver.a in "$ASL_LOCATION"; 
check --with-amplsolver option!!!]))


AC_SUBST([LIBS])
AC_SUBST([ASL_LOCATION])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
