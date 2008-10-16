# automake/autoconf startup stuff
AC_PREREQ(2.57)
AC_INIT(crux, 1.0, aklammer@u.washington.edu)
AC_CONFIG_SRCDIR(src/c/ion_series.c)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(crux, 1.01)
AM_CONFIG_HEADER([config.h])

# Check for required programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(PYTHON, python)

# Checks for required header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([ \
  fcntl.h \
  inttypes.h \
  limits.h \
  malloc.h netinet/in.h \
  stddef.h \
  stdint.h \
  stdlib.h \
  string.h \
  sys/time.h \
  unistd.h \
])

# Check for required library functions.
AC_FUNC_CLOSEDIR_VOID
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([ \
  getpagesize \
  gettimeofday \
  memmove \
  memset \
  mkdir \
  munmap \
  pow \
  rmdir \
  sqrt \
  strcasecmp \
  strchr \
  strdup \
  strrchr \
  strstr \
])


# Checks for required math functions.
AC_CHECK_LIB(
  [m], 
  [sqrt, log, exp, fabs, pow],
  [], 
  [AC_MSG_WARN([Math library not found])]
)

# Set platform-specific settings
case $host in
   *cygwin*)
      AC_DEFINE(Cygwin, [], [Using cygwin])
      MACHID="cygwin"
      ;;
   *darwin*)
      AC_DEFINE(MACOSX, [], [Using MacOS X])
      MACHID="darwin"
      ;;
   *linux*)
      AC_DEFINE(Linux, [], [Using Linux])
      MACHID="linux"
      ;;
esac

####################################
# Check for configuration switches
####################################

# set percolator location
AC_ARG_WITH(
    percolator-dir,
    [  --with-percolator-dir=DIR  Set path to percolator source. \
      Default: ${top_srcdir}/percolator],
    [PERCOLATOR_DIR=$withval],
    [PERCOLATOR_DIR="\${top_srcdir}/../percolator"])

# These variables will be available in the Makefiles.
AC_SUBST(MACHID)
AC_SUBST(VERSION)
AC_SUBST(PERCOLATOR_DIR)

# List all the Makefiles for the project.
AC_CONFIG_FILES([\
    Makefile \
    doc/Makefile \
    doc/user/Makefile \
    src/Makefile \
    src/c/Makefile \
])

PERCOLATOR_HEADER_STATUS=""
PERCOLATOR_LIB_STATUS=""
if test ! -e "${PERCOLATOR_DIR}/PercolatorCInterface.h"
then
  PERCOLATOR_HEADER_STATUS="$PERCOLATOR_DIR/PercolatorCInterface.h not found!"
fi
if test ! -e "${PERCOLATOR_DIR}/libpercolator.la"
then
  PERCOLATOR_LIB_STATUS="$PERCOLATOR_DIR/libpercolator.la not found!"
fi

AM_CONDITIONAL(PERCOLATOR, \
  test [ -z "$PERCOLATOR_HEADER_STATUS" -a -z "$PERCOLATOR_LIB_STATUS" ] \
)

AC_OUTPUT

echo \
"
        ========================
        Configuration parameters
        ========================

  Install path:               ${prefix}
  Version:                    ${VERSION}
  C compiler:                 ${CC}
  C compiler flags:           ${CFLAGS}
  percolator directory:       ${PERCOLATOR_DIR} ${PERCOLATOR}

  Run the following commands to compile and install crux:

        make   
        make install

  Run the following command to uninstall crux:

        make uninstall
"

if test -n "${PERCOLATOR_HEADER_STATUS}"
then
  echo
  echo "Warning: $PERCOLATOR_HEADER_STATUS"
fi
if test -n "${PERCOLATOR_LIB_STATUS}"
then
  echo
  echo "Warning: $PERCOLATOR_LIB_STATUS"
fi
if test -n "${PERCOLATOR_LIB_STATUS}" -o -n "$PERCOLATOR_HEADER_STATUS"
then
  echo
  echo "Have you downloaded and built Percolator?"
  echo 
  echo "Percolator support disabled!"
fi
