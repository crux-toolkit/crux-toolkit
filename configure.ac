# automake/autoconf startup stuff
AC_PREREQ(2.57)
AC_INIT(crux, 1.0, aklammer@u.washington.edu)
AC_CONFIG_SRCDIR(src/c/ion_series.cpp)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(crux, 1.31)
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

# Use double precision floating point through out.
doubles=false
AC_ARG_ENABLE(
    doubles,
    [  --enable-doubles  Use double precision floating point througout. Default: disabled],
    [doubles="true"])
AM_CONDITIONAL(USEDOUBLES, test "x$doubles" = xtrue )

# The search-for-matches functionality depends on SEQUEST intellectual property.
# By default we will include this in the build and distribution, but we allow 
# for disabling it, and not including it in the distribution.
search="true"
AC_ARG_ENABLE(
    search,
    [  --disable-search  Disable search-for-matches functionality. Default: enabled],
    [search="false"])
AM_CONDITIONAL(SEARCH, test "x$search" = xtrue )

# These variables will be available in the Makefiles.
AC_SUBST(MACHID)
AC_SUBST(VERSION)

# List all the Makefiles for the project.
AC_CONFIG_FILES([\
    Makefile \
    doc/Makefile \
    doc/user/Makefile \
    src/Makefile \
    src/c/Makefile \
    src/c/ArrayLibrary/Makefile \
    src/c/MSToolkit/Makefile \
    src/c/percolator/Makefile \
    src/c/q-ranker/Makefile \
    src/c/xlink/Makefile \
    src/python/Makefile \
])

AM_CONDITIONAL(DARWIN, \
  test [ "x$MACHID" = "xdarwin" ] \
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
  System:                     ${MACHID}
  Use doubles:                ${doubles}

  Run the following commands to compile and install crux:

        make   
        make install

  Run the following command to uninstall crux:

        make uninstall
"
if test "$search" = "false"
then
echo "Peptide-spectrum database search is disabled."
fi
