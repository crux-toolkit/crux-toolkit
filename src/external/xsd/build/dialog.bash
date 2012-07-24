# file      : build/dialog.bash
# author    : Boris Kolpackov <boris@codesynthesis.com>
# copyright : Copyright (c) 2004-2011 Code Synthesis Tools CC
# license   : GNU GPL v2; see accompanying LICENSE file

# bld_root  - build root
# MAKE      - GNU make executable
#

#@@ can't I do exec 1>&2 or something?
#
function echo_error ()
{
  echo "$*" 1>&2
}

echo=echo_error

function abspath_old ()
{
  local r=`eval echo $1` # for tilde-expansion

  if [ "`echo $r | sed -e 's%^/.*$%%'`" ]; then

    r=`pwd`/$r
  fi

  # remove trailing `/'
  #
  r=`echo $r | sed -e 's%^\(.\+\)/$%\1%'`

  echo $r
}

# Get rid of jobserver info. See GNU make bug #12229.
#
MAKEFLAGS=`echo $MAKEFLAGS | sed -e 's/--jobserver-fds=[^ ]* *//g'`
MAKEFLAGS=`echo $MAKEFLAGS | sed -e 's/-j *//g'`


function abspath ()
{
  local r=`eval echo $1` # for tilde-expansion

  $MAKE --no-print-directory -f $bld_root/abspath.make $r
}

# $1  Default answer ('y' or 'n'). [optional]
#
function read_y_n ()
{
  local r

  while [ -z "$r" ]; do

    read -p "[$1]: " r

    if [ -z "$r" ]; then
      r=$1
    fi

    case "$r" in
      y) ;;
      n) ;;
      *) r= ;;
    esac
  done

  echo $r
}

# $1  Default answer. [optional]
#
function read_value ()
{
  local r

  read -e -p "[$1]: " r

  if [ -z "$r" ]; then
    r=$1
  fi

  echo $r
}

# $1  Space-separated list of options.
# $2  Default answer from the list of options. [optional]
#
function read_option ()
{
  local d o r i

  if [ "$2" ]; then

    d=1

    for o in $1; do

      if [ "$o" = "$2" ]; then break; fi

      d=$(($d + 1))
    done
  fi

  while [ -z "$r" ]; do

    read -p "[$d]: " r

    if [ -z "$r" ]; then

      r=$d
    fi

    i=1

    for o in $1; do

      if [ "$i" = "$r" ]; then

	echo $o
	return 0
      fi

      i=$(($i + 1))
    done

    # User must have entered some junk so let's ask her again.
    #
    r=
  done
}


# $1  Default answer. [optional]
#
# --command   - path to be read is a command (implies --exist)
#
# --directory - path to be read is a directory
#
# --exist     - path to be read should exist
#
#
function read_path ()
{
  local command=0
  local directory=0
  local exist=0

  # Parse options
  #
  while [ "$1" ]; do
    case $1 in

      --command)

	command=1
	exist=1
	shift
	;;

      --directory)

	directory=1
	shift
	;;

      --exist)

	exist=1
	shift
	;;

      *)

	break
	;;

    esac
  done

  local r tmp

  while [ -z "$r" ]; do

    read -e -p "[$1]: " r

    if [ -z "$r" ]; then

      r=$1
    fi

    if [ -z "$r" ]; then

      continue
    fi


    if [ $command == 1 ]; then

      # Do tilde expansion.
      #
      r=`eval echo $r`

      # Extract first word.
      #
      local cmd=`echo $r | cut -d " " -f 1`

      type -p $cmd 1>/dev/null 2>&1

      if [ $? != 0 ]; then

	$echo "$r: command not found"
	r=

      fi

      continue
    fi

    # abspath & checks
    #
    r=`abspath $r`

    if [ $exist == 1 ]; then

      if [ $directory == 1 ]; then

	if [ ! -d "$r" ]; then

	  $echo "$r: no such directory"
	  r=
	  continue
	fi

      else # something else

	if [ ! -e "$r" ]; then

	  $echo "$r: no such file or directory"
	  r=
	  continue
	fi
      fi
    fi
  done

  echo $r
}
