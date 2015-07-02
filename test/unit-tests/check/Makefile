CC = gcc
CFLAGS=-I.
CHECK_OBJ = \
  check.o \
  check_error.o \
  check_list.o \
  check_log.o \
  check_msg.o \
  check_pack.o \
  check_print.o \
  check_run.o \
  check_str.o

check: $(CHECK_OBJ)
	ar rcs check.ar $(CHECK_OBJ)


clean:
	rm -f *.o *.ar

