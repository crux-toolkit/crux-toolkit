// Benjamin Diament
#include <gflags/gflags.h>
#include "max_mz.h"

DEFINE_double(max_bin, 0, "During search, ignore peaks with higher m/z than "
	      "this. Zero means use all peaks in the dataset. During indexing,"
	      " 0 means to include all theoretical peaks.");

MaxBin MaxBin::global;
