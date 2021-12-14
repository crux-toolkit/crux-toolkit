import argparse
import sys
import os

PROG_NAME = "Python file comparator"
PRECISION_VAR = "COMPARE_DIGITS"

def assertAndDie(condition, message):
    if not condition:
        print ( PROG_NAME + ":" + message )
        sys.exit(1)

def diffFilter(d):
    return (d["expected"] != d["actual"]) & ((pan.isna(d["expected"]) | pan.isna(d["actual"]))^True)

def roundColumnToDecDigits(col, digits):
    exp = numpy.zeros(len(col))
    significand = col.copy()
    while any(significand > 10):
        exp += significand > 10
        significand[significand > 10] = significand[significand > 10] / 10
    return numpy.round(significand, decimals=digits - 1) * 10 ** exp

def compareSorted(dataExpected, dataActual, precision):
    for colName in dataExpected.columns:
        assertAndDie(dataExpected[colName].dtype == dataActual[colName].dtype, "data types do not match for column {0}".format(colName))

        if dataExpected[colName].dtype == float:
            roundExpected = roundColumnToDecDigits(dataExpected[colName], precision)
            roundActual = roundColumnToDecDigits(dataActual[colName], precision)
            #
            diff = (roundExpected - roundActual).abs().loc[lambda d : d > 0]
            msg = "column {0} differs by more than {1} decimal digits at the line {2}"
            assertAndDie(len(diff) == 0, msg.format(colName, precision, diff.index[0] if len(diff.index) > 0 else 0 ))
        else:
            compFrame = pan.DataFrame({"expected": dataExpected[colName], "actual": dataActual[colName]})
            diff = compFrame.loc[diffFilter]
            msg = "column {0} doesn't match at the line {1}"
            assertAndDie(len(diff) == 0, msg.format(colName, diff.index[0] if len(diff.index) > 0 else 0 ))

class ComparableRow(object):
    def __init__(self, row, precision):
        self.row_data = row
        comp_data = []
        for val in row:
            if "dtype" in dir(val) and val.dtype == float:
                comp_data.append(self.roundToDecDigit(val, precision))
            else:
                comp_data.append(val)
        self.row_hash = tuple(comp_data).__hash__()
    #
    def roundToDecDigit(self, x, digits):
        exp = 0
        significand = x
        while significand > 10:
            exp += 1
            significand /= 10
        return round(significand, digits - 1) * 10 ** exp


    def __hash__(self):
        return self.row_hash
    #
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


def compareUnsorted(dataExpected, dataActual, precision):
    expected_only = set()
    actual_only = set()
    for i in dataExpected.index:
        #compare directly first
        compFrame = pan.DataFrame({"expected":dataExpected.loc[i, :], "actual":dataActual.loc[i, :]})
        if len(compFrame.loc[diffFilter]) == 0: continue
        er = ComparableRow(dataExpected.loc[i, :], precision)
        ar = ComparableRow(dataActual.loc[i, :], precision)
        if er == ar: continue
        if ar in expected_only:
            expected_only.remove(ar)
        else:
            actual_only.add(ar)
        if er in actual_only:
            actual_only.remove(er)
        else:
            expected_only.add(er)
    msg = "There are {0} rows in the actual data that do not match {1} rows in the expected data."
    assertAndDie(not expected_only and not actual_only, msg.format(len(expected_only), len(actual_only)))




parser = argparse.ArgumentParser("Comparator for crux output files for testing purposes.")
parser.add_argument("-p", "--check_python", action="store_true", dest="checkPython", help="Verify that python executable is available.")
parser.add_argument("-d", "--check_pandas", action="store_true", dest="checkPandas", help="Verify that pandas package is installed.")
parser.add_argument("-u", "--unsorted", action="store_true", dest="isUnsorted", help="Do unsorted rows comparison.")
parser.add_argument("-r", "--precision", type=int, dest="precision", metavar="p", help="Number of digits after the decimal point to use for numeric comparisons.")
parser.add_argument("files", action="store", nargs="*")

args = parser.parse_args()

if args.checkPython:
    exit(0)
if args.checkPandas:
    import pandas as pan
    exit(0)
    
import numpy
import pandas as pan

precision = 5
if os.environ.get(PRECISION_VAR) != None:
    precision = os.environ.get(PRECISION_VAR)
if args.precision != None:
    precision = args.precision

assertAndDie(len(args.files) == 2, "2 files are required")
assertAndDie(os.path.isfile(args.files[0]), "{0} is not a file".format(args.files[0]))
assertAndDie(os.path.isfile(args.files[1]),  "{0} is not a file".format(args.files[1]))

dataExpected = pan.read_csv(args.files[0], sep="\t")
dataActual = pan.read_csv(args.files[1], sep="\t")

assertAndDie( (dataExpected.columns == dataActual.columns).all(), "column names do not match.")
assertAndDie( len(dataExpected.index) == len(dataActual.index), "row counts do not match.")

if args.isUnsorted:
    compareUnsorted(dataExpected, dataActual, precision)
else:
    compareSorted(dataExpected, dataActual, precision)

exit(0)
