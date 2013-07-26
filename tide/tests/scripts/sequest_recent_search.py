import os.path, testutil

testutil.ImportPredefined(globals())

TIDE_DIR = testutil.DirWalkUp(WORKING_DIR, "tide")
DATA_DIR = os.path.join(TIDE_DIR, "data")
SEQUEST_BINDIR = "/net/noble/vol1/home/bdiament/sequest/recent"
SEQUEST_BIN = "sequest_20091120"

class Task:
  XCORR_ONLY = True
  MASS_WINDOW = 3.0
  MODS = False
  INDEX = False
  MONO_PRECURSOR = False
  REMOVE_PRECURSOR = False
  RequiredVariables = ['ORGANISM', 'MS2_FILE', 'PARTIAL_DIGEST']

  def RequiredFiles(self):
    params_template = open(os.path.join(SCRIPTS_DIR,
                                        "sequest_recent.params.template")).read()
    environ = Task.__dict__.copy()
    environ.update(self.__dict__)
    params_text = testutil.Eval(params_template, environ)

    files = [os.path.join(SEQUEST_BINDIR, SEQUEST_BIN),
             {"sequest.params": params_text},
             os.path.join(DATA_DIR, self.ORGANISM, self.MS2_FILE),
             os.path.join(DATA_DIR, self.ORGANISM, self.ORGANISM + ".fasta")]
    return files

  def CommandLine(self):
    return "./%s %s" % (SEQUEST_BIN, os.path.basename(self.MS2_FILE))
