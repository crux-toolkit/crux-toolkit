import os.path, testutil

testutil.ImportPredefined(globals())

TIDE_DIR = testutil.DirWalkUp(WORKING_DIR, "tide")
DATA_DIR = os.path.join(TIDE_DIR, "data")
SEQUEST_DIR = "/net/noble/vol1/home/bdiament/sequest"
SEQUEST_BIN = "search28_mine"


class Task:
  MASS_WINDOW = 3.0
  MODS = False
  RequiredVariables = ['ORGANISM', 'DTA_FILE_PATTERN']

  def RequiredFiles(self):
    params_template = open(os.path.join(SCRIPTS_DIR,
                                        "sequest_orig.params.template")).read()
    environ = Task.__dict__.copy()
    environ.update(self.__dict__)
    params_text = testutil.Eval(params_template, environ)

    files = [os.path.join(SEQUEST_DIR, SEQUEST_BIN),
             {"sequest.params": params_text},
             os.path.join(SEQUEST_DIR, "dta"),
             os.path.join(DATA_DIR, self.ORGANISM, self.ORGANISM + ".fasta")]
    return files

  def CommandLine(self):
    return "./%s dta/%s-%s.*.dta" % (SEQUEST_BIN, self.ORGANISM, self.DTA_FILE_PATTERN)
