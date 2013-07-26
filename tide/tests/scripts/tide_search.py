import os.path, testutil

testutil.ImportPredefined(globals())

TIDE_DIR = testutil.DirWalkUp(WORKING_DIR, "tide")
BINROOT=os.path.join(TIDE_DIR, "bin-opt-%s" % ARCH)
INDEXROOT=os.path.join(TIDE_DIR, "index")
SPECTRAROOT=os.path.join(TIDE_DIR, "data")

class Task:
  MASS_WINDOW = 3
  RESULTS = "text"
  RequiredVariables = ['ORGANISM', 'PEPTIDE_SET', 'DB', 'SPECTRUM_FILE']

  def RequiredFiles(self):
    SPECRECS = not PARAMS.has_key('ms2')
    self.executable = "search"
    self.spec_suffix = SPECRECS and ".spectrumrecords" or ".ms2"
    self.peptides_file = "%s.pepix" % self.DB
    self.spectrum_file = self.ORGANISM + "-" + self.SPECTRUM_FILE + self.spec_suffix

    files = [os.path.join(BINROOT, self.executable),
             os.path.join(INDEXROOT, self.ORGANISM, "raw_proteins.protix"),
             os.path.join(INDEXROOT, self.ORGANISM, self.PEPTIDE_SET, "local", self.peptides_file),
             os.path.join(SPECTRAROOT, self.ORGANISM, self.spectrum_file)]
    return files

  def CommandLine(self):
    PROFILER = ""
    if PARAMS.has_key('prof'):
      PROFILER = "CPUPROFILE=./cpuprofile "
    return (PROFILER + "./" + self.executable + " " +
            "--proteins=raw_proteins.protix " +
            ("--peptides=%s " % self.peptides_file) +
            "--spectra=" + self.spectrum_file +
            (" --mass_window=%s" % self.MASS_WINDOW) +
            " --results=%s" % self.RESULTS)
