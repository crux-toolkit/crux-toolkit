#!/usr/bin/env python

# Authors: Benjamin Diament, Sean McIlwain
# General framework for running benchmarks and correctness tests.
#
# Create a directory to work in, such as
# /scratch/bdiament/tide_tests/benchmarks/2010-09-02/, called the "working
# directory".  Results will be created in this directory.
#
# In the new directory, create a file task_list_file (or specify another file
# with --tasks=/some/other/task/file) with a list of tests and/or
# benchmarks. These tests will look like function calls to files found in
# tide/tests/scripts. For example, task_list_file may have the following two
# lines:
#
# tide_search(ORGANISM='worm', PEPTIDE_SET='partially_tryptic', DB='with_peaks' SPECTRUM_FILE='worm-06-10000', MASS_WINDOW=3.0)
# tide_search(ORGANISM='yeast', PEPTIDE_SET='fully_tryptic', DB='with_peaks' SPECTRUM_FILE='yeast-02-10000', MASS_WINDOW=3.0)
#
# which refer to tide_search.py in tide/tests/scripts.  Those are the tests that
# will be run, and the calling format is discussed below (see TASK SCRIPTS).
# Any line in task_list_file beginning with a hash (#) is a comment.
#
# From the new working directory, run this script. (Alternatively, run this
# script from anywhere, specifying
# --wd=/scratch/bdiament/tide_tests/benchmarks/2010-09-02/ on the command line.)
#
# Our working directory will now have one subdirectory for each test or
# benchmark. In it, all needed files will be linked, and a test or benchmark
# executed. If the execution succeeded, a file execution_time will be created
# there, with the timing results of running the test. Standard output will be
# placed in task_stdout. For now the stderr of the task will appear in the same
# file as execution_time.  TODO: make stderr go to its own file.  TODO: collect
# timing results from all tasks.
#
# If some test fails, make the needed correction to the code and simply rerun
# runner.py, which will rerun/retime only those tasks that haven't completed
# successfully.
#
# TASK SCRIPTS
#
# As explained above, each line in task_list_file refers to a script in
# SCRIPTS_DIR with parameters. The script must define a class Task, with
# methods RequiredVariables, RequiredFiles and CommandLine. RequiredVariables
# defines a minimal list of variables that must be set in the
# task_list_file. RequiredFiles returns a list of required files to link
# to. (Each entry in the list may be a pair, with the second element being the
# name of the symlink in the task directory.) CommandLine just indicates what to
# run. The Task class will be instantiated by the framework, and when it runs
# RequiredFiles and CommandLine, it will have access to the following instance
# variables:
#
# self.ARCH: machine architecture as supplied environment variable ARCH
#            if defined, or output of "uname -m" otherwise
# self.PARAMS: a dictionary of command-line flags as passed to runner.py
# self.WORKING_DIR: as above
# self.SCRIPTS_DIR: as above
# TODO: make these available even upon object construction
#
# TODO: Document advice on task script organization

DONE_TOKEN_NAME = "execution_time"

import sys, os, os.path, subprocess, re, testutil

def ParseParameters():
  param_dict = {}
  args = [term for term in sys.argv[1:] if term [:2] != '--']
  params = [param[2:] for param in sys.argv[1:] if param[:2] == '--']
  params = [param.split('=') for param in params]
  for param in params:
    if len(param) == 2:
      name, val = param
      param_dict[name] = val
    elif len(param) == 1:
      param_dict[param[0]] = None
  return param_dict, args


PARAMS, ARGS = ParseParameters()
WORKING_DIR = testutil.NormPath(PARAMS.get('wd', '.'))
SCRIPTS_DIR = testutil.NormPath(PARAMS.get('scripts_dir', sys.path[0] or WORKING_DIR))
ARCH = os.environ.get('ARCH') or subprocess.Popen(["uname", "-m"], stdout=subprocess.PIPE).communicate()[0].strip()

testutil.g_predefined.update({
  'ARCH': ARCH,
  'PARAMS': PARAMS.copy(),
  'ARGS': [i for i in ARGS],
  'WORKING_DIR': WORKING_DIR,
  'SCRIPTS_DIR': SCRIPTS_DIR
  })


def CreateDir(dir_name):
  # returns:
  # 0 if run already done
  # 1 if previously started, but failed
  # 2 if brand new
  assert(os.path.sep not in dir_name) # test directories must appear immediately under the Working dir.
  dir_name = os.path.join(WORKING_DIR, dir_name)
  if os.path.isdir(dir_name):
    if os.path.exists(os.path.join(dir_name, DONE_TOKEN_NAME)):
      return 0
    return 1
  os.mkdir(dir_name)
  return 2

class TaskSpecException(Exception):
  pass


class Task:
  def __init__(self, match_dict):
    self.task_name = match_dict['task_name']
    self.params = eval('(lambda **kwargs: kwargs)(%s)' % match_dict['params'], {}, {})
    self.count = int('0' + match_dict['count']) or 1
    temp_sys_path = sys.path
    self.task_class = None
    try:
      sys.path = [SCRIPTS_DIR] + sys.path
      script = __import__(self.task_name)
      self.task_class = script.Task
    finally:
      sys.path = temp_sys_path
    if not self.task_class:
      raise TaskSpecException("Can't find script %s/%s.py with class Task" %
                              (SCRIPTS_DIR, self.task_name))
    required_vars = self.task_class.RequiredVariables
    for req in required_vars:
      if req not in self.params.keys():
        raise TaskSpecException("Required parameter %s not specified." % req)
    allowable_vars = required_vars + self.task_class.__dict__.keys()
    for param in self.params.keys():
      if param not in allowable_vars:
        raise TaskSpecException("Unrecognized parameter %s." % param)

  def GetName(self):
    return self.task_name + ''.join(["_%s_%s" % (p, self.params[p]) for p in self.params])

  def GetInstance(self):
    instance = self.task_class()
    instance.__dict__.update(self.params)
    return instance

  def Instances(self):
    name = self.GetName()
    suffixes = [""] + ["_%d" % (i+2) for i in range(self.count - 1)]
    return [(name + suffix, self.GetInstance()) for suffix in suffixes]

    
def GetTaskList():
  pattern = re.compile(r'^(?P<task_name>\w+)[(](?P<params>[^)]*)[)]\s*(?P<count>\d*)$')
  task_list_file = PARAMS.get('tasks', 'task_list_file')
  tasks = [line.strip() for line in open(task_list_file).readlines()]
  tasks = [line for line in tasks if line[0] != '#']
  tasks = [Task(pattern.match(line).groupdict()) for line in tasks]
  return [i for task in tasks for i in task.Instances()]


def LinkRequiredFiles(task_name, required_files):
  for file in required_files:
    if type(file) == dict:
      dest = file.keys()[0]
      text = file[dest]
      dest = os.path.join(WORKING_DIR, task_name, dest)
      if os.path.exists(dest):
        os.remove(dest) # in case there was a symlink here.
      open(dest, "w").write(text)
    else:
      if type(file) == tuple:
        file, dest = file
      else:
        head, dest = os.path.split(file)
      dest = os.path.join(WORKING_DIR, task_name, dest)
      if os.path.islink(dest):
        os.remove(dest)
      if os.path.exists(dest):
        raise Exception("Failed to create link to %s at %s: %s exists and is not a symlink.")
      os.symlink(file, dest)


def TimeRun(task_name, command_line, subtask_num):
  task_dir = os.path.join(WORKING_DIR, task_name)
  subtask_prefix = "subtask%d_" % subtask_num if subtask_num else ""
  done_token_name = os.path.join(task_dir, "%s%s" % (subtask_prefix, DONE_TOKEN_NAME))
  task_stdout_name = os.path.join(task_dir, "%stask_stdout" % subtask_prefix)
  timing_info_name = os.path.join(task_dir, "%stask_timing" % subtask_prefix)
  subtask_name = "%s%s" % (task_name, " subtask %d" % subtask_num if subtask_num else "")
  if os.path.exists(os.path.join(task_dir, done_token_name)):
    sys.stderr.write("%s: ALREADY RUN. Skipping.\n" % subtask_name)
    return
  if os.path.exists(os.path.join(task_dir, timing_info_name)):
    sys.stderr.write("%s: RERUNNING..." % subtask_name)
  else:
    sys.stderr.write("%s: RUNNING..." % subtask_name)
  task_stdout = open(task_stdout_name, "w")
  timing_info = open(timing_info_name, "w")
  open(os.path.join(task_dir, "%scommand_line" % subtask_prefix), "w").write("%s\n" % command_line)
  retval = None
  if PARAMS.has_key("n"):
    # don't actually do the run -- just show the command line
    sys.stderr.write("\n")
    print command_line
  else:
    retval = subprocess.call("time " + command_line, shell=True, cwd=task_dir,
                             stdout=task_stdout, stderr=timing_info)
  task_stdout.close()
  timing_info.close()
  if retval == None:
    pass
  elif retval == 0:
    sys.stderr.write("SUCCEEDED\n")
    os.rename(timing_info_name, done_token_name)
  else:
    sys.stderr.write("FAILED WITH EXIT CODE %d\n" % retval)


def RunTest(task_name, task_object):
  LinkRequiredFiles(task_name, task_object.RequiredFiles())
  commands = task_object.CommandLine()
  if type(commands) == type([]):
    for command_line, subtask_num in zip(commands, range(len(commands))):
      TimeRun(task_name, command_line, subtask_num + 1)
  else:
    TimeRun(task_name, commands, None)


def RunTests():
  os.chdir(WORKING_DIR)
  task_list = GetTaskList()
  for task_name, task_object in task_list:
    CreateDir(task_name)
    RunTest(task_name, task_object)


# main
RunTests()
