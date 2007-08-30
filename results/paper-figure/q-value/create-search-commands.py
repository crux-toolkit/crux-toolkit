import sys
import os
# simple python script to submit to the cluster
MATCH_SEARCH_PATH="/nfs/gs/home/aklammer/crux/bin/match_search"

if len(sys.argv) < 3:
  raise SystemExit, "%s: <mass-range-file> [mass-analysis args]*" % sys.argv[0]

mass_range_file = sys.argv[1]
match_analysis_args = sys.argv[2:]
fh = open(mass_range_file)
cwd = os.getcwd()
for line in fh:
  start, end = line.rstrip("\n").split()
  cmdFile = os.getcwd() + os.sep + "match_analysis_%s_%s" % (start, end)
  outFh = open(cmdFile, "w")
  cmd = "%s --spectrum-min-mass %.3f --spectrum-max-mass %.3f %s \n" \
     % (MATCH_SEARCH_PATH, float(start), float(end), " ".join(match_analysis_args))
  outFh.write(cmd)
  outFh.close()
  clusterCmd = "qsub -o %s -e %s %s" % (cwd, cwd, cmdFile)
  print clusterCmd

