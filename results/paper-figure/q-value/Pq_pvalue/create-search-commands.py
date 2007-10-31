import sys
import os
# simple python script to submit to the cluster
MATCH_SEARCH_PATH=os.getcwd() + os.sep + ".." + os.sep + "match_search"

if len(sys.argv) < 3:
  raise SystemExit, "%s: <mass-range-file> [mass-search args]*" % sys.argv[0]

mass_range_file = sys.argv[1]
fasta_file = sys.argv[3]
base = os.path.basename(fasta_file).split(".")[0]
match_search_args = sys.argv[2:]
fh = open(mass_range_file)
cwd = os.getcwd()
for line in fh:
  start, end = line.rstrip("\n").split()
  tag = "%s_%s_%s" % (start, end, base)
  cmdFile = os.getcwd() + os.sep + "match_search_%s" % tag
  outFh = open(cmdFile, "w")
  cmd = "cd %s; %s --output-mode binary --number-decoy-set 0 --spectrum-min-mass %.3f --spectrum-max-mass %.3f %s >& %s.err \n" \
     % (os.getcwd(), MATCH_SEARCH_PATH, float(start), float(end), " ".join(match_search_args), tag)
  outFh.write(cmd)
  outFh.close()
  clusterCmd = "qsub -cwd -o %s -e %s %s" % (cwd, cwd, cmdFile)
  print clusterCmd

