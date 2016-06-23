This directory houses executables to perform benchmark testing on Barista. A user
should run dispatch.sh to dispatch a numbe rof uge jobs. After the jobs have finished,
the user can run summarize.py to accumulate the output and generate graphs.
Theoretically, one can compare these results to stored plots (which have yet to be
generated) by opening barista-benchmark.html

TODO:
- Look through various todos in the scripts
- either configure this to find an recent build of crux or move src/crux to .
- Conduct a sample run, put plots in stored plots
- change the outputs of the scripts to correspond to this directory. Currently they
  output to ~/benchmark_tests. Not sure how to fix this. It may not be a problem if
  summarize merely consumes the temp files in ~/benchmark_tests and writes the graphs
  to trunk/test/barista/plots/

Files:
- dispatch.sh
  - main executable, dispatches uge jobs
- memory_job.sh
  - UGE handle for analysis of barista memory usage. Calls setup.sh and runs barista
- time_job.sh
  - UGE handle for analysis of barista time usage. Calls setup.sh and runs barista
- setup.sh
  - Shared code between *_job.sh. Moves files to a temp directory, unzips them, runs
    tide index and search
- summarize.py
  - Collects output and generates graphs
- barista-benchmark.html
  - Display output and compare to previous runs
- plots/
  - Where summarize.py writes plots
- stored-plots/
  - A known valid version of plots/
