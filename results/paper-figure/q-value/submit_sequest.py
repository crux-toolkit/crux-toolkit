#!/usr/bin/env python
import sys
import os
import glob

zsPerJob=3000
scriptPath = "/nfs/gs/home/aklammer/crux/results/paper-figure/q-value/"
thisScript = "%s/submit_sequest.py" % scriptPath

def writeScript(scriptName,cmd):
  job_file = open(scriptName, "w")
  job_file.write("#!/bin/csh -efx\n")
  job_file.write("#$ -l h_cpu=6:00:00\n");
  job_file.write("hostname\n");
  job_file.write("cd %s\n" % os.getcwd())
  job_file.write("%s\n" % cmd)
  job_file.write("rm %s\n" % scriptName)
  job_file.close()


def submitJob(pos,n,tmpDir,ms2file):
  jobFileName = "%s/sequest-%04d.csh" % (tmpDir,n)
  cmd = "%s %s %d %d" % (thisScript,ms2file,pos,n)
  writeScript(jobFileName,cmd)
  jname = "s-%04d" % n
  command = "qsub -cwd -N %s %s" % (jname,jobFileName)

  # Run the command, or submit the job.
  sys.stderr.write("Running (%s)\n" % command)
  pipe=os.popen(command)
  # extract jobID
  line=pipe.readline()
  jid=line.split()[2]
  return jid

def makeSqt(ms2File,sqtFile,tmpDir):
  file = open (sqtFile,"w")
  file.write("H\tSQTGenerator SEQUEST\nH\tSQTGeneratorVersion\t2.7\n")
  file.close
#  cmd="grep '^H' %s  >> %s" % (ms2File,sqtFile)
#  os.system(cmd)
  fileList=glob.glob("%s/seq-*.sqt" % tmpDir)
  fileList.sort()
  i=0
  missing=0
  for file in fileList:
    i+=1
    fno=int(file[-8:-4])
    if(fno!=i):
      print "File %s/seq-%04d.sqt is missing" % (tmpDir,i)
      missing+=1
      i=fno
    os.system("cat %s >> %s" % (file,sqtFile))
  if (missing==0):
    os.system("rm -r %s" % tmpDir)


if (len(sys.argv)<2 or len(sys.argv)>4):
  print "Wrong number of arguments"
  sys.exit(-1)

ms2File = sys.argv[1]
tmpDir="tmp-" + ms2File[:-4]

if (len(sys.argv)==2):
  # Submit scripts mode
  mode=0
  submited=[]
elif (len(sys.argv)==4):
  # Execution mode, process a part of the ms2 with search27
  mode=1
  firstpos=float(sys.argv[2])
  pos=firstpos-100
  if pos<0: pos=0
  n=float(sys.argv[3])
  outFileName="%s/seq-%04d.out" % (tmpDir,n)
  sqtFileName="%s/seq-%04d.sqt" % (tmpDir,n)
elif (len(sys.argv)==3):
  # concatinate all the different parts into one sqt-file
  mode=2
  makeSqt(ms2File,sys.argv[2],tmpDir)
  sys.exit(0)

if mode==0:
  os.mkdir(tmpDir)

file = open(ms2File,"r")
if mode==1:
  file.seek(pos)
noZLines=0
noJobs=0
sscan,escan,z,info="0","0",0,[]
while 1:
  pos=file.tell()
  line=file.readline()
  if not line or len(line)==0:
    file.close()
    break
  if line[0]=='S' and mode==1:
    word = line.split()
    sscan = word[1]
    escan = word[2]
  if line[0]=='Z':
    if noZLines%zsPerJob==0 and mode==0:
      noJobs+=1
      submited+=[submitJob(pos,noJobs,tmpDir,ms2File)]
    if mode==1 and pos>=firstpos:
      word = line.split()
      z = word[1]
      info += [(sscan,escan,z)]
      os.system("%s/search27 %s:%d >> %s" % (scriptPath,ms2File,pos,outFileName))
      if noZLines+1>=zsPerJob:
        file.close()
        break
    noZLines+=1
if mode==0:
  jobFileName = "%s/sequest-end.csh" % tmpDir
  cmd= "%s %s %s.sqt\n" % (thisScript,ms2File,ms2File[:-4])
  writeScript(jobFileName,cmd)
  dep=",".join(submited)
  os.system("qsub -cwd -N s-end -hold_jid %s %s" % (dep,jobFileName))
if mode==1:
  # convert out to sqt
  outFile = open(outFileName, "r")
  lines = outFile.readlines()
  outFile.close()
  sqtFile = open(sqtFileName, "w")
  ix=0
  for line in lines:
    if line[0]=='S':
      line="\t".join(["S",info[ix][0],info[ix][1],info[ix][2]] + line.split()[4:]) + "\n"
      ix+=1
    sqtFile.write(line)
  sqtFile.close()

