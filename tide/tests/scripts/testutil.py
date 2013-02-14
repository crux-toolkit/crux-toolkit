import os, os.path, re, types

g_predefined = {}

def ImportPredefined(d):
  d.update(g_predefined)

def Eval(s, locals):
  return re.sub('{{(.*?)}}', lambda mo:str(eval(mo.group(1), globals(), locals)), s)

def NormPath(path):
  if path:
    return os.path.normpath(os.path.abspath(os.path.expandvars(os.path.expanduser(path))))
  return None

def DirWalkUp(start, match):
  if type(match) != types.ListType:
    match = [match]
  head = start
  tail = ''
  while head != '/' and tail not in match:
    head, tail = os.path.split(head)
  if tail not in match:
    return None
  return os.path.join(head, tail)
