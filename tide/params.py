import sys

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
