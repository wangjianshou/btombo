import io
from os import path

class F5BytesIO(io._io.BytesIO):
  def __init__(self, name='', is_new=False, initial_bytes=None):
    #assert isinstance(name, str), "name must be string"
    self.name = name
    if is_new:
      super(__class__, self).__init__(initial_bytes)
    else:
      with io.open(name, 'rb') as f:
        super(__class__, self).__init__(initial_bytes=f.read())
  def __add__(self, value):
    if isinstance(value, str):
      return self.name + value
    elif isinstance(value, __class__):
      return self.name + value.name
  def __radd__(self, value):
    if isinstance(value, str):
      return value + self.name
    elif isinstance(value, __class__):
      return value.name + self.name


def getF5fromName(names, f5s):
  # names must be list or tuble or set
  return [i for i in f5s if i.name in names]
