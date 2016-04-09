try:
  import _landau
except:
  raise Exception("Need to build the package first!")

__all__ = ['landau']

from landau import gauss, landau, gauss_landau

