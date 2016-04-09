import numpy

# _landau is the C extension that does the heavy lifting
import _landau

def gauss(x_vals, gauss_mean, gauss_sigma, normed=False):
    '''
    This function applies the gaussian (normal) distribution
    to x_vals.  It always returns such that the sum,
    properly weighted with the bin widths, is 1.
    However, use the normed keyword to force the output
    to sum to 1 regardless of the bin widths.
    '''
    y_vals = _landau.gauss(x_vals,gauss_mean,gauss_sigma)
    if normed:
      y_vals /= numpy.sum(y_vals)
      return y_vals
    else:
      return y_vals



def landau(x_vals, landau_mu, landau_sigma, scale=1.0):
    '''
    This function applies the landau distribution
    to x_vals.  The mu parameter is approximately the 
    most probable value of the distribution, and the 
    sigma is a horizontal scale factor.  The true sigma
    of a landau distribution is undefined.
    '''
    y_vals = _landau.landau(x_vals,landau_mu,landau_sigma)
    if scale != 1.0:
      y_vals *= scale
      return y_vals
    else:
      return y_vals

def gauss_landau(x_vals, landau_mu, landau_sigma, gauss_sigma, scale=1.0):
    '''
    This function applies the landau distribution, convolved
    with a gaussian of width gauss_sigma, to x_vals.
    The mu parameter is approximately the most
    probable value of the distribution, and the 
    sigma is a horizontal scale factor.  The true sigma
    of a landau distribution is undefined.
    '''
    y_vals = _landau.gausslandau(x_vals,landa_mu,landau_sigma, gauss_sigma)
    if scale != 1.0:
      y_vals *= scale
      return y_vals
    else:
      return y_vals