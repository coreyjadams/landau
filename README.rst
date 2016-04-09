Landau
------

This package contains c based routines to apply probability functions to python array-like objects

The functions included convert to numpy array structures underneath, so all return types are as numpy
arrays that match the input.

Currently, multidimensional arrays are only supported if it's a numpy array:

a = [ [1,2,3], [4,5,6] ] will not work, but

a = numpy.zeros((2,3)) will work.  Also, only contiguous numpy arrays will work reliably.

You can also use single elements in the functions and that works.

There are three functions implemented:

gauss(x_vals, gauss_mean, gauss_sigma, normed=False) applies the normal distribution to your input

landau(x_vals, landau_mu, landau_sigma, scale=1.0) applies a landau function

gauss_landau(x_vals, landau_mu, landau_sigma, gauss_sigma, scale=1.0) applies a landau function convolved with a gaussian.


To use the package, do:


python setup.py build

python setup.py install
