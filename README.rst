Landau
------

This package contains c based routines to apply probability functions to python array-like objects

The functions included convert to numpy array structures underneath, so all return types are as numpy
arrays that match the input.

Currently, multidimensional arrays are only supported if it's a numpy array:

a = [ [1,2,3], [4,5,6] ] will not work, but

a = numpy.zeros((2,3)) will work.

You can also use single elements in the functions


To use the package, do:

