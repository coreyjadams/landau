#include <Python.h>
#include <numpy/arrayobject.h>
#include "landau_funcs.h"

static char module_docstring[] =
  "This module provides an interface for making distribution functions from numpy arrays.";
static char gauss_docstring[] =
  "Calculate the gaussian of an array (also accepts single point)";
static char landau_docstring[] =
  "Calculate the landau distribution of an array (also accepts single point)";
static char gausslandau_docstring[] =
  "Calculate the landau distribution, convolved with a gaussian, of an array (also accepts single point)";

static PyObject *landau_gauss(PyObject *self, PyObject *args);
static PyObject *landau_landau(PyObject *self, PyObject *args);
static PyObject *landau_gausslandau(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
  {"gauss", landau_gauss, METH_VARARGS, gauss_docstring},
  {"landau", landau_landau, METH_VARARGS, landau_docstring},
  {"gausslandau", landau_gausslandau, METH_VARARGS, gausslandau_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_landau(void)
{
  PyObject *m = Py_InitModule3("_landau", module_methods, module_docstring);
  if (m == NULL)
    return 0;

  /* Load `numpy` functionality. */
  import_array();
  return 0;
}

static PyObject *landau_gauss(PyObject *self, PyObject *args)
{
  double mean, sigma;
  PyObject *x_obj;

  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "Odd", &x_obj, &mean, &sigma))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  if (x_array == NULL) {
    Py_XDECREF(x_array);
    return NULL;
  }

  PyObject *out_array = PyArray_NewLikeArray(x_array, NPY_ANYORDER, NULL, 0);

  /* If that didn't work, throw an exception. */
  if (out_array == NULL) {
    // Py_XDECREF decreases the reference counter only if the object is not null
    Py_XDECREF(out_array);
    return NULL;
  }

  /* How many data points are there? */
  int NDIM = (int) PyArray_NDIM(x_array);
  int N = 0;
  if (NDIM == 0) {
    // This is a 1D array, which means it can't be indexed
    // and you can't access the length of the dimensions
    N = 1;
  }
  else {
    N = (int) PyArray_Size(x_array);
  }

  /* Get pointers to the data as C-types. */
  double *x    = (double*)PyArray_DATA(x_array);
  double *y    = (double*)PyArray_DATA(out_array);

  /* Call the external C function to compute the chi-squared. */
  gauss(x, y, N, mean, sigma);

  /* Clean up. */
  Py_DECREF(x_array);
  // Py_DECREF(out_array);

  // if (value < 0.0) {
  //   PyErr_SetString(PyExc_RuntimeError,
  //                   "Chi-squared returned an impossible value.");
  //   return NULL;
  // }

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("O", out_array);
  return ret;
}


static PyObject *landau_landau(PyObject *self, PyObject *args)
{
  double mu, sigma;
  PyObject *x_obj;

  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "Odd", &x_obj, &mu, &sigma))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  if (x_array == NULL) {
    Py_XDECREF(x_array);
    return NULL;
  }

  PyObject *out_array = PyArray_NewLikeArray(x_array, NPY_ANYORDER, NULL, 0);

  /* If that didn't work, throw an exception. */
  if (out_array == NULL) {
    // Py_XDECREF decreases the reference counter only if the object is not null
    Py_XDECREF(out_array);
    return NULL;
  }

  /* How many data points are there? */
  int NDIM = (int) PyArray_NDIM(x_array);
  int N = 0;
  if (NDIM == 0) {
    // This is a 1D array, which means it can't be indexed
    // and you can't access the length of the dimensions
    N = 1;
  }
  else {
    N = (int) PyArray_Size(x_array);
  }

  /* Get pointers to the data as C-types. */
  double *x    = (double*)PyArray_DATA(x_array);
  double *y    = (double*)PyArray_DATA(out_array);

  /* Call the external C function to compute the chi-squared. */
  landau(x, y, N, mu, sigma, true);

  /* Clean up. */
  Py_DECREF(x_array);
  // Py_DECREF(out_array);

  // if (value < 0.0) {
  //   PyErr_SetString(PyExc_RuntimeError,
  //                   "Chi-squared returned an impossible value.");
  //   return NULL;
  // }

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("O", out_array);
  return ret;
}

static PyObject *landau_gausslandau(PyObject *self, PyObject *args)
{
  double mu_landua, sigma_landau, sigma_gauss;
  PyObject *x_obj;

  /* Parse the input tuple */
  if (!PyArg_ParseTuple(args, "Oddd", &x_obj, &mu_landua, &sigma_landau, &sigma_gauss))
    return NULL;

  /* Interpret the input objects as numpy arrays. */
  PyObject *x_array = PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_IN_ARRAY);

  if (x_array == NULL) {
    Py_XDECREF(x_array);
    return NULL;
  }

  PyObject *out_array = PyArray_NewLikeArray(x_array, NPY_ANYORDER, NULL, 0);

  /* If that didn't work, throw an exception. */
  if (out_array == NULL) {
    // Py_XDECREF decreases the reference counter only if the object is not null
    Py_XDECREF(out_array);
    return NULL;
  }

  /* How many data points are there? */
  int NDIM = (int) PyArray_NDIM(x_array);
  int N = 0;
  if (NDIM == 0) {
    // This is a 1D array, which means it can't be indexed
    // and you can't access the length of the dimensions
    N = 1;
  }
  else {
    N = (int) PyArray_Size(x_array);
  }

  /* Get pointers to the data as C-types. */
  double *x    = (double*)PyArray_DATA(x_array);
  double *y    = (double*)PyArray_DATA(out_array);

  /* Call the external C function to compute the chi-squared. */
  conv_gauss_landau(x, y, N, sigma_gauss, sigma_landau, mu_landua);

  /* Clean up. */
  Py_DECREF(x_array);
  // Py_DECREF(out_array);

  // if (value < 0.0) {
  //   PyErr_SetString(PyExc_RuntimeError,
  //                   "Chi-squared returned an impossible value.");
  //   return NULL;
  // }

  /* Build the output tuple */
  PyObject *ret = Py_BuildValue("O", out_array);
  return ret;
}