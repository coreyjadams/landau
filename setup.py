# from distutils.core import setup, Extension
# import numpy.distutils.misc_util

# setup(
#     ext_modules=[Extension("_landau", ["landau.c", "landau_funcs.c"])],
#     include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
# )

from setuptools import setup, Extension
import numpy.distutils.misc_util


def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='landau',
    ext_modules=[Extension("_landau", ["src/landau.c", "src/landau_funcs.c"])],
    include_dirs=[numpy.get_include(), "inc/"],
    version='0.1',
    description='Landau functions with fast gaussian convolutions.',
    url='http://github.com/storborg/funniest',
    author='Corey Adams',
    author_email='corey.adams@yale.edu',
    license='MIT',
    packages=['landau'],
    install_requires=[
        'numpy'
    ],
    zip_safe=False)
