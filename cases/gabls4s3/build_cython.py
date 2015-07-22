# To build the Cython modules. Run as:
# python build_cython.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension("trilin", ["trilin.pyx"])]

setup(
  name = 'Cython Kernels',
  cmdclass = {'build_ext': build_ext},
  include_dirs = [np.get_include()],
  ext_modules = ext_modules
)
