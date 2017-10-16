# usage: python setup.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import os

os.environ["CC"] = "g++-7"
os.environ["CXX"] = "g++-7"

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ib_tools_cython",
                             sources=["ib_tools.pyx"],
                             language="c++",
                             extra_compile_args=["-Ofast", "-march=native", "-mtune=native", "-fno-wrapv"],
                             include_dirs=[numpy.get_include()])],
)
