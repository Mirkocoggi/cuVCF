from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools

__version__ = '0.0.1'

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """
    def __str__(self):
        import pybind11
        return pybind11.get_include()

ext_modules = [
    Extension(
        'VCFparser_mt_col_struct',
        ['src/VCFparser_bindings.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include()
        ],
        language='c++',
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp']
    ),
]

setup(
    name='vcfparser',
    version=__version__,
    author='Mirko Coggi, Riccardo Fiorentini',
    description='Python binding of VCFparser_mt_col_struct',
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.5.0'],
    cmdclass={'build_ext': build_ext},
    zip_safe=False,
)
