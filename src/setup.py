from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        'vcfparser',
        ['vcfparser_bindings.cpp'],
        include_dirs=[pybind11.get_include()],
        language='c++'
    ),
]

setup(
    name='vcfparser',
    version='0.1',
    author='Mirko Coggi, Riccardo Fiorentini',
    description='Binding di VCFparser_mt_col_struct in Python',
    ext_modules=ext_modules,
    zip_safe=False,
)
