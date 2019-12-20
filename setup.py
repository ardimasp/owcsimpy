
import os.path
import sys
from setuptools import setup, Extension, find_packages

import numpy
import sys

import platform

# Using clang is slower
import os
# os.environ["CC"] = "clang"
# os.environ["CC"] = "icc"
# os.environ["LD_SHARED"] = "icc -shared"

try:
    from Cython.Build import cythonize
    has_cython = True
except ImportError:
    has_cython = False

with open("README.rst", "r") as fh:
    long_description = fh.read()

# info = {}
# VERSION = info['__version__']

print ('Argument List:', str(sys.argv))

use_cython = '--cython' in sys.argv or '--with-cython' in sys.argv
if '--cython' in sys.argv:
    sys.argv.remove('--cython')
if '--with-cython' in sys.argv:
    sys.argv.remove('--with-cython')

if use_cython and not has_cython:
    raise RuntimeError('Cython is required to build owcsimpy.')

if use_cython:
    suffix = '.pyx'
# else:
#     suffix = '.c'

ext_modules = []

if platform.system() == 'Linux':
    """
    Need an extra link arg '-L/usr/lib/x86_64-linux-gnu/' as 
    -lpthread and -lc not found issue

    I don't need this for Ubuntu shell in Windows 10, but 
    I need it for standalone Ubuntu.

    """
    ext_modules.append(
        Extension("owcsimpy.geoutils.cutils", 
            ["owcsimpy/geoutils/cutils.pyx","owcsimpy/geoutils/utils.c"],
            include_dirs=[numpy.get_include(),
            "owcsimpy/geoutils/"],
            extra_compile_args=["-DCYTHON_WITHOUT_ASSERTIONS"],
            extra_link_args=['-L/usr/lib/x86_64-linux-gnu/']
            )
        )
    ext_modules.append(
        Extension("owcsimpy.cir.cirutils", 
            ["owcsimpy/cir/cirutils.pyx","owcsimpy/cir/cirutils_c.c",
            "owcsimpy/geoutils/utils.c"],
            include_dirs=[
            numpy.get_include(),
            "owcsimpy/cir/","owcsimpy/geoutils/"],
            extra_compile_args=["-DCYTHON_WITHOUT_ASSERTIONS","-fopenmp"],
            extra_link_args=["-fopenmp",'-L/usr/lib/x86_64-linux-gnu/']
            )
        )
else:
    # Tested on 'Darwin' (MacOS) 
    ext_modules.append(
        Extension("owcsimpy.geoutils.cutils", 
            ["owcsimpy/geoutils/cutils.pyx","owcsimpy/geoutils/utils.c"],
            include_dirs=[numpy.get_include(),
            "owcsimpy/geoutils/"],
            extra_compile_args=["-DCYTHON_WITHOUT_ASSERTIONS"]
            )
        )
    ext_modules.append(
        Extension("owcsimpy.cir.cirutils", 
            ["owcsimpy/cir/cirutils.pyx","owcsimpy/cir/cirutils_c.c",
            "owcsimpy/geoutils/utils.c"],
            include_dirs=[
            numpy.get_include(),
            "owcsimpy/cir/","owcsimpy/geoutils/"],
            extra_compile_args=["-DCYTHON_WITHOUT_ASSERTIONS","-fopenmp"],
            extra_link_args=["-fopenmp"]
            )
        )



if use_cython:
    try:
        from Cython.Compiler.Options import get_directive_defaults
        directive_defaults = get_directive_defaults()
    except ImportError:
        # for Cython < 0.25
        # from Cython.Compiler.Options import directive_defaults
        raise RuntimeError('Failed cythonizing!')
    directive_defaults['embedsignature'] = True
    directive_defaults['binding'] = True
    ext_modules = cythonize(ext_modules,annotate=True,compiler_directives={'language_level':3})
    # ext_modules = cythonize(ext_modules,annotate=True)


setup(
    name="owcsimpy",
    version="0.0.1",
    author="ardimasp",
    author_email="ardimasandipurwita@outlook.com",
    description="Python simulator for Optical Wireless Communications",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="NaN",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    ext_modules=ext_modules,
    zip_safe=False,
)
