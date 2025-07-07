import os
import re
import subprocess
import sys



######################################
# This approach assumes compiled Python C++ API and only check it is available.
######################################
from setuptools import setup, find_packages
import sys

try:
    import flow123d_python_api
except ImportError:
    sys.exit("Error: flow123d_python_api module not found. Please ensure the compiled library is installed.")


######################################
# Following approach tries to compile flow123d through CMake.
######################################


from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

__version__ = '1.0.1'



# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # In this example, we pass in the version to C++. You might not need to.
        cmake_args += [f"-DEXAMPLE_VERSION_INFO={self.distribution.get_version()}"]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja  # noqa: F401

                    ninja_executable_path = os.path.join(ninja.BIN_DIR, "ninja")
                    cmake_args += [
                        "-GNinja",
                        f"-DCMAKE_MAKE_PROGRAM:FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:

            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        build_temp = os.path.join(self.build_temp, ext.name)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)



######################################
# Main configuration part.
######################################

# The information here can also be placed in setup.cfg - better separation of
# logic and declaration, and simpler if you include description/version in a file.
setup(
    name="flow123d",
    version=__version__,
    author='Jan Brezina',
    author_email='jan.brezina@tul.cz',    
    description="Flow123d extension",
    long_description="",
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("flow123d")],
#    cmdclass={"build_ext": CMakeBuild},
#    zip_safe=False,
#    extras_require={"test": ["pytest>=6.0"]},
    python_requires=">=3.6",
    license='GPL 3.0',
    url='https://github.com/flow123d/bih',
    packages=find_packages("src"),
    package_dir={"": "src"},
    entry_points={
        "console_scripts": [
            "flow123d=flow123d.__main__:main",
            # other user scripts
        ]
    },
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers        
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python :: 3',        
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Scientific/Engineering',
    ],
        
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
)

#############################################################################
# ORIG
#############################################################################


# import glob
# import setuptools


# class get_pybind_include(object):
#     """Helper class to determine the pybind11 include path
#     The purpose of this class is to postpone importing pybind11
#     until it is actually installed, so that the ``get_include()``
#     method can be invoked. """
# 
#     def __init__(self, user=False):
#         self.user = user
# 
#     def __str__(self):
#         #print("CWD bind:", os.getcwd())
#         import pybind11
#         return pybind11.get_include(self.user)
# 
# 
# 
# def get_sources(root, patterns):
#     #print("CWD :", os.getcwd())
#     sources = []
#     for p in patterns:
#         for path in glob.glob(os.path.join(root, p)):
#             print("Path: ", path)
#             sources.append(path)
#     return sources
# 
# ext_modules = [
#     setuptools.Extension(
#         'flowpy',
#         get_sources('src/fields', ['python_field_base.cc']),
#         include_dirs=[
#             'src',
#             # Path to pybind11 headers
#             get_pybind_include(),
#             get_pybind_include(user=True)
#         ],
#         language='c++'
#     ),
# ]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
# def has_flag(compiler, flagname):
#     """Return a boolean indicating whether a flag name is supported on
#     the specified compiler.
#     """
#     import tempfile
#     with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
#         f.write('int main (int argc, char **argv) { return 0; }')
#         try:
#             compiler.compile([f.name], extra_postargs=[flagname])
#         except setuptools.distutils.errors.CompileError:
#             return False
#     return True
# 
# 
# def cpp_flag(compiler):
#     """Return the -std=c++[11/14] compiler flag.
#     The c++14 is prefered over c++11 (when it is available).
#     """
#     if has_flag(compiler, '-std=c++14'):
#         return '-std=c++14'
#     elif has_flag(compiler, '-std=c++11'):
#         return '-std=c++11'
#     else:
#         raise RuntimeError('Unsupported compiler -- at least C++11 support '
#                            'is needed!')
# 
# 
# class BuildExt(build_ext):
#     """A custom build extension for adding compiler-specific options."""
#     c_opts = {
#         'msvc': ['/EHsc'],
#         'unix': [],
#     }
# 
#     if sys.platform == 'darwin':
#         c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
# 
#     def build_extensions(self):
#         ct = self.compiler.compiler_type
#         opts = self.c_opts.get(ct, [])
#         if ct == 'unix':
#             opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
#             opts.append(cpp_flag(self.compiler))
#             if has_flag(self.compiler, '-fvisibility=hidden'):
#                 opts.append('-fvisibility=hidden')
#         elif ct == 'msvc':
#             opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
#         for ext in self.extensions:
#             ext.extra_compile_args = opts
#         build_ext.build_extensions(self)


        
# with open("README.md", "r") as fh:
#     long_description = fh.read()


# setuptools.setup(
#     packages = setuptools.find_packages('src'),
#     package_dir={'': 'src'},
#     py_modules=[os.path.splitext(os.path.basename(path))[0] for path in glob.glob('src/*.py')],
#     package_data={
#         # If any package contains *.txt or *.rst files, include them:
#         #'': ['*.txt', '*.rst'],
#         # And include any *.msg files found in the 'hello' package, too:
#         #'hello': ['*.msg'],
#     },
# 
#     # include automatically all files in the template MANIFEST.in
#     include_package_data=True,
#     zip_safe=False,
#     install_requires=['pybind11>=2.2'],
#     python_requires='>=3',
#     extras_require={
#         # eg:
#         #   'rst': ['docutils>=0.11'],
#         #   ':python_version=="2.6"': ['argparse'],
#     },
#     # entry_points={
#     #     'console_scripts': [
#     #         'nameless = nameless.cli:main',
#     #     ]
#     # },
# 
#     ext_modules=ext_modules,
#     cmdclass={'build_ext': BuildExt}
#     #test_suite='test.pytest_bih'
# )        
        
        
        
