from distutils.core import setup, Extension
from Cython.Distutils import build_ext

extensions = [ Extension("pycolibri",
                ["unordered_map.pxd", "pycolibri_classes.pxd", "pycolibri_wrapper.pyx"],
                language='c++',
                include_dirs=['/home/proycon/local/include/colibri/', '/home/proycon/local/include/', '/usr/include/libxml2' ],
                library_dirs=['/home/proycon/local/lib/','/usr/include/lib'],
                libraries=['colibri'],
                extra_compile_args=['--std=c++0x'],
                ) ]

setup(
    name = 'pycolibri',
    ext_modules = extensions,
    cmdclass = {'build_ext': build_ext},
)
