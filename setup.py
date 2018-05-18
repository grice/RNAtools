from setuptools import setup, Extension
from Cython.Distutils import build_ext
import numpy as np

ext_1 = Extension('RNAtools.partAlign', 
                  ['RNAtools/partAlign.pyx'],
                  libraries=[],
                  include_dirs=[np.get_include()])


EXTENSIONS = [ext_1]
setup(name='RNAtools',
      version='0.8.3',
      description='Tools for working with RNA files.',
      author='Gregg M. Rice',
      author_email='gmr@unc.edu',
      license='MIT',
      packages=['RNAtools'],
      include_package_data=True,
      zip_safe=False,
      cmdclass={"build_ext": build_ext},
      ext_modules=EXTENSIONS,
      entry_points={
        'console_scripts': ['dot2ct=RNAtools.scripts.dot2ct:main',
                            'colorRNA=RNAtools.scripts.colorRNA:main',
                            'arcsRNA=RNAtools.scripts.arcsRNA:main',
                            'cleanCT=RNAtools.scripts.cleanCT:main']
      })
