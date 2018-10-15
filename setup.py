import os
from setuptools import setup


here = os.path.abspath(os.path.dirname(__file__))
readme_file = os.path.join(here, 'README.rst')
mf_dir = os.path.join(here, 'microfem')
version_file = os.path.join(mf_dir, 'version.py')

mf_name = 'microfem'
mf_description = 'Simple FEM Package for Plate Structures'
mf_license = 'MIT'
mf_author = 'Steven Moore'
mf_author_email = 'steven.ian.moore@gmail.com'
#mf_url = 'http://github.com/simoore/micro-fem'
mf_install_requires = ['numpy', 'scipy', 'matplotlib']
mf_pacakges = ['microfem']

with open(version_file, encoding='utf-8') as f:
    exec(f.read())  
    mf_version = __version__

with open(readme_file, encoding='utf-8') as f:
    mf_long_description = f.read()

#tt_classifiers = [
#    'Environment :: Console',
#    'License :: OSI Approved :: MIT License',
#    'Operating System :: MacOS :: MacOS X',
#    'Operating System :: Microsoft :: Windows',
#    'Operating System :: POSIX',
#    'Programming Language :: Python :: Implementation :: CPython',
#    'Programming Language :: Python :: 2.7',
#    'Programming Language :: Python :: 3.3',
#    'Programming Language :: Python :: 3.4',
#    'Programming Language :: Python :: 3.5',
#    'Programming Language :: Python :: 3.6',
#    'Topic :: Utilities'
#]

setup(name=mf_name,
      version=mf_version,
      description=mf_description,
      long_description=mf_long_description,
      #url=mf_url,
      #classifiers=mf_classifiers,
      author=mf_author,
      author_email=mf_author_email,
      license=mf_license,
      packages=mf_pacakges,
      install_requires=mf_install_requires,
      zip_safe=False)
