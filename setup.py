from setuptools import setup

with open('README.rst', 'r') as fh:
    mf_long_description = fh.read()

mf_name = 'microfem'
mf_description = 'Simple FEM Package for Plate Structures'
mf_license = 'MIT'
mf_author = 'Steven Moore'
mf_author_email = 'steven.ian.moore@gmail.com'
mf_install_requires = ['numpy', 'scipy', 'matplotlib']
mf_pacakges = ['microfem']
mf_version = '0.1.0'


setup(name=mf_name,
      version=mf_version,
      description=mf_description,
      author=mf_author,
      long_description=mf_long_description,
      author_email=mf_author_email,
      license=mf_license,
      packages=mf_pacakges,
      install_requires=mf_install_requires,
      zip_safe=False)
