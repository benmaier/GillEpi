from setuptools import setup

setup(name='GillEpi',
      version='0.0.6',
      description='Provides classes to simulate epidemics on (potentially time varying) networks using a Gillespie stochastic simulation algorithm and the classic agent based algorithm.',
      url='https://www.github.com/benmaier/GillEpi',
      author='Benjamin F. Maier',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['GillEpi'],
      install_requires=[
          'numpy',
          'networkx',
      ],
      zip_safe=False)
