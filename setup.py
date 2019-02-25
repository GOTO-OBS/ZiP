from setuptools import setup

setup(name='zogyp',
      version='1.4.5',
      description='A parallel version of ZOGY image subtraction for astronomy, with Proper Co-Addition and Image alignment',
      url='https://github.com/ryanc123/RyZiP',
      author='Ryan Cutter',
      author_email='R.Cutter@wawrcik.ac.uk',
      license='MIT',
      packages=['zogyp'],
      package_data={'zogyp': ['configfls/*.txt']},
      include_package_data=True,
      zip_safe=False)
