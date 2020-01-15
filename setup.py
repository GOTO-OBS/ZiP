from setuptools import setup

setup(name='zogyp',
      version='1.4.5',
      description='A parallel version of ZOGY image subtraction for astronomy, with Proper Co-Addition and Image alignment',
      url='https://zogy-in-parallel.readthedocs.io/en/latest/index.html',
      author='Ryan Cutter',
      author_email='R.Cutter@warwick.ac.uk',
      license='MIT',
      packages=['zogyp'],
      package_data={'zogyp': ['configfls/*.txt']},
      include_package_data=True,
      entry_points={
        'console_scripts': [
            'ZOGY=ZOGY:run'
        ]
      }
      zip_safe=False)
