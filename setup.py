from setuptools import setup

setup(name='zogyp',
      version='1.6.3',
      description='A parallel version of ZOGY image subtraction for astronomy, with Proper Co-Addition and Image alignment',
      url='https://zogy-in-parallel.readthedocs.io/en/latest/',
      author='Ry Cutter',
      author_email='R.Cutter@warwick.ac.uk',
      license='General Public License v3.0',
      packages=['zogyp'],
      package_data={'zogyp': ['configfls/*.txt']},
      include_package_data=True,
      entry_points={
        'console_scripts': [
            'ZOGY=zogyp.ZOGY:main'
        ]
      },
      zip_safe=False)
