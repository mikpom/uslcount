from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='uslcount',
      version='0.1.1.14',
      description='Assigning reads from BAM alignments to gene features',
      long_description=readme(),
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU GPL v3.0',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries'],
      keywords= ['Bioinformatics', 'NGS', 'GeneExpression'],
      url='https://gitlab.lji.org/mikhail/uslcount.git',
      author='Mikhail Pomaznoy',
      author_email='mikhail@lji.org',
      license='GNU GPL v3.0',
      packages=['uslcount'],
      zip_safe=False,
      install_requires=[
          'numpy',
          'pysam',
          'bx-python'],
      test_suite='uslcount.main_tests',
      tests_require=['setuptools'])
