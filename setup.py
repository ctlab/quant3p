from setuptools import setup

setup(name='quant3p',
      version='0.1',
      description='Tools to do 3\' RNA-seq quantification',
      url='http://github.com/ctlab/quant3p',
      author='Alexey Sergushichev',
      author_email='alserg@rain.ifmo.ru',
      license='Apache2',
      packages=['quant3p'],
      install_requires=[
          'pysam',
          'pybedtools',
          'macs2 >= 2.1.0'
          ],
      scripts=['bin/macs2-stranded', 'bin/quant3p'],
      entry_points = {
          'console_scripts': [
              'fix-mm=quant3p.fixmm:main',
              'gtf-extend=quant3p.gtfextend:main'],
          },
      zip_safe=False)
