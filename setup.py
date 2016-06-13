import os
import sys
from setuptools import setup, find_packages

from orio import __version__


if sys.argv[-1] == 'publish_test':
    os.system('python setup.py sdist upload -r https://testpypi.python.org/pypi')
    os.system('python setup.py bdist_wheel upload -r https://testpypi.python.org/pypi')
    sys.exit()

if sys.argv[-1] == 'publish_production':
    os.system('python setup.py sdist upload')
    os.system('python setup.py bdist_wheel upload')
    sys.exit()


def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='orio',
    version=__version__,
    description='ORIO',
    long_description=readme(),
    url='https://github.com/shapiromatron/orio',
    author='Andy Lavender, Andy Shapiro',
    author_email='add@nobody.com',
    license='TBD',
    packages=find_packages(exclude=['tests']),
    install_requires=[
        'click',
        'clint',
        'docopt',
        'scipy',
        'numpy',
        'requests'
    ],
    include_package_data=True,
    zip_safe=False
)
