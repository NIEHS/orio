import os
import re
import sys
from setuptools import setup, find_packages


# get version from __init__, mirroring request library:
with open('orio/__init__.py', 'r') as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                        fd.read(), re.MULTILINE).group(1)

if sys.argv[-1] == 'publish_test':
    os.system('python setup.py sdist upload -r https://testpypi.python.org/pypi')
    os.system('python setup.py bdist_wheel upload -r https://testpypi.python.org/pypi')
    sys.exit()

if sys.argv[-1] == 'publish_production':
    os.system('python setup.py sdist upload')
    os.system('python setup.py bdist_wheel upload')
    sys.exit()

if sys.argv[-1] == 'tag':
    os.system("git tag -a %s -m 'version %s'" % (version, version))
    os.system("git push --tags")
    sys.exit()


def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='orio',
    version=version,
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
        'requests',
        'pymysql'
    ],
    include_package_data=True,
    zip_safe=False
)
