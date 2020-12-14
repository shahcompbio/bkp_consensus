from setuptools import setup, find_packages

setup(
    name='bkpconsensus',
    version='0.1',
    description='Breakpoint consensus',
    packages=find_packages(),
    entry_points = {
        'console_scripts': ['bkpconsensus=bkpconsensus.consensus:main'],
    }
)

