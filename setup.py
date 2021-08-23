try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

json_files = package_files('hippounit/tests/stimuli')
default_NMDAr = package_files('hippounit/tests/default_NMDAr')

setup(
    name='hippounit',
    version='1.3.5.4',
    author='Sara Saray, Szabolcs Kali, Christian Rossert, Andrew Davison, Shailesh Appukuttan',
    author_email='saray.sara@koki.mta.hu, kali@koki.hu, christian.rossert@epfl.ch, andrew.davison@unic.cnrs-gif.fr, shailesh.appukuttan@unic.cnrs-gif.fr',
    packages=['hippounit', 'hippounit.tests', 'hippounit.capabilities', 'hippounit.scores'],
    package_data={'hippounit': json_files + default_NMDAr},
    url='http://github.com/kalilab/hippounit',
    license='MIT',
    description='A SciUnit library for data-driven validation testing of models of hippocampus.',
    long_description="",
    install_requires=['sciunit>=0.2.1', 'efel']
)
