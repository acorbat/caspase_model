from setuptools import setup, find_packages

setup(
    name='caspase_model',
    version='0.1.0',
    packages=find_packages(),
    url='https://github.com/acorbat/caspase_model/tree/master/',
    license='MIT',
    author='Agustin Corbat',
    author_email='acorbat@df.uba.ar',
    description='Caspase model modified at Grecco Lab.',
    install_requires=['pysb',
                      'earm @ git+ssh://git@github.com/acorbat/earm@v1.0#egg=earm']
)
