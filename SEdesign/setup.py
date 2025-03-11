from setuptools import setup, find_packages

setup(
    name='se_design',
    version='0.1',
    packages=find_packages(),
    description='DNA sequence design tools.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Antimonyleo/SEdesign',
    author='Hao Liu and Matthew Sample',
    author_email='matsample1@gmail.com',
    license='MIT',
    install_requires=[
    'numba',
    'scipy',
    'matplotlib',
    'ipywidgets',
    'ipykernel',
    'pandas',
    'tqdm',
    'pyarrow',
    'pytest'
    # ... other dependencies ...
    ],
    # dependencies can be listed under install_requires
)
