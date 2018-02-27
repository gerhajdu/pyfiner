from distutils.core import setup

setup(
    name='PyFiNeR',
    version='1.0',
    author='Gergely Hajdu',
    author_email='hajdu.gergely86@gmail.com',
    packages=['pyfiner'],
    #scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE',
    description='Fitting Near-infrared RR Lyrae light curves',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.12",
        "scipy >= 1.0",
        "matplotlib >= 2.1",
    ],
)

