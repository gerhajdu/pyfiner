from distutils.core import setup

setup(
    name='pyfiner',
    version='1.0',
    author='Gergely Hajdu',
    author_email='hajdu.gergely86@gmail.com',
    packages=['pyfiner'],
    scripts=['scripts/pyfiner.py','scripts/pyfiner'],
    url='https://github.com/gerhajdu/pyfiner',
    license='LICENSE',
    description='Fitting Near-infrared RR Lyrae light curves',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.12",
        "scipy >= 1.0",
        "matplotlib >= 2.1",
    ],
)

