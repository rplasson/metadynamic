import setuptools
import os
import sys

from glob import glob

_here = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(_here, "README.md")) as f:
        long_description = f.read()
else:
    with open(os.path.join(_here, "README.md"), encoding="utf-8") as f:
        long_description = f.read()

version = {}
with open(os.path.join(_here, "metadynamic", "version.py")) as f:
    exec(f.read(), version)

setuptools.setup(
    name="metadynamic",
    version=version["__version__"],
    description=("Metadynamic Gillespie simulations."),
    long_description=long_description,
    author="Raphaël Plasson",
    author_email="raphael.plasson@univ-avignon.fr",
    # url='https://github.com/bast/somepackage',
    license="GPL V3.0",
    packages=setuptools.find_packages(),
    # no dependencies in this example
    install_requires=["setuptools>=33.1.1", "numpy>=1.12.1"],
    # install_requires=[
    #       'dependency==1.2.3',
    # ],
    # no scripts in this example
    scripts=["bin/metarun"],
    include_package_data=True,
    data_files=[
        ("share/man/man1", ["docs/man/metarun.1"]),
        ("share/doc/python3-metadynamic/apidocs", glob("docs/apidocs/*")),
        ("share/doc/python3-metadynamic/examples", glob("docs/examples/*")),
        ("share/doc/python3-metadynamic/tests", glob("docs/tests/*")),
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
