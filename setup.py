# This file is part of ENSO.
#
# Copyright (C) 2020 Sebastian Ehlert
#
# ENSO is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ENSO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ENSO. If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="enso",
    version="2.0.1",
    description="Energetic sorting of CREST CRE for automated NMR calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/grimme-lab/enso",
    license="LGPL3",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: LGPL3",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    package_dir={"": "."},
    py_modules=["enso"],
    entry_points={"console_scripts": ["enso=enso:main"]},
    install_requires=[
        "argparse",
    ],
)
