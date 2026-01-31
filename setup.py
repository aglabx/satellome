import sys
import os

from setuptools import find_packages, setup

# Import custom install command
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))
try:
    from satellome.install_command import PostInstallCommand
    cmdclass = {'install': PostInstallCommand}
except ImportError:
    # Fallback if import fails (e.g., during package build)
    cmdclass = {}

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read version from pyproject.toml (single source of truth)
version = None
with open("pyproject.toml", "r") as f:
    for line in f:
        if line.strip().startswith("version"):
            version = line.split("=")[1].strip().strip('"')
            break
assert version, "Could not read version from pyproject.toml"

setup(
    name="satellome",
    version=version,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={"": ["README.md"]},
    python_requires=">=3.6",
    include_package_data=True,
    scripts=[],
    license="MIT",
    url="https://github.com/aglabx/Satellome",
    author="Aleksey Komissarov",
    author_email="ad3002@gmail.com",
    description="Satellome: a tool for satellite DNA analysis in T2T assemblies.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "PyExp",
        "numpy",
        "pyyaml",
        "simplejson",
        "intervaltree",
        "tqdm",
        "plotly",
        "editdistance",
        "networkx",
        "matplotlib",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
        'console_scripts': [
            'satellome = satellome.main:main',
        ],
    },
    cmdclass=cmdclass,
)
