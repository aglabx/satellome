import re

from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = None
for line in open("./src/satellome/__init__.py"):
    m = re.search("__version__\s*=\s*(.*)", line)
    if m:
        version = m.group(1).strip()[1:-1]  # quotes
        break
assert version

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
        "editdistance",
        "tqdm",
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
)
