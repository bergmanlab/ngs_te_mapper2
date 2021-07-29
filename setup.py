import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ngs_te_mapper2",
    version="1.0.2",
    author="Shunhua Han",
    author_email="hanshunhua0829@gmail.com",
    description="A program to identify transposable element insertions using next generation sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bergmanlab/ngs_te_mapper2",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    entry_points={"console_scripts": ["ngs_te_mapper2 = ngs_te_mapper2.ngs_te_mapper2:main"]},
    classifiers=(
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
    ),
)
