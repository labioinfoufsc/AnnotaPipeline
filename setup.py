from setuptools import find_packages, setup

requirements = [
    "biopython>=1.73",
    "pandas>=0.24.1",
    "PyYAML>=6.0",
    #'shutil',
    #'argparse',
    #'configparser',
    #'logging',
    #'os',
    #'pathlib',
    #'subprocess',
    #'sys',
    #'re',
    #'warnings'
]

# default python libraries
# https://docs.python.org/3/library/

with open("README.md") as rm:
    long_description = rm.read()

setup(
    name="AnnotaPipeline",
    packages=find_packages(),
    py_modules=["annotation_io", "blast_classification"],
    package_dir={"": "scripts"},
    version="2.0.0",
    description="Pipeline to predict and annotate eukaryotic genes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    author="Guilherme Augusto Maia, Eric Kazuo Kawagoe, Vilmar Benetti Filho, \
            Tatiany Aparecida Soratto, Renato Moreira Simões, Glauber Wagner",
    author_email="labinfo.ufsc@gmail.com",
    zip_safe=False,
    scripts=[
        "scripts/annotapipeline.py",
        "scripts/blastp_parser.py",
        "scripts/fasta_simple.py",
        "scripts/fastatogff.py",
        "scripts/funcannotation_parser.py",
        "scripts/gfftofasta_parser.py",
        "scripts/info_parser.py",
        "scripts/kallisto_parser.py",
        "scripts/percolator_parser.py",
        "scripts/summary_parser.py",
    ],
    include_package_data=True,
    python_requires=">=3.6.9",
    url="https://github.com/labioinfoufsc/AnnotaPipeline.git",
    download_url="https://github.com/labioinfoufsc/AnnotaPipeline/archive/refs/heads/v1.0.zip",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Natural Language :: Portuguese (Brazilian)Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
