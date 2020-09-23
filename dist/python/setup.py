import sys

from setuptools import setup, find_packages
        
def build_native(spec):
    cmd = ["cargo", "build", "--manifest-path", "../../Cargo.toml", "--release"]

    build = spec.add_external_build(cmd=cmd, path=".")

    rtld_flags = ["NOW"]
    if sys.platform == "darwin":
        rtld_flags.append("NODELETE")
        
    spec.add_cffi_module(
        module_path="pcon._lowlevel",
        dylib=lambda: build.find_dylib("pcon", in_path="../../target/release"),
        header_filename=lambda: build.find_header("pcon.h", in_path="../c/"),
        rtld_flags=rtld_flags,
    )

CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Rust",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

CLASSIFIERS.append("Development Status :: 4 - Beta")

SETUP_METADATA = {
    "name": "pcon",
    "description": "lib to perform kmer counting in memory (k <= 31)",
    "url": "https://github.com/natir/pcon",
    "author": "Pierre Marijon",
    "author_email": "marijon@hhu.de",
    "license": "BSD 3-clause",
    "packages": find_packages(),
    "setup_requires": [
                       "setuptools>=38.6.0",
                       "cffi>=1.0.0",
                       "milksnake",
                       ],
    "install_requires": [
                         "cffi>=1.0.0",
                         ],
    "version": "0.1",
    "platforms": "any",
    "classifiers": CLASSIFIERS,
    "milksnake_tasks": [build_native],
}

setup(**SETUP_METADATA)

