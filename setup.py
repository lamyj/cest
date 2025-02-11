import glob
import os
import subprocess
import sys
import tempfile

import setuptools
import setuptools.command.build

here = os.path.abspath(os.path.dirname(__file__))

class BuildCMake(setuptools.Command, setuptools.command.build.SubCommand):
    def __init__(self, *args, **kwargs):
        setuptools.Command.__init__(self, *args, **kwargs)
        setuptools.command.build.SubCommand.__init__(self, *args, **kwargs)
        self.build_lib = None
        self.editable_mode = False
        
        self.sources = []
        self.stanc_options = ""
        self.cxxflags = ""
    
    def initialize_options(self):
        pass
        
    def finalize_options(self):
        self.sources = [*sorted(glob.glob(f"cest/*.cpp"))]
        self.set_undefined_options("build_py", ("build_lib", "build_lib"))
    
    def run(self):
        with tempfile.TemporaryDirectory() as build_dir:
            subprocess.check_call(
                [
                    "cmake", f"-DPython_EXECUTABLE={sys.executable}",
                    "-DCMAKE_BUILD_TYPE=Release",
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY="
                        f"{os.path.join(here, self.build_lib, 'cest')}", 
                    "-S", here, "-B", build_dir])
            
            subprocess.check_call(
                [
                    "cmake", "--build", build_dir,
                    "--target", "pycest",
                    "--config", "Release", "--parallel"])
    
    def get_source_files(self):
        return self.sources

setuptools.command.build.build.sub_commands.append(("build_cmake", None))

long_description = open(os.path.join(here, "README.md")).read()
setuptools.setup(
    name="cest",
    version="0.1.0",
    
    description="CEST Toolbox",
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    url="https://github.com/lamyj/cest",
    
    author="Julien Lamy",
    author_email="lamy@unistra.fr",
    
    license="MIT",
    
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    
    keywords = ["cest", "mri"],
    
    cmdclass={"build_cmake": BuildCMake},

    packages=setuptools.find_packages(where="."),
    package_dir={"": "."},
    
    install_requires=["nibabel", "numpy", "spire"],
)
