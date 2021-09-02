from setuptools import setup, Extension

# Getting Description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# C Extension
cutilsmodule = Extension('cutils',
                    sources = [ 'aleph/cutils.c',],
                    include_dirs = ['aleph'],
                    #define_macros=[ ('CUTILS', None) ],
                    # Removed '-march=native' for now.
                    extra_compile_args=['-O3','-Wno-unknown-pragmas', '-fPIC', '-shared'], # '-DCUTILS', ghash_arg, '-fstrict-aliasing', '-std=c99', '-D_GNU_SOURCE'
                    #extra_link_args=extra_link_args,
                    language = 'c',
                    )

setup(
      name="aleph",
      version="0.0.1",
      author="Jose Pena Z.",
      author_email="jpena@das.uchile.cl",
      description="ALeRCE's Ephemeris",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/josepenaz/ephs",
      #project_urls={
      #    "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
      #},
      #classifiers=[
      #    "Programming Language :: Python :: 3",
      #    "License :: OSI Approved :: MIT License",
      #    "Operating System :: OS Independent",
      #],
      package_dir={"aleph": "aleph"},
      packages=['aleph'],
      python_requires=">=3.6",
      install_requires=['numpy>=1.18.1',
                        'astropy>=4.2',
                        'pandas>=1.1.3',
                        'rebound>=3.12.2',
                        ],
      package_data = {'aleph':['cutils.h']},
      ext_modules = [cutilsmodule],

)
