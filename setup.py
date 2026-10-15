from setuptools import setup

setup(
    cffi_modules=["aleph/compile_cutils.py:ffibuilder"],
)