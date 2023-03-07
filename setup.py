from setuptools import setup
from os.path import join, dirname

setup(
    name="diversity_tools",
    version="1.0",
    packages=["diversity_tools"],
    long_description=open(join(dirname(__file__), "README.MD")).read(),
)