from setuptools import setup, find_packages

setup(
    name="pyAir",
    version="0.2.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_packages_data=True,
)
