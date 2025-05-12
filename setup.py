import setuptools
from glob import glob

# Will load the README.md file into a long_description of the package
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
# Load the requirements file
with open('requirements.txt') as f:
    required = f.read().splitlines()

if __name__ == "__main__":
    setuptools.setup(
        name='chquant',
        version='1.0.0',
        author='Jaiming Chung, Bienfait Isamura',
        author_email='jaiming.chung@manchester.ac.uk, bienfait.isamura@postgrad.manchester.ac.uk',
        description="A package for performing convex hull calculations on molecular dynamics trajectories",
        long_description=long_description,
        long_description_content_type="text/markdown",
        license='MIT',
        python_requires='>=3.7.3',
        install_requires=required,
        zip_safe= False,
        package_dir={"": "src"},
        packages=setuptools.find_packages(where='src'),
        include_package_data=True,
    )
