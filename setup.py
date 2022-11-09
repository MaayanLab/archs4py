import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="archs4py",
    version="0.0.11",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="ARCHS4 python package supporting data loading and data queries.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maayanlab/archs4py",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_data={
        "archs4py": ["data/*"]
    },
    include_package_data=True,
    install_requires=list(map(str.strip, open('requirements.txt', 'r').readlines())),
    python_requires='>=3.6',
)