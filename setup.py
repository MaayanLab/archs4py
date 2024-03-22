import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="archs4py",
    version="0.2.13",
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
    install_requires=[
        'h5py==3.8.0',
        'numpy==1.24.3',
        'pandas==2.0.2',
        'qnorm==0.8.1',
        'setuptools==67.8.0',
        'tqdm==4.65.0',
        'wget==3.2',
        's3fs==2023.12.2',
        'biomart==0.9.2',
        'xalign==0.1.74'

    ],
    python_requires='>=3.7',
)