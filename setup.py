from setuptools import setup, find_packages

setup(
    name="mbd-framework",
    version="1.0.0",
    description="A formal validation toolkit calculating Many-Body Dispersion bounds connecting geometric theorems derived in Lean.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/edqa/MBD-Theoretical-Framework",
    author="Edwin Maina",
    license="MIT",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pyscf>=2.0.0",
    ],
    entry_points={
        "console_scripts": [
            "mbd-compute=mbd_framework.compute_volumes:main",
            "mbd-crystal=mbd_framework.crystal_validation:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
)
