import setuptools
import codecs

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

exec(codecs.open("MisER/version.py").read())


INSTALL_REQUIRES = [
    'numpy',
    'pandas',
    'pysam',
    'parasail>=1.1.17'
]

setuptools.setup(
    name='MisER',
    version=VERSION,
    author='Zhen Liu',
    author_email='liuzhen_sirius@outlook.com',
    description='Find and fix missed small exons.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/zhenLiuXplr/MisER-project',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    package_dir={"MisER": "MisER"},
    install_requires=INSTALL_REQUIRES,
    python_requires='>=3',
    entry_points={
        'console_scripts': ['MisER=MisER.run_miser:main',],
    }
)
