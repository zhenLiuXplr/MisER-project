import setuptools
import codecs

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

exec(codecs.open("MisSER/version.py").read())


INSTALL_REQUIRES = [
    'numpy',
    'pandas',
    'pysam',
    'parasail>=1.1.17'
]

setuptools.setup(
    name='MisSER',
    version=VERSION,
    author='Zhen Liu',
    author_email='liuzhen2018@sibs.ac.cn',
    description='Find and fix missed small exons.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/zhenLiuExplr/MisSER-project',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    package_dir={"MisSER": "MisSER"},
    install_requires=INSTALL_REQUIRES,
    python_requires='>=3',
    entry_points={
        'console_scripts': ['MisSER=MisSER.run_misser:main',],
    }
)
