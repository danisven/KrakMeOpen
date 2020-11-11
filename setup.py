from setuptools import setup, find_packages
from krakmeopen import __version__

setup(
    name="KrakMeOpen",
    version=__version__,
    url="https://github.com/danisven/krakmeopen",
    description="A Kraken2 downstream analysis toolkit.",
    license="MIT",

    # Author details
    author='Daniel Svensson',
    author_email='daniel.svensson@umu.se',

    keywords="Bioinformatics NGS kraken2 Metagenomics",
    classifiers=[
        'Development Status :: 5 - Beta',
        'License :: OSI Approved :: MIT',
        'Programming Language :: Python :: 3'],
    include_package_data=True,
    install_requires=['pyyaml', 'pandas'],
    python_requires='>=3.7',
    packages=find_packages(exclude=['contrib', 'docs', 'test*'], include=['krakmeopen']),
    entry_points={'console_scripts': ['krakmeopen=krakmeopen.krakmeopen:krakmeopen']})
