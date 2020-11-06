from setuptools import setup, find_packages
from krakmeopen.krakmeopen import __version__

setup(
    name="KrakMeOpen",
    version=__version__,
    url="https://github.com/danisven/krakmeopen",
    description="A Kraken2 downstream analysis toolkit.",
    license="MIT",

    # Author details
    author='Daniel Svensson',
    author_email='daniel.svensson@umu.se',

    keywords="Bioinformatics NGS kraken2",
    classifiers=[
        'Development Status :: 5 - Beta',
        'License :: OSI Approved :: MIT',
        'Programming Language :: Python :: 3'
        ],
    install_requires=['dataclasses'],
    packages=find_packages(exclude=['contrib', 'docs', 'test*'], include=['stringmeup']),
    entry_points={'console_scripts': [  'stringmeup=stringmeup.stringmeup:stringmeup',
#                                        'kraken2-taxonomy=kraken2_confidence_recal.taxonomy:main',

]})
