from setuptools import setup, Extension
import re

description = 'A python library to watch over the Cambridge STARS code.'

try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = description

metadata = {"version": "",
            "author": "",
            "email": ""
           }

metadata_file = open("kaitiaki/_metadata.py", "rt").read()

for item in metadata.keys():
    version_regex = rf"^__{item}__ = ['\"]([^'\"]*)['\"]"

    match = re.search(version_regex, metadata_file, re.M)

    if match:
        metadata[item] = match.group(1)

setup(name='kaitiaki',
      license           = 'MIT License',
      version           = metadata['version'],
      description       = description,
      long_description  = long_description,
      author            = metadata['author'],
      author_email      = metadata['email'],
      packages          = ['kaitiaki'],
      zip_safe          = False,
      homepage          = 'https://github.com/Krytic/kaitiaki',
      install_requires  = ['numpy',
                           'matplotlib',
                           'pandas',
                          ]
    )
