import os
import re
import codecs
from setuptools import setup, find_packages

NAME = 'cmonkey2'
PACKAGES = ['cmonkey', 'cmviewer', 'config']
DESCRIPTION = 'cmonkey2 is an implementation of the cmonkey biclustering method in Python'
LICENSE = 'LGPL V3'
URI = 'https://github.com/baliga-lab/cmonkey2'
AUTHOR = 'Baliga Lab, Institute for Systems Biology'
VERSION = '1.1.0'

KEYWORDS = ['class', 'cmonkey2']
CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Scientists",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU Lesser Public License",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Software Development :: Libraries :: Python Modules"
    ]
INSTALL_REQUIRES = ['numpy', 'scipy',
                        'rpy2', 'beautifulsoup4', 'configparser', 'pandas', 'biopython', 'lxml',
                        'cherrypy', 'jinja2', 'routes']

PACKAGE_DATA = {
    '': ['*.json', '*.ini', '*.fa', 'KEGG_taxonomy', 'proteome2taxid'],
    'cmviewer': ['static/images/*',
                 'static/javascripts/*', 'static/javascripts/modules/*',
                 'static/javascripts/jquery-ui-1.11.2.custom/*',
                 'static/javascripts/jquery-ui-1.11.2.custom/external/jquery/*',
                 'static/javascripts/jquery-ui-1.11.2.custom/images/*',
                 'static/stylesheets/*',
                 'templates/*']
    }

if __name__ == '__main__':
    setup(name=NAME, description=DESCRIPTION,
          license=LICENSE,
          url=URI,
          version=VERSION,
          author=AUTHOR,
          author_email='wwu@systemsbiology.net',
          maintainer=AUTHOR,
          maintainer_email='wwu@systemsbiology.net',
          keywords=KEYWORDS,
          packages=PACKAGES,
          zip_safe=False,
          classifiers=CLASSIFIERS,
          install_requires=INSTALL_REQUIRES,
          include_package_data=True, package_data=PACKAGE_DATA,
          scripts=['bin/cmonkey2.py', 'bin/cm2view.py'])
