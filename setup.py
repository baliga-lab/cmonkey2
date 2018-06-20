import os
import re
import codecs
from setuptools import setup

NAME = 'cmonkey2'
PACKAGES = ['cmonkey', 'cmonkey.cmviewer', 'cmonkey.tools', 'cmonkey.meme']
DESCRIPTION = 'cmonkey2 is an implementation of the cmonkey biclustering method in Python'
LICENSE = 'LGPL V3'
URI = 'https://github.com/baliga-lab/cmonkey2'
AUTHOR = 'Baliga Lab, Institute for Systems Biology'
VERSION = '1.2.9'

KEYWORDS = ['class', 'cmonkey2']

# See trove classifiers
# https://testpypi.python.org/pypi?%3Aaction=list_classifiers

CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Libraries :: Python Modules"
    ]
INSTALL_REQUIRES = ['numpy', 'scipy',
                    'rpy2<=2.8.6', 'beautifulsoup4', 'configparser', 'pandas', 'biopython', 'lxml',
                    'cherrypy', 'cherrypy_cors', 'jinja2', 'routes', 'svgwrite', 'sqlalchemy', 'sqlalchemy-utils']

PACKAGE_DATA = {
    'cmonkey': ['default_config/*.json', 'default_config/*.ini', 'default_config/*.fa',
                'default_config/KEGG_taxonomy', 'default_config/proteome2taxid'],
    'cmonkey.cmviewer': ['static/images/*',
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
          scripts=['bin/cmonkey2', 'bin/cm2view', 'bin/cm2plot', 'bin/cm2export'])
