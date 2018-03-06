#!/usr/bin/env python
from distutils.core import setup

pname = "libsncompress"
setup(name=pname, version="0.0.3",
      description="Compress JLA-like supernova data",
      author="Cong Ma",
      author_email="cong.ma@obspm.fr",
      url="https://gitlab.com/congma/libsncompress/",
      packages=[pname],
      scripts=["scripts/jlacompress"],
      requires=["six", "numpy (>=1.6.0)", "scipy (>=0.11.0)", "astropy",
                "cachetools"],
      provides=[pname])
