#!/usr/bin/env python
import sys
import os
import os.path
from glob import glob
import errno
from setuptools import setup, find_packages
import six


if six.PY3:
    def open_utf8(path, mode):
        return open(path, mode, encoding="utf-8")
else:
    open_utf8 = open


def get_rst_text_from_md(path):
    with open_utf8(path, "r") as f:
        orig_text = f.read()
    try:
        import pypandoc
    except (OSError, ImportError):
        # This is the case when pypandoc isn't available.
        return None
    # Use PyPandoc to convert the markdown, and overwrite the .rst file
    return pypandoc.convert_text(orig_text, "rst", "md",
                                 extra_args=["--reference-links",
                                             "--strip-comments"])


def should_convert_p(from_path, to_path):
    try:
        stat_to = os.stat(to_path)
    except OSError as err:
        if err.errno == errno.ENOENT:
            # "to_path" does not exist, therefore we should convert.
            return True
        # Other kinds of OSError are definitely an error.
        six.raise_from(IOError, err)
    stat_from = os.stat(from_path)
    # If "from_path" is more recent than "to_path", we should convert.
    return stat_from.st_mtime > stat_to.st_mtime


from_doc = "README.md"
to_doc = "README.rst"
if should_convert_p(from_doc, to_doc):
    rst_text = get_rst_text_from_md(from_doc)
    if rst_text is not None:
        if six.PY2:
            rst_text = rst_text.encode("utf-8")
        with open(to_doc, "w") as f:
            f.write(rst_text)
    else:
        sys.stderr.write("WARNING: Automatic README.rst generation failed.\n"
                         "         Using placeholder (not written to file).\n")
        rst_text = ("Please refer to the documentation at the `homepage`_.\n"
                    ".. _homepage: https://gitlab.com/congma/libsncompress/\n"
                    "(You are seeing this because automatic generation of "
                    "``README.rst`` failed.)\n")
else:
    with open_utf8(to_doc, "r") as f:
        rst_text = f.read()


pname = "libsncompress"
setup(name=pname, version="0.0.8",
      description="Compress JLA-like supernova data",
      long_description=rst_text,
      author="Cong Ma",
      author_email="cong.ma@obspm.fr",
      url="https://gitlab.com/congma/libsncompress/",
      packages=find_packages("src"),
      package_dir={"": "src"},
      scripts=["scripts/jlacompress"],
      data_files=[("share/libsncompress", glob("testdata/m16/table_*.txt"))],
      python_requires=", ".join([">=2.7",
                                 "!=3.0.*",
                                 "!=3.1.*",
                                 "!=3.2.*",
                                 "!=3.3.*",
                                 "!=3.4.*",
                                 "!=3.5.*",
                                 "<3.8"]),
      install_requires=["six", "numpy >= 1.6.0", "scipy >= 0.11.0", "astropy",
                        "cachetools"],
      setup_requires=["six", "pypandoc"],
      tests_require=["numpy >= 1.14.0", "scipy >= 1.0.0", "coverage >= 4.2",
                     "pytest >= 3.2.0", "tox >= 2.8.0"],
      provides=[pname],
      license="BSD",
      classifiers=["Development Status :: 4 - Beta",
                   "Environment :: Console",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: BSD License",
                   "Operating System :: OS Independent",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3.6",
                   "Programming Language :: Python :: 3.7",
                   "Programming Language :: Python :: Implementation :: "
                   "CPython",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   "Topic :: Scientific/Engineering :: Information Analysis",
                   "Topic :: Software Development :: Libraries :: "
                   "Python Modules"])
