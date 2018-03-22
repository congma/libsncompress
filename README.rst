libsncompress
=============

|pipeline| |coverage| |pypi| |license| |pyversions|

Summary
-------

``libsncompress`` – efficient and reproducible Python utility for
compressing supernova cosmological data.

Introduction
------------

The Python package ``libsncompress`` implements the linear compression
method described in the paper “Application of Bayesian graphs to SN Ia
data analysis and compression” (C. Ma, P.-S. Corasaniti, &
B. A. Bassett, 2016, `MNRAS, submitted`_, “M16”; accepted version: 2016
MNRAS, 463, 1651, DOI: `☞`_\ ``10.1093/mnras/stw2069``, BibCode:
`☞ <http://adsabs.harvard.edu/abs/2016MNRAS.463.1651M>`__\ ``2016MNRAS.463.1651M``).
It is designed for use with the `JLA`_ dataset, but can be easily
extended for other similar datasets.

It also includes a Python executable script,
`☞ <https://gitlab.com/congma/libsncompress/blob/master/scripts/jlacompress>`__\ ``jlacompress``,
that serves as a command-line user interface.

The intended usage includes the following tasks:

-  To obtain a compressed sample of the JLA dataset, obtaining the
   location (mean) vector and scatter (covariance) matrix of the
   luminosity distance modulus, possibly for subsequent cosmological
   analysis.
-  To perform cross-validation on portions of the dataset.
-  To study the posterior distribution of SN standardization parameters
   ``alpha``, ``beta``, and ``delta_M``.
-  To replicate the results in the above-cited research paper that is
   built on the preceding tasks.
-  To help developing new methods of data analysis and compression with
   current and future data based on the current work.

The programs work with both Python 2.7 and 3.6.

Installation
------------

The most convenient to install the package depends on your intended
usage.

To simply use the package and compression script, just install the
latest version fro the PyPI with ``pip``:

::

    pip install -U libsncompress

Additional Python packages are required at runtime (please refer to the
section “`Dependencies`_” for the list of dependencies of the current
version). The recommended installation method is to use the above
command, which will make sure that the supporting packages are installed
automatically.

To contribute to the package development, run the tests with real data,
or verify the reproducibility of the research, it is necessary to clone
the repository using ``git``:

::

    git clone https://gitlab.com/congma/libsncompress.git

The full repository includes also the necessary files for testing and
verification. Please refer to the section “`Testing and Development`_”
for further details.

It is also possible to use the library package ``libsncompress`` without
installation, for example, by including them directly in your own
project.

Using the JLA Compression Script
--------------------------------

This utility comes with an executable script ``jlacompress`` that is
tailored to the compression of the `JLA`_ dataset, as done in our `M16`_
paper.

Synopsis
~~~~~~~~

::

    jlacompress [-h] [-d DIR] [-t FILE] [-p PREFIX] [-c z1 z2 [...]] [-n] [-v]

The script requires JLA data files to run. See the section “`Data
Files`_” for details.

Command-Line Options
~~~~~~~~~~~~~~~~~~~~

-  ``-h``, ``--help``: show help message and exit
-  ``-d DIR``, ``--covdir DIR``: path to directory of FITS covariance
   files (default: ``./covmat``)
-  ``-t FILE``, ``--table FILE``: path to JLA data table file (default:
   ``./jla_lcparams.txt``)
-  ``-p PREFIX``, ``--output-prefix PREFIX``: prefix to output file
   names
-  ``-c z1 z2 [...]``, ``--controls z1 z2 [...]``: locations of control
   points (as redshift). At least two control points are required. If
   unspecified, use the default control points in the `JLA paper`_.
-  ``-n``, ``--no-logdet``: don’t use the correct conditional
   probability (default: off; *warning:* use at your own risk)
-  ``-v``, ``--verbose``: turn on verbose output (default: off)

Output
~~~~~~

The script writes three output files when it solves the optimization
problem successfully:

-  Mean (approximate, actually posterior-maximizing) compression
   parameters in the order of (``alpha``, ``beta``, ``delta_M``,
   ``mu1``, ``mu2``, …), where ``muN`` is the value of distance modulus
   at the ``N``-th control point. It is therefore a list of ``N + 3``
   numbers
-  Covariance matrix (symmetric, and with all elements filled) of the
   parameters. It is of shape ``N + 3`` × ``N + 3``, and the first 3
   rows and columns correspond to the three standardization parameters
   ``alpha``, ``beta``, and ``delta_M``.
-  Redshifts of the ``N`` control points.

The default output file names are ``mean.txt``, ``cov.txt`` and
``redshift.txt`` respectively, but they can be prefixed by arbitrary
strings specified by the user with the ``-p``/``--prefix`` option.

The slash (``/``) character in the prefix will be understood as
directory separators. If the directory part of a resulting path does not
exist, it will be created if possible, and nested directories may be
created by this process.

The output files will have suffix ``-no-logdet`` appended to the path
but before the ``.txt`` extension, if ``-n`` or ``--no-logdet`` is
specified.

When verbose output is enabled by ``-v`` or ``--verbose``, additional
text will be written to the standard error.

Exit Status
~~~~~~~~~~~

The script exits with ``0`` for success. Any other value indicates
error.

Example Usage
~~~~~~~~~~~~~

Assuming the data files are in their default locations, the following
command reproduces the default compression results in the `JLA paper`_.

::

    jlacompress -n

Data Files
----------

The JLA data files are *required* for using the package. However, we
cannot distribute them with the source package. Please read the `JLA
readme`_ page for details about the data files.

The following *two* files must be downloaded:

1. The file
   `☞ <http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v6.tgz>`__\ ``jla_likelihood_v6.tgz``,
   compressed archive containing the file ``data/jla_lcparams.txt``.
   This file contains the supernova sample catalogue. The other files in
   this archive are not necessary.
2. The FITS files containing the components of data covariance, in the
   compressed archive
   `☞ <http://supernovae.in2p3.fr/sdss_snls_jla/covmat_v6.tgz>`__\ ``covmat_v6.tgz``.
   The non-FITS files in this archive are not necessary.

If the JLA data archives are already downloaded, you simply need to
extract the required files and specify their locations when using the
``jlacompress`` script, as described `above`_.

The Git repository includes a shell script to download and extract these
files:
`☞ <https://gitlab.com/congma/libsncompress/blob/master/download_jla.sh>`__\ ``download_jla.sh``.
This script is meant to be run manually, and it is not distributed with
the source package on PyPI.

To use the download script, simply invoking the script in the repository
directory

::

    ./download_jla.sh

will suffice – this will populate the ``testdata`` directory with the
necessary files and check the file integrity. Doing so also ensures that
the tests can run.

Testing and Development
-----------------------

Using ``libsncompress`` in Your Project
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use the package directly in your own Python project, simply

.. code:: python

    import libsncompress

This will import three classes from its sub-modules into the
``libsncompress`` namespace:

-  ``BinnedSN``: data-file loader and pre-processor
-  ``BinCollection``: redshift binning and sanitizer; not very useful on
   its own
-  ``CovEvaluator``: the actual compressor

The first thing you need to do is to specify a list (or ``numpy`` array)
of control points, by their *base-10 logarithm* values. Currently, due
to development legacy, the “binning” class and methods are not
particularly efficient. This is usually not a problem because it will be
used only once.

This list or array of control points must be encapsulate in *another*
container (list, array, or tuple, etc.) before passing to the
initializer of ``libsncompress.BinnedSN`` class. The instance can be
initialized by

.. code:: python

    binned_sn = libsncompress.BinnedSN(basedirpath,
                                       tablepath,
                                       logbins=control_points)

Here ``basedirpath`` is the path to the directory containing the FITS
covariance data files, ``tablepath`` the path to the text file
containing the JLA dataset table, and ``logbins`` is the nested list of
control points just obtained.

After this, we can initialize the evaluator
``libsncompress.CovEvaluator`` class, which implements the evaluation of
probability log-density functions and their first 2 derivatives, like
this:

.. code:: python

    ev = libsncompress.CovEvaluator(binned_sn, withlogdet=True)

The optional argument ``withlogdet`` controls whether the full effect of
parameter-dependent covariance matrix is taken into account. It is so
named due to the ubiquitous presence of “ln det Cov” term. It defaults
to ``True`` but can be set to ``False``, which will evaluate the
functions as if the customary chi-squared method were used.

The ``CovEvaluator`` instance, ``ev``, provides a method ``minimize``,
which is a wrapper of ``scipy.optimize.minimize``. Additional positional
and keyword arguments are passed over to that function. The recommended
optimization algorithm is ``trust-ncg`` which fully utilizes the Hessian
matrix. This is the default minimization algorithm if left unspecified,
and other algorithms supported by
`☞ <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html>`__\ ``scipy.optimize.minimize``
can be passed as the optional keyword parameter ``method``.

The return value of ``CovEvaluator.minimize`` method is simply that of
the underlying ``scipy`` function, but with results suitably scaled.

The Hessian of log-PDF function can be obtained, then, at the minimizing
point in the parameter space. This can be used for constructing the
approximate covariance of compression parameters.

Please notice that this implementation here is not a general, abstract
implementation of the linear compression method detailed in `our
paper`_. It specifically implements the sawtooth-basis compression,
which is compatible with the original `JLA one`_. The implementation
details, as well as the exposed API, are likely to see significant
revisions in the future.

Setting Up the Testing Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run the tests (including the reproducibility tests), it is necessary
to set up the environment with supporting packages and data.

As described in the `preceding section`_, “Data Files”, the recommended
way is to clone the Git repository and populate the ``testdata``
directory in the repository with the necessary files, which can be done
using the ``download_jla.sh`` script.

After obtaining the data files, it is recommended to use the recent
version of `☞ <https://tox.readthedocs.io/>`__\ ``tox`` to manage the
testing environments.

::

    pip install 'tox >= 2.8.0'

Although not strictly necessary for running the tests themselves *per
se*, it is recommended to install the
`☞ <http://pandoc.org/>`__\ ``pandoc`` program (please consult your
operating system documentation) and the
`☞ <https://github.com/bebraw/pypandoc>`__\ ``pypandoc`` Python package.

Running the Tests
~~~~~~~~~~~~~~~~~

If you have both Python 2.7 and 3.6 installed, simply invoking

::

    tox

will create the source distribution and run the tests under both Python
variants. The default configuration will pull the latest supporting
packages from PyPI specified in the file ``devel-requirements.txt``.

If you have only one working variant of Python, for example Python 2.7,
you can run

::

    tox -e py2,coverage-report

and skip the unavailable test environment setting.

Reproducibility Tests
~~~~~~~~~~~~~~~~~~~~~

One important goal of the test suits in this repository is to ensure
that the results of JLA SNIa compression are always reproducible.

First, as we have shown in `M16`_, the `JLA`_ compression results (their
Tables F.1 and F.2), especially the covariance matrix, are “very close”
to the ones obtained using this program on the `JLA data release`_, but
with the (highly discouraged) ``withlogdet=False`` option enabled for
``libsncompress.CovEvaluator``.

Second, the compression results produced by this program on the released
JLA data must match those presented in `M16`_, Tables A1 and A2.

The reproducibility tests check that these constraints are satisfied by
all revisions to the codebase. These tests are included in the
``tests/test_reprod.py`` script and are run by ``tox`` by default.

Dependencies
------------

-  `☞ <https://pythonhosted.org/six/>`__\ ``six`` (unknown version), for
   Python 2 and 3 compatibility;
-  `☞ <http://www.numpy.org/>`__\ ``numpy`` (``>= 1.6.0``), for array
   data structure and basic operations;
-  `☞ <https://www.scipy.org/>`__\ ``scipy`` (``>= 0.11.0``), for linear
   algebra and numerical optimization;
-  `☞ <https://www.astropy.org/>`__\ ``astropy`` (unknown version), for
   loading FITS files with the ``astropy.io.fits`` module, which
   replaces the dependence on
   `☞ <https://pythonhosted.org/pyfits/>`__\ ``pyfits`` in earlier
   versions;
-  `☞ <https://pythonhosted.org/cachetools/>`__\ ``cachetools`` (unknown
   version), for caching partial evaluation results, which is essential
   for compression speed.

Performance Notes
-----------------

Performance is mostly determined by the following two conditions:

1. Underlying BLAS/LAPACK libraries used by ``numpy``/``scipy``,
   especially the “linear solver by Cholesky decomposition”,
   ``(D)POTRS`` function of LAPACK. For `NetLib LAPACK`_, this in turn
   is largely determined by the speed of the level-3 BLAS triangular
   solver, ``(D)TRSM``. The NetLib reference implementation is rather
   naive, and an optimized implementation of BLAS is likely to boost the
   performance.
2. Choice of initial value and scaling for numerical optimization. If
   they are suitably chosen, the number of iterations required to
   achieve convergence is reduced.

The script
`☞ <https://gitlab.com/congma/libsncompress/blob/master/scripts/jlacompress>`__\ ``jlacompress``
attempts to automatically create acceptable initial value and scaling
that is optimized for the *default* compression used in the `JLA
paper`_. The automatic initial value and scaling are not optimized for
any other usage cases.

Reporting Bugs
--------------

Please report problems via the `issue tracker`_.

Bibliography
------------

If you use this program in your research, we would like to suggest you
cite the following paper (“M16”):

Ma, C., Corasaniti, P.-S., & Bassett, B. A. 2016, MNRAS, 463, 1651,
`☞`_\ ``doi: 10.1093/mnras/stw2069``

The following BibTeX entry could be useful in a LaTeX document:

::

    @ARTICLE{2016MNRAS.463.1651M,
       author = {{Ma}, C. and {Corasaniti}, P.-S. and {Bassett}, B.~A.},
        title = "{Application of Bayesian graphs to SN Ia data analysis and compression}",
      journal = {MNRAS},
    archivePrefix = "arXiv",
       eprint = {1603.08519},
         year = 2016,
        month = dec,
       volume = 463,
        pages = {1651-1665},
          doi = {10.1093/mnras/stw2069}
    }

.. _MNRAS, submitted: https://arxiv.org/abs/1603.08519
.. _☞: https://doi.org/10.1093/mnras/stw2069
.. _JLA: https://arxiv.org/abs/1401.4064
.. _Dependencies: #dependencies
.. _Testing and Development: #testing-and-development
.. _M16: https://arxiv.org/abs/1603.08519
.. _Data Files: #data-files
.. _JLA paper: https://arxiv.org/abs/1401.4064
.. _JLA readme: http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html
.. _above: #command-line-options
.. _our paper: https://arxiv.org/abs/1603.08519
.. _JLA one: https://arxiv.org/abs/1401.4064
.. _preceding section: #data-files
.. _JLA data release: http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html
.. _NetLib LAPACK: http://www.netlib.org/lapack/
.. _issue tracker: https://gitlab.com/congma/libsncompress/issues

.. |pipeline| image:: https://gitlab.com/congma/libsncompress/badges/master/pipeline.svg
   :target: https://gitlab.com/congma/libsncompress/commits/master
.. |coverage| image:: https://gitlab.com/congma/libsncompress/badges/master/coverage.svg
   :target: https://gitlab.com/congma/libsncompress/commits/master
.. |pypi| image:: https://img.shields.io/pypi/v/libsncompress.svg
   :target: https://pypi.org/project/libsncompress/
.. |license| image:: https://img.shields.io/pypi/l/libsncompress.svg?longCache=true
   :target: https://gitlab.com/congma/libsncompress/blob/master/COPYING
.. |pyversions| image:: https://img.shields.io/pypi/pyversions/libsncompress.svg
   :target: #introduction
