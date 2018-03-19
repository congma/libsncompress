"""Test infrastructure: setup, utility functions, fixtures."""
import os
import os.path
from collections import defaultdict
import warnings
import pytest
import six
import numpy
from scipy.linalg import cho_factor
import libsncompress
from libsncompress.evaluator import (lmultiply_inv,
                                     quadratic_inv,
                                     logdet_cholesky)


def kldgaussian(from_mean, from_cov, to_mean, to_cov):
    """Kullback--Leibler divergence of two Gaussian variables of the same
    dimension.
    """
    d = from_mean.shape[0]
    vdiff = to_mean - from_mean
    from_covfac, to_covfac = map(cho_factor, (from_cov, to_cov))
    kl = logdet_cholesky(to_covfac) - logdet_cholesky(from_covfac)
    kl += numpy.trace(lmultiply_inv(to_covfac, from_cov)) - d
    kl += quadratic_inv(to_covfac, vdiff, vdiff)
    kl *= 0.5
    return kl


def search_for(directory, *basenames, **kwargs):
    """Search for file names in ``basename`` descending from ``directory``.
    The name must be an exact match.  Returns a list of arbitrary matches,
    where ``None`` indicates a non-matching element.

    If "fordirs=True" is specified in the keyword arguments, look for the names
    in directory names instead. (Default: look up in plain file names).
    """
    # Return early on trivial case.
    fordirs = kwargs.get("fordirs", False)
    if not basenames:
        return []
    uniqnames = defaultdict(list)
    for idx, name in enumerate(basenames):
        uniqnames[name].append(idx)
    n = len(uniqnames)
    matches = {}
    long_matches = [None] * len(basenames)
    for root, subdirs, files in os.walk(os.path.abspath(directory)):
        # Skip as much as we can in case of no base files under root.
        target = subdirs if fordirs else files
        if not target:
            continue  # with next directory-walk step
        # No need to shrink the candidate set as we match; we only need
        # membership test which is O(1).
        target = frozenset(target)
        for name in uniqnames:
            if name in target:
                matches[name] = os.path.join(root, name)
        if len(matches) >= n:
            break  # from the directory-walker early, as we've done our job.
    for name, path in six.iteritems(matches):
        for long_idx in uniqnames[name]:
            long_matches[long_idx] = path
    return long_matches


def data_basedir(envkey="LSNZ_TESTDATA_BASE"):
    basedir = os.getenv(envkey)
    if not basedir or basedir is None:
        # Environment not properly set.
        warnings.warn("Test data base-directory not set "
                      "by environment variable '%s'" % envkey,
                      RuntimeWarning)
        # Try to assume we're in the repo directory and guess basedir
        basedir = "testdata"
    return basedir


TESTDATA_ROOT = data_basedir()


@pytest.fixture(scope="session")
def jla_target(vname="jla_mub.txt", covname="jla_mub_covmatrix.dat"):
    """Returns the JLA (B14) compression redshift, mean vector, and covariance
    matrix.
    """
    vpath, covpath = search_for(TESTDATA_ROOT, vname, covname)
    varray = numpy.loadtxt(vpath)
    z = varray[:, 0]
    v = varray[:, 1]
    carray = numpy.loadtxt(covpath, skiprows=1)
    cov = carray.reshape(len(z), -1)
    return z, v, cov


@pytest.fixture(scope="session")
def m16_target(vname="table_A1.txt", covname="table_A2.txt"):
    """Returns the M16 compression redshift, mean vector, and covariance
    matrix.
    """
    vpath, covpath = search_for(TESTDATA_ROOT, vname, covname)
    varray = numpy.loadtxt(vpath)
    z = varray[:, 0]
    v = varray[:, 1]
    cov = numpy.loadtxt(covpath)
    return z, v, cov


@pytest.fixture(scope="session")
def jla_full_paths():
    """Returns the JLA data file paths needed by the computation.
    """
    fits_cov_path, = search_for(TESTDATA_ROOT, "covmat", fordirs=True)
    table_path, = search_for(TESTDATA_ROOT, "jla_lcparams.txt")
    return fits_cov_path, table_path


@pytest.fixture(scope="session")
def outdir(tmpdir_factory):
    out_dir = tmpdir_factory.mktemp("lsnz_test_output")
    yield out_dir
    out_dir.remove(rec=True)


@pytest.fixture(scope="session")
def ref_binned_sn(jla_full_paths):
    binned_sn = libsncompress.BinnedSN(*jla_full_paths)
    return binned_sn


@pytest.fixture(scope="session")
def ref_ev(ref_binned_sn):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    ev.minimize()
    assert ev.res.success
    return ev


@pytest.fixture(scope="session")
def ref_ev_nologdet(ref_binned_sn):
    ev = libsncompress.CovEvaluator(ref_binned_sn, withlogdet=False)
    ev.minimize()
    assert ev.res.success
    return ev
