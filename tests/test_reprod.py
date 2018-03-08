"""Testing whether libsncompress creates numerically reproducible outcome.

Comparisons are made against the output of B14 that should be reproduced by
setting ``withlogdet=False``, and that of M16.
"""


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
                      warnings.ResourceWarning)
        # Try to assume we're in the repo directory and guess basedir
        basedir = "testdata"
    return basedir


@pytest.fixture(scope="session")
def jla_target(vname="jla_mub.txt", covname="jla_mub_covmatrix.dat"):
    """Returns the JLA (B14) compression redshift, mean vector, and covariance
    matrix.
    """
    vpath, covpath = search_for(data_basedir(), vname, covname)
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
    vpath, covpath = search_for(data_basedir(), vname, covname)
    varray = numpy.loadtxt(vpath)
    z = varray[:, 0]
    v = varray[:, 1]
    cov = numpy.loadtxt(covpath)
    return z, v, cov


@pytest.fixture(scope="session")
def jla_full_paths():
    """Returns the JLA data file paths needed by the computation.
    """
    fits_cov_path, = search_for(data_basedir(), "covmat", fordirs=True)
    table_path, = search_for(data_basedir(), "jla_lcparams.txt")
    return fits_cov_path, table_path


def test_reprod_jla(jla_target, jla_full_paths):
    """Test if we reproduce JLA's compression results when given
    withlogdet=False.
    """
    their_z, their_v, their_cov = jla_target
    fitsdir, tablepath = jla_full_paths
    # NOTE: Don't use their z list, which is truncated.
    # They MEANT that their z's should be equidistant in log-z space
    # and this is the default.
    binned_sn = libsncompress.BinnedSN(fitsdir, tablepath)
    sinned_ev = libsncompress.CovEvaluator(binned_sn, withlogdet=False)
    minres = sinned_ev.minimize(method="trust-ncg", options={})
    assert minres.success   # Must converge
    hess = sinned_ev.chisqhess(minres.x)
    hessfac = cho_factor(hess, check_finite=False)
    cov = lmultiply_inv(hessfac, numpy.eye(sinned_ev.dimension))
    cov = cov + cov.T
    our_v = minres.x[3:]
    our_cov = cov[3:, 3:]
    kld = kldgaussian(their_v, their_cov, our_v, our_cov)
    assert 0.0 <= kld <= 0.0001     # Difference in nats
    # Notice that B14's convergence is not well-understood, and directly
    # comparing the optimizing vectors may not be meaningful.


def test_reprod_m16(m16_target, jla_full_paths):
    """Test if we reproduce the tables in M16."""
    their_z, their_v, their_cov = m16_target
    fitsdir, tablepath = jla_full_paths
    binned_sn = libsncompress.BinnedSN(fitsdir, tablepath,
                                       [numpy.log10(their_z)])
    ev = libsncompress.CovEvaluator(binned_sn)
    minres = ev.minimize(method="trust-ncg", options={})
    assert minres.success
    hess = ev.chisqhess(minres.x)
    hessfac = cho_factor(hess, check_finite=False)
    cov = lmultiply_inv(hessfac, numpy.eye(ev.dimension))
    cov = cov + cov.T
    # Check in the truncated parameter space of only distances, without
    # the corresponding standardization parameters (the first 3 ones).
    # Reason: in M16 the exact values are only available for the compressed
    # distance vectors as electronic table (Table A1).
    # In this case we need element-level almost exact match for this subspace.
    assert numpy.allclose(their_v, minres.x[3:])
    # Notice that the full-dimension covariance table is available in exact
    # values (Table A2).
    assert numpy.allclose(their_cov, cov)
    # Check that the full optimization result matches.  This is done by
    # concatenating the target vector array with the three values from Table 3
    # of the M16 publication on the head.
    their_full_v = numpy.concatenate(((0.125, -2.58, -0.052), their_v))
    kld = kldgaussian(their_full_v, their_cov, minres.x, cov)
    # In this case, because the first 3 target values are not exact, we check
    # the KL divergence in nats.
    assert 0.0 <= kld <= 0.004
