"""Testing usage of libsncompress.evaluator"""
import pytest
import numpy
import libsncompress
from lsnz_test_infra import jla_full_paths, kldgaussian


@pytest.fixture(scope="module")
def binned_sn(jla_full_paths):
    binned_sn = libsncompress.BinnedSN(*jla_full_paths)
    return binned_sn


@pytest.fixture(scope="module")
def ref_ev(binned_sn):
    ref = libsncompress.CovEvaluator(binned_sn)
    ref.minimize()
    assert ref.res.success
    return ref


def isveryclose_p(ref, alt, tol=0.0001, check_positive=False):
    kld = kldgaussian(ref.res.x, ref.compressed_cov,
                      alt.res.x, alt.compressed_cov)
    return (kld <= tol) and (kld >= 0.0 if check_positive else True)


def test_min_method_newton_cg(binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(binned_sn)
    ev.minimize(method="Newton-CG", options=dict(xtol=1e-7))
    assert ev.res.success
    assert isveryclose_p(ref_ev, ev)


@pytest.mark.filterwarnings("ignore:Method bfgs does not use Hessian")
def test_min_method_bfgs(binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(binned_sn)
    # BFGS does not use Hessian at all, we can even pass in garbage
    ev.minimize(method="bfgs", hess=lambda x: None, options=dict(gtol=1e-6))
    assert ev.res.success
    assert isveryclose_p(ref_ev, ev)


def test_min_method_trust_regions(binned_sn, ref_ev):
    for method in ("dogleg", "trust-exact"):
        ev = libsncompress.CovEvaluator(binned_sn)
        ev.minimize(method=method)
        assert ev.res.success
        assert isveryclose_p(ref_ev, ev)


def test_min_alter_x0(binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(binned_sn)
    x0 = numpy.zeros(3 + binned_sn.bins.ncontrolpoints)
    x0[3:] = 10.0
    ev.minimize(x0=x0)
    assert ev.res.success
    assert isveryclose_p(ref_ev, ev)


def test_min_alter_scalings(binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(binned_sn)
    xs = numpy.ones(3 + binned_sn.bins.ncontrolpoints)
    ev.minimize(xscalings=xs)
    assert ev.res.success
    assert isveryclose_p(ref_ev, ev)
