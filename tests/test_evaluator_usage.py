"""Testing usage of libsncompress.evaluator"""
import pytest
import numpy
import libsncompress
from lsnz_test_infra import jla_full_paths, ref_binned_sn, ref_ev, kldgaussian


def veryclose_p(ref, alt, tol=0.0001, demand_positive=False):
    kld = kldgaussian(ref.res.x, ref.compressed_cov,
                      alt.res.x, alt.compressed_cov)
    truth = kld <= tol
    if demand_positive:
        return truth and kld >= 0.0
    return truth and (numpy.allclose(kld, 0.0) if kld < 0.0 else True)


def test_evaluator_compressed_ev(ref_ev):
    assert ref_ev.res.compressed_cov is ref_ev.compressed_cov


def test_unready_evaluator(ref_binned_sn):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    assert ev.res is None
    with pytest.raises(ValueError):
        _ = ev.compressed_cov


def test_min_method_newton_cg(ref_binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    ev.minimize(method="Newton-CG", options=dict(xtol=1e-7))
    assert ev.res.success
    assert veryclose_p(ref_ev, ev)


@pytest.mark.filterwarnings("ignore:Method bfgs does not use Hessian")
def test_min_method_bfgs(ref_binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    # BFGS does not use Hessian at all, we can even pass in garbage
    ev.minimize(method="bfgs", hess=lambda x: None, options=dict(gtol=1e-6))
    assert ev.res.success
    assert veryclose_p(ref_ev, ev)


@pytest.mark.filterwarnings("ignore:Method Powell does not use Hessian")
@pytest.mark.filterwarnings("ignore:Method Powell does not use gradient")
def test_min_method_powell(ref_binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    ev.minimize(method="Powell")
    assert ev.res.success
    # Do no check for closeness -- Powell's method has very poor convergence.
    # This is just to check that the routines can go without the Jacobian and
    # not break.


@pytest.mark.parametrize("method", ["dogleg", "trust-exact"])
def test_min_method_trust_regions(ref_binned_sn, ref_ev, method):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    ev.minimize(method=method)
    assert ev.res.success
    assert veryclose_p(ref_ev, ev)


def test_min_alter_x0(ref_binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    x0 = numpy.zeros(3 + ref_binned_sn.bins.ncontrolpoints)
    x0[3:] = 10.0
    ev.minimize(x0=x0)
    assert ev.res.success
    assert veryclose_p(ref_ev, ev)


def test_min_alter_scalings(ref_binned_sn, ref_ev):
    ev = libsncompress.CovEvaluator(ref_binned_sn)
    xs = numpy.ones(3 + ref_binned_sn.bins.ncontrolpoints)
    ev.minimize(xscalings=xs)
    assert ev.res.success
    assert veryclose_p(ref_ev, ev)


def test_min_sorted(jla_full_paths, ref_ev):
    binned_sn_sorted = libsncompress.BinnedSN(*jla_full_paths, sort_by_z=True)
    ev = libsncompress.CovEvaluator(binned_sn_sorted)
    ev.minimize()
    assert ev.res.success
    assert numpy.allclose(ref_ev.res.x, ev.res.x)
    assert numpy.allclose(ref_ev.compressed_cov, ev.compressed_cov)
