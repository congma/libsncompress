"""Testing whether libsncompress creates numerically reproducible outcome.

Comparisons are made against the output of B14 that should be reproduced by
setting ``withlogdet=False``, and that of M16.
"""


import numpy
import libsncompress
from lsnz_test_infra import (jla_target, m16_target,
                             jla_full_paths,
                             ref_binned_sn, ref_ev, ref_ev_nologdet,
                             kldgaussian)


def test_reprod_jla(jla_target, ref_ev_nologdet):
    """Test if we reproduce JLA's compression results as they appeared in the
    B14 paper with our default evaluator given withlogdet=False.
    """
    their_z, their_v, their_cov = jla_target
    # NOTE: Don't use their z list, which is truncated.
    # They MEANT that their z's should be equidistant in log-z space
    # and this is the default.
    # We use the default evaluator object instead, already present as fixture
    assert ref_ev_nologdet.res.success   # Fixture is useable.
    our_v = ref_ev_nologdet.res.x[3:]
    our_cov = ref_ev_nologdet.compressed_cov[3:, 3:]
    kld = kldgaussian(their_v, their_cov, our_v, our_cov)
    assert 0.0 <= kld <= 0.0001     # Difference in nats
    # Notice that B14's convergence is not well-understood, and directly
    # comparing the optimizing vectors may not be meaningful.


def test_reprod_m16(m16_target, ref_ev):
    """Test if our default evaluator reproduces the tables in M16 in an
    element-wise compatible way.
    """
    their_z, their_v, their_cov = m16_target
    assert ref_ev.res.success
    cov = ref_ev.compressed_cov
    # Check in the truncated parameter space of only distances, without
    # the corresponding standardization parameters (the first 3 ones).
    # Reason: in M16 the exact values are only available for the compressed
    # distance vectors as electronic table (Table A1).
    # In this case we need element-level almost exact match for this subspace.
    assert numpy.allclose(their_v, ref_ev.res.x[3:])
    # Notice that the full-dimension covariance table is available in exact
    # values (Table A2).
    assert numpy.allclose(their_cov, cov)
    # Check that the full optimization result matches.  This is done by
    # concatenating the target vector array with the three values from Table 3
    # of the M16 publication on the head.
    their_full_v = numpy.concatenate(((0.125, -2.58, -0.052), their_v))
    kld = kldgaussian(their_full_v, their_cov, ref_ev.res.x, cov)
    # In this case, because the first 3 target values are not exact, we check
    # the KL divergence in nats.
    assert 0.0 <= kld <= 0.004
