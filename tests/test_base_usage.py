"""Testing usage of libsncompress.base"""
import os
import os.path
import pytest
import six
import six.moves as sm
import numpy
from numpy.random import shuffle
import libsncompress
from lsnz_test_infra import jla_full_paths, outdir


@pytest.fixture
def extra_file(jla_full_paths):
    fits_dir = jla_full_paths[0]
    fpath = os.path.join(fits_dir, "test_tmp.dat")
    fh = open(fpath, "wb")
    fh.close()
    yield fpath
    try:
        os.remove(fpath)
    except OSError:
        pass


def test_invalid_fits_dir(outdir, jla_full_paths):
    emptydir = "%s" % outdir.mkdir("invalid")
    with pytest.raises(ValueError):
        base = libsncompress.BinnedSN(emptydir, jla_full_paths[1])


def test_extra_file_in_fits_dir(extra_file, jla_full_paths):
    os.stat(extra_file)
    base = libsncompress.BinnedSN(*jla_full_paths)


def test_binning_that_selects_nothing(jla_full_paths):
    with pytest.raises(ValueError):
        base = libsncompress.BinnedSN(*jla_full_paths,
                                      logbins=[numpy.log10([1.5, 2.0])])


@pytest.mark.parametrize("size", list(sm.range(2, 36)))
def test_bins_sizes(jla_full_paths, size):
    rndlogz = numpy.linspace(-2.0, numpy.log10(1.3), num=size)
    shuffle(rndlogz)
    base = libsncompress.BinnedSN(*jla_full_paths, logbins=[rndlogz])
    # Sanity check for the number of bins.
    assert len(base.binidcontent) == size - 1
    assert base.bins.nbins == size - 1
    # Each data point is mapped (i.e. binned).
    revlookup_range_size = sum(len(v) for v in
                               six.itervalues(base.binidcontent))
    assert revlookup_range_size == base.datadimension
    for lgz in base.logredshifts:
        assert base.bins.searchenum(lgz)[0] is not None
