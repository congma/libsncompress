"""Testing usage of libsncompress.base"""
import os
import os.path
import pytest
import numpy
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
    base = libsncompress.BinnedSN(*jla_full_paths)


def test_binning_that_selects_nothing(jla_full_paths):
    with pytest.raises(ValueError):
        base = libsncompress.BinnedSN(*jla_full_paths,
                                      logbins=[numpy.log10([1.5, 2.0])])
