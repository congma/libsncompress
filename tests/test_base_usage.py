"""Testing usage of libsncompress.base"""
import pytest
import libsncompress
from lsnz_test_infra import jla_full_paths, outdir


def test_invalid_fits_dir(outdir, jla_full_paths):
    emptydir = "%s" % outdir.mkdir("invalid")
    with pytest.raises(ValueError):
        base = libsncompress.BinnedSN(emptydir, jla_full_paths[1])
