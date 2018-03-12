"""Testing invocation of the script ``jlacompress``."""
import subprocess
import pytest
from lsnz_test_infra import jla_full_paths, outdir, ref_ev, ref_ev_nologdet


@pytest.fixture(scope="module")
def basic_cmd(jla_full_paths, outdir):
    dpath, fpath = jla_full_paths
    return ["jlacompress", "-d", dpath, "-t", fpath, "-v", "-p",
            "%s/test_invoke_script-" % outdir]


def test_run_script_help():
    cmd = ["jlacompress", "-h"]
    status = subprocess.call(cmd)
    assert status == 0


@pytest.mark.parametrize("additional_args", [[], ["-n"]])
def test_run_script(basic_cmd, additional_args):
    status = subprocess.call(basic_cmd + additional_args)
    assert status == 0
