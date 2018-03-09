"""Testing invocation of the script ``jlacompress``."""
import subprocess
from lsnz_test_infra import jla_full_paths


def test_run_script_help():
    cmd = ["jlacompress", "-h"]
    status = subprocess.call(cmd)
    assert status == 0


def test_run_script(jla_full_paths):
    dpath, fpath = jla_full_paths
    cmd = ["jlacompress", "-d", dpath, "-t", fpath,
           "-v", "-p", "test_run_script"]
    status = subprocess.call(cmd)
    assert status == 0
