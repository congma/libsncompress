"""Testing invocation of the script ``jlacompress``."""
import os
import subprocess
import six
import pytest
from lsnz_test_infra import jla_full_paths, outdir, ref_ev, ref_ev_nologdet


def config_to_cmdlist(config_dict):
    cmdlist = []
    for option, value in six.iteritems(config_dict):
        cmdlist.append(option)
        if value is not None:
            cmdlist.append(value)
    return cmdlist


@pytest.fixture(scope="module")
def cmd_config(jla_full_paths, outdir):
    dpath, fpath = jla_full_paths
    config = {"-d": dpath,
              "-t": fpath,
              "-v": None,
              "-p": "%s/test_invoke_script-" % outdir}
    return config


def test_run_script_help():
    cmd = ["jlacompress", "-h"]
    status = subprocess.call(cmd)
    assert status == 0


@pytest.mark.parametrize("additional_args", [[], ["-n"]])
def test_run_script(cmd_config, additional_args):
    cmd = ["jlacompress"] + config_to_cmdlist(cmd_config) + additional_args
    status = subprocess.call(cmd)
    assert status == 0


def test_new_dir_in_prefix(cmd_config, outdir):
    d = dict(cmd_config)
    test_dir = "%s/nested/new_dir/" % outdir
    d["-p"] = test_dir
    cmd = ["jlacompress"] + config_to_cmdlist(d)
    status = subprocess.call(cmd)
    assert status == 0
    dircontent = frozenset(os.listdir(test_dir))
    for name in ("redshift", "mean", "cov"):
        assert ("%s.txt" % name) in dircontent
