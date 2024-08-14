# These tests are meant to be run inside the docker container as they
# use the cmaq executable and reference paths only available in the docker container.
#
# Run `make tests/test-data/mcip` before running these tests


import subprocess
import dotenv
import os
from pathlib import Path
import pytest


@pytest.fixture
def test_data_dir():
    return Path(__file__).parents[1] / "tests" / "test-data"


def prep_run_dir(run_dir: Path):
    output_dir = run_dir / "output"
    chkpnt_dir = run_dir / "chkpnt"

    output_dir.mkdir()
    chkpnt_dir.mkdir()

    return chkpnt_dir, output_dir


def run_cmaq(cmd, env_file, overrides):
    # Variables to copy across from current environment
    current_env = {
        k: os.environ[k]
        for k in [
            "HOSTNAME",
            "LD_LIBRARY_PATH",
            "PATH",
            "RUN_DIR",
            "DATA_DIR",
        ]
    }
    env_values = {
        **current_env,
        **dotenv.dotenv_values(env_file),
        **overrides,
    }

    res = subprocess.run(
        cmd,
        shell=True,
        executable="/bin/csh",
        capture_output=True,
        text=True,
        env=env_values,
        cwd=current_env["RUN_DIR"],
    )

    print(res.stdout)
    print(res.stderr)

    if os.path.exists(env_values["LOGFILE"]):
        with open(env_values["LOGFILE"], "r") as f:
            print(f.read())

    return res


def test_run_cmaq_fwd(monkeypatch, tmp_path, test_data_dir):
    run_cmd = "/opt/cmaq/cmaq_adj/BLD_fwd_CH4only/ADJOINT_FWD"

    monkeypatch.setenv("RUN_DIR", str(tmp_path))
    monkeypatch.setenv("DATA_DIR", str(test_data_dir))

    prep_run_dir(tmp_path)

    res = run_cmaq(run_cmd, test_data_dir / "env_fwd_2022-07-23.txt", {})

    assert res.returncode == 0


def test_run_cmaq_fwd_mp(monkeypatch, tmp_path, test_data_dir):
    run_cmd = "mpirun -np 2 /opt/cmaq/cmaq_adj/BLD_fwd_CH4only/ADJOINT_FWD"

    monkeypatch.setenv("RUN_DIR", str(tmp_path))
    monkeypatch.setenv("DATA_DIR", str(test_data_dir))

    prep_run_dir(tmp_path)

    res = run_cmaq(
        run_cmd,
        test_data_dir / "env_fwd_2022-07-23.txt",
        {"NPCOL_NPROW": "2 1"},
    )

    assert res.returncode == 0
