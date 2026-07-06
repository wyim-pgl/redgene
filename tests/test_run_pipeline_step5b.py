"""Dispatch tests for step 5b wiring in run_pipeline.py."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import run_pipeline  # noqa: E402


def test_step_order_has_5b_after_5():
    order = run_pipeline.STEP_ORDER
    assert "5b" in order
    assert order.index("5b") == order.index("5") + 1


def test_step_scripts_has_5b():
    assert run_pipeline.STEP_SCRIPTS["5b"] == "scripts/s05b_construct_assembly.py"


def test_parse_steps_excludes_5b_from_1_5_but_includes_in_1_5b():
    assert "5b" not in run_pipeline.parse_steps("1-5")
    assert "5b" in run_pipeline.parse_steps("1-5b")
    assert run_pipeline.parse_steps("5b") == ["5b"]


def _sample_cfg():
    return {
        "host_reference": "db/host.fa",
        "construct_reference": "db/construct.fa",
        "reads": {"r1": "x_R1.fq.gz", "r2": "x_R2.fq.gz"},
    }


def test_build_step_cmd_5b_argv(tmp_path):
    cmd = run_pipeline.build_step_cmd(
        "5b", "tsample", _sample_cfg(),
        outdir=tmp_path / "results", threads=4, base_dir=Path("/repo"),
        no_remote_blast=True,
    )
    assert "scripts/s05b_construct_assembly.py" in " ".join(cmd)
    assert "--contigs-all" in cmd and "--s05-dir" in cmd
    assert "--sample-name" in cmd and "tsample" in cmd
    assert "--remote-blast" not in cmd


def test_build_step_cmd_5b_remote_when_enabled(tmp_path):
    cmd = run_pipeline.build_step_cmd(
        "5b", "tsample", _sample_cfg(),
        outdir=tmp_path / "results", threads=4, base_dir=Path("/repo"),
        no_remote_blast=False,
    )
    assert "--remote-blast" in cmd


def test_spades_memory_gb_from_slurm(monkeypatch):
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "98304")  # 96 GiB in MB
    # 80% of 96 GiB = 76 GiB
    assert run_pipeline._spades_memory_gb(16) == 76


def test_spades_memory_gb_fallback_heuristic(monkeypatch):
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    assert run_pipeline._spades_memory_gb(16) == 32  # 16*4//2


def test_spades_memory_gb_bad_env_falls_back(monkeypatch):
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "not-a-number")
    assert run_pipeline._spades_memory_gb(8) == 16  # 8*4//2
