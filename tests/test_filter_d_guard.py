"""A2 — Filter D must not swallow a BLAST failure.

`_check_construct_host_coverage` returns `is_fp`.  Returning `False` means "this
site is NOT a construct+host false positive", so a crashed `blastn` that returns
`False` silently lets the site through to CANDIDATE.  Filters A and B already
log their failures (`filter_a_host.py`, `filter_b_flank.py`); D was missed when
that hardening pass landed.
"""
from __future__ import annotations

import subprocess

import pytest

from scripts.s05 import filter_d_altlocus
from scripts.s05.filter_d_altlocus import _check_construct_host_coverage


@pytest.fixture()
def inputs(tmp_path):
    insert = tmp_path / "insert.fa"
    insert.write_text(">contig_1\n" + "ACGT" * 100 + "\n")
    db = tmp_path / "elements.fa"
    db.write_text(">elem\n" + "ACGT" * 50 + "\n")
    return insert, db, tmp_path


def test_blast_failure_is_logged_and_fails_open(inputs, monkeypatch, capsys):
    insert, db, workdir = inputs

    def _boom(cmd, **kwargs):
        return subprocess.CompletedProcess(
            cmd, returncode=2, stdout="", stderr="BLAST Database error: no alias",
        )

    monkeypatch.setattr(filter_d_altlocus.subprocess, "run", _boom)

    is_fp, construct_frac, host_frac, combined_frac = _check_construct_host_coverage(
        insert, db, host_fraction=0.4, host_bp=160, insert_len=400, n_count=0,
        workdir=workdir,
    )

    assert (is_fp, construct_frac) == (False, 0.0)
    # host evidence measured by Filter A must survive the Filter D failure
    assert host_frac == 0.4
    assert combined_frac == 0.4

    err = capsys.readouterr().err
    assert "Filter D" in err
    assert "rc=2" in err
    assert "BLAST Database error" in err


def test_blast_failure_stderr_is_not_discarded(inputs, monkeypatch):
    """The subprocess must be invoked with stderr captured, not DEVNULL."""
    insert, db, workdir = inputs
    seen: dict = {}

    def _spy(cmd, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr(filter_d_altlocus.subprocess, "run", _spy)

    _check_construct_host_coverage(
        insert, db, host_fraction=0.0, host_bp=0, insert_len=400, n_count=0,
        workdir=workdir,
    )

    assert seen.get("stderr") == subprocess.PIPE
    assert seen.get("stderr") is not subprocess.DEVNULL


def test_no_num_threads_flag_in_subject_mode(inputs, monkeypatch):
    """`threads` is accepted for signature parity but must not reach blastn.

    Filter D BLASTs with `-subject`, not `-db`.  blastn 2.16.0 answers a
    `-num_threads` in that mode with

        Warning: [blastn] 'num_threads' is currently ignored when 'subject'
        is specified.

    so passing it buys nothing and pollutes stderr on every site.
    """
    insert, db, workdir = inputs
    seen: dict = {}

    def _spy(cmd, **kwargs):
        seen["cmd"] = cmd
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr(filter_d_altlocus.subprocess, "run", _spy)

    _check_construct_host_coverage(
        insert, db, host_fraction=0.0, host_bp=0, insert_len=400, n_count=0,
        workdir=workdir, threads=11,
    )

    assert "-subject" in seen["cmd"]
    assert "-num_threads" not in seen["cmd"]


def test_successful_blast_still_computes_fractions(inputs, monkeypatch):
    insert, db, workdir = inputs

    def _hit(cmd, **kwargs):
        # 200 bp of the 400 bp insert covered at 95% identity
        return subprocess.CompletedProcess(
            cmd, returncode=0, stdout="1\t200\t95.0\n", stderr="",
        )

    monkeypatch.setattr(filter_d_altlocus.subprocess, "run", _hit)

    is_fp, construct_frac, host_frac, combined_frac = _check_construct_host_coverage(
        insert, db, host_fraction=0.6, host_bp=240, insert_len=400, n_count=0,
        workdir=workdir,
    )

    assert construct_frac == pytest.approx(0.5)
    assert combined_frac == pytest.approx(1.0)
    assert is_fp is True
