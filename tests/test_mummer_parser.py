"""Tests for MummerChecker.parse_alignment_results.

These don't run nucmer — they hand-write a minimal nucmer .delta file
and verify it parses correctly.
"""
import pytest

from multiplex_checker.MultiplexCheckerMummer import MummerChecker


def write_delta(path, lines):
    """Write a .delta file. `lines` is an iterable of literal lines."""
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def test_mummer_parse_forward_hit(tmp_path):
    checker = MummerChecker()
    prefix = tmp_path / "aln"
    write_delta(tmp_path / "aln.delta", [
        "/path/to/ref /path/to/primers",
        "NUCMER",
        ">chr1 P1 1000 20",
        "100 119 1 20 0 0 0",
        "0",
    ])
    results = checker.parse_alignment_results(str(prefix))
    assert len(results) == 1
    r = results[0]
    assert r["primer_name"] == "P1"
    assert r["chromosome"] == "chr1"
    assert r["ref_start"] == 100
    assert r["ref_end"] == 119
    # primer_start=1 < primer_end=20 → forward
    assert r["strand"] == "forward"
    assert r["is_forward_strand"] is True


def test_mummer_parse_reverse_hit(tmp_path):
    checker = MummerChecker()
    prefix = tmp_path / "aln"
    # On the reverse strand nucmer flips the query coordinates: 20 → 1
    write_delta(tmp_path / "aln.delta", [
        "/path/to/ref /path/to/primers",
        "NUCMER",
        ">chr1 P1 1000 20",
        "100 119 20 1 0 0 0",
        "0",
    ])
    results = checker.parse_alignment_results(str(prefix))
    assert len(results) == 1
    assert results[0]["strand"] == "reverse"
    assert results[0]["is_forward_strand"] is False


def test_mummer_parse_multiple_blocks(tmp_path):
    checker = MummerChecker()
    prefix = tmp_path / "aln"
    write_delta(tmp_path / "aln.delta", [
        "/path/to/ref /path/to/primers",
        "NUCMER",
        ">chr1 P1 1000 20",
        "100 119 1 20 0 0 0",
        "0",
        ">chr1 P2 1000 20",
        "200 219 1 20 0 0 0",
        "0",
        ">chr2 P1 1000 20",
        "300 319 1 20 0 0 0",
        "0",
    ])
    results = checker.parse_alignment_results(str(prefix))
    assert len(results) == 3
    # All three hits captured
    starts = sorted(r["ref_start"] for r in results)
    assert starts == [100, 200, 300]


def test_mummer_parse_indel_distance_lines_ignored(tmp_path):
    # Within a coordinate record, nucmer emits per-indel offsets as single
    # integers. The parser must not treat these as coordinate lines.
    checker = MummerChecker()
    prefix = tmp_path / "aln"
    write_delta(tmp_path / "aln.delta", [
        "/path/to/ref /path/to/primers",
        "NUCMER",
        ">chr1 P1 1000 20",
        "100 119 1 20 0 0 0",
        "5",
        "-3",
        "0",
    ])
    results = checker.parse_alignment_results(str(prefix))
    assert len(results) == 1
    assert results[0]["ref_start"] == 100


def test_mummer_parse_missing_delta_raises(tmp_path):
    checker = MummerChecker()
    with pytest.raises(FileNotFoundError):
        checker.parse_alignment_results(str(tmp_path / "nonexistent"))


def test_mummer_parse_empty_alignments(tmp_path):
    # No alignment blocks — just the file header
    checker = MummerChecker()
    prefix = tmp_path / "aln"
    write_delta(tmp_path / "aln.delta", [
        "/path/to/ref /path/to/primers",
        "NUCMER",
    ])
    results = checker.parse_alignment_results(str(prefix))
    assert results == []
