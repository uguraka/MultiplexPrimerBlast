"""Tests for BlastChecker.parse_alignment_results.

These don't run blastn — they hand-write a BLAST -outfmt 7 tabular file
and verify it parses correctly.
"""
import pytest

from multiplex_checker.MultiplexCheckerBlast import BlastChecker


BLAST_COLUMNS = (
    "qseqid", "sseqid", "qstart", "qend",
    "sstart", "send", "length", "pident",
    "evalue", "sstrand",
)


def write_blast_tab(path, rows, include_header_comment=True):
    """Write a -outfmt 7 file. Each row is a tuple matching BLAST_COLUMNS."""
    with open(path, "w") as f:
        if include_header_comment:
            f.write("# BLASTN output\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def test_blast_parse_forward_hit(tmp_path):
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    write_blast_tab(blast_file, [
        ("P1", "chr1", 1, 20, 100, 119, 20, 100.0, 1e-5, "plus"),
    ])
    results = checker.parse_alignment_results(str(blast_file))
    assert len(results) == 1
    r = results[0]
    assert r["primer_name"] == "P1"
    assert r["chromosome"] == "chr1"
    assert r["ref_start"] == 100
    assert r["ref_end"] == 119
    assert r["strand"] == "forward"
    assert r["is_forward_strand"] is True


def test_blast_parse_reverse_hit(tmp_path):
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    # On the minus strand, BLAST emits sstart > send. The parser should
    # normalise these into (min, max).
    write_blast_tab(blast_file, [
        ("P1", "chr1", 1, 20, 200, 181, 20, 100.0, 1e-5, "minus"),
    ])
    results = checker.parse_alignment_results(str(blast_file))
    assert len(results) == 1
    r = results[0]
    assert r["strand"] == "reverse"
    assert r["is_forward_strand"] is False
    assert r["ref_start"] == 181  # normalised min
    assert r["ref_end"] == 200    # normalised max


def test_blast_parse_filters_low_identity(tmp_path):
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    write_blast_tab(blast_file, [
        ("P1", "chr1", 1, 20, 100, 119, 20, 100.0, 1e-5, "plus"),
        ("P1", "chr1", 1, 20, 200, 219, 20,  75.0, 1e-3, "plus"),  # < 80%
    ])
    results = checker.parse_alignment_results(str(blast_file))
    assert len(results) == 1
    assert results[0]["ref_start"] == 100


def test_blast_parse_identity_boundary_inclusive(tmp_path):
    # Exactly 80% should be kept (>= 80)
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    write_blast_tab(blast_file, [
        ("P1", "chr1", 1, 20, 100, 119, 20, 80.0, 1e-5, "plus"),
    ])
    results = checker.parse_alignment_results(str(blast_file))
    assert len(results) == 1


def test_blast_parse_empty_file(tmp_path):
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    blast_file.write_text("# header only, no hits\n")
    results = checker.parse_alignment_results(str(blast_file))
    assert results == []


def test_blast_parse_missing_file(tmp_path):
    # The parser logs an error and returns [] rather than raising
    checker = BlastChecker()
    results = checker.parse_alignment_results(str(tmp_path / "missing.tab"))
    assert results == []


def test_blast_parse_numeric_chromosome_stays_string(tmp_path):
    # Numeric-looking sseqid values must not be coerced to int
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    write_blast_tab(blast_file, [
        ("P1", "12345", 1, 20, 100, 119, 20, 100.0, 1e-5, "plus"),
    ])
    results = checker.parse_alignment_results(str(blast_file))
    assert results[0]["chromosome"] == "12345"
    assert isinstance(results[0]["chromosome"], str)


def test_blast_parse_multiple_hits_preserved_order(tmp_path):
    checker = BlastChecker()
    blast_file = tmp_path / "out.tab"
    write_blast_tab(blast_file, [
        ("P1", "chr1", 1, 20, 100, 119, 20, 100.0, 1e-5, "plus"),
        ("P2", "chr1", 1, 20, 500, 519, 20,  95.0, 1e-4, "plus"),
        ("P3", "chr2", 1, 20, 300, 281, 20,  90.0, 1e-3, "minus"),
    ])
    results = checker.parse_alignment_results(str(blast_file))
    assert len(results) == 3
    assert [r["primer_name"] for r in results] == ["P1", "P2", "P3"]
