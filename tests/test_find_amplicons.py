"""Tests for the core amplicon-discovery algorithm.

`find_amplicons` is MPB's unique algorithmic contribution. These tests
cover the algorithm itself — they do NOT exercise alignment, Tm
calculation, or any file I/O.
"""
from tests.conftest import make_hit


def test_empty_input(checker):
    assert checker.find_amplicons([]) == []


def test_only_forward_hits(checker):
    # No reverse hits → nothing to pair with
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("F2", "chr1", 200, 220, "forward"),
    ]
    assert checker.find_amplicons(hits) == []


def test_only_reverse_hits(checker):
    hits = [
        make_hit("R1", "chr1", 100, 120, "reverse"),
        make_hit("R2", "chr1", 200, 220, "reverse"),
    ]
    assert checker.find_amplicons(hits) == []


def test_different_chromosomes_dont_pair(checker):
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("R1", "chr2", 280, 300, "reverse"),
    ]
    assert checker.find_amplicons(hits) == []


def test_happy_path_single_amplicon(checker):
    # forward starts at 100, reverse ends at 300 → size = 300 - 100 + 1 = 201
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("R1", "chr1", 280, 300, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 1
    a = amplicons[0]
    assert a.chromosome == "chr1"
    assert a.forward_primer == "F1"
    assert a.reverse_primer == "R1"
    assert a.start == 100
    assert a.end == 300
    assert a.size == 201


def test_reverse_upstream_of_forward_not_paired(checker):
    # Reverse hit ends before forward starts — would be negative size, must be excluded
    hits = [
        make_hit("F1", "chr1", 500, 520, "forward"),
        make_hit("R1", "chr1", 100, 120, "reverse"),
    ]
    assert checker.find_amplicons(hits) == []


def test_size_at_min_boundary_is_included(checker):
    # min_amplicon_size=50. Forward at 100, reverse end at 149 → size = 50
    hits = [
        make_hit("F1", "chr1", 100, 110, "forward"),
        make_hit("R1", "chr1", 140, 149, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 1
    assert amplicons[0].size == 50


def test_size_below_min_is_excluded(checker):
    # size = 49, below min of 50
    hits = [
        make_hit("F1", "chr1", 100, 110, "forward"),
        make_hit("R1", "chr1", 140, 148, "reverse"),
    ]
    assert checker.find_amplicons(hits) == []


def test_size_at_max_boundary_is_included(checker):
    # max_amplicon_size=1000. Forward at 100, reverse end at 1099 → size = 1000
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("R1", "chr1", 1080, 1099, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 1
    assert amplicons[0].size == 1000


def test_size_above_max_is_excluded(checker):
    # size = 1001
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("R1", "chr1", 1080, 1100, "reverse"),
    ]
    assert checker.find_amplicons(hits) == []


def test_one_forward_multiple_reverses(checker):
    # F1 should pair with both R1 and R2
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("R1", "chr1", 280, 300, "reverse"),
        make_hit("R2", "chr1", 480, 500, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 2
    sizes = sorted(a.size for a in amplicons)
    assert sizes == [201, 401]


def test_multiple_forwards_one_reverse(checker):
    # Both F1 and F2 should pair with R1
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("F2", "chr1", 200, 220, "forward"),
        make_hit("R1", "chr1", 380, 400, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 2
    sizes = sorted(a.size for a in amplicons)
    # F1+R1: 400-100+1 = 301
    # F2+R1: 400-200+1 = 201
    assert sizes == [201, 301]


def test_multiple_chromosomes_independent(checker):
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("R1", "chr1", 280, 300, "reverse"),
        make_hit("F2", "chr2", 500, 520, "forward"),
        make_hit("R2", "chr2", 680, 700, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 2
    chroms = sorted(a.chromosome for a in amplicons)
    assert chroms == ["chr1", "chr2"]


def test_sliding_window_handles_interleaved_hits(checker):
    # Three forwards and three reverses on one chromosome.
    # F1@100 — pairs with R1(end=300, size=201), R2(end=500, size=401), R3(end=900, size=801)
    # F2@400 — pairs with R2(end=500, size=101), R3(end=900, size=501)
    # F3@600 — pairs with R3(end=900, size=301)
    # Total: 6
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward"),
        make_hit("F2", "chr1", 400, 420, "forward"),
        make_hit("F3", "chr1", 600, 620, "forward"),
        make_hit("R1", "chr1", 280, 300, "reverse"),
        make_hit("R2", "chr1", 480, 500, "reverse"),
        make_hit("R3", "chr1", 880, 900, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 6


def test_self_pair_is_reported(checker):
    # The same primer name can appear as both forward and reverse hit
    # (it bound the + strand somewhere and the - strand elsewhere).
    # MPB makes no special case for this — pair is reported.
    hits = [
        make_hit("P1", "chr1", 100, 120, "forward"),
        make_hit("P1", "chr1", 280, 300, "reverse"),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 1
    a = amplicons[0]
    assert a.forward_primer == "P1"
    assert a.reverse_primer == "P1"


def test_tm_values_propagated_to_amplicon(checker):
    hits = [
        make_hit("F1", "chr1", 100, 120, "forward", tm=58.4),
        make_hit("R1", "chr1", 280, 300, "reverse", tm=61.2),
    ]
    amplicons = checker.find_amplicons(hits)
    assert len(amplicons) == 1
    assert amplicons[0].forward_tm == 58.4
    assert amplicons[0].reverse_tm == 61.2
