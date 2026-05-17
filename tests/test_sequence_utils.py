import pytest


def test_complement_canonical(checker):
    assert checker.complement("ATGC") == "TACG"


def test_complement_with_n(checker):
    assert checker.complement("ATGCN") == "TACGN"


def test_complement_lowercase_input(checker):
    assert checker.complement("atgc") == "TACG"


def test_complement_invalid_base_raises(checker):
    with pytest.raises(ValueError, match="Invalid base"):
        checker.complement("ATGX")


def test_reverse_complement(checker):
    assert checker.reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_palindrome(checker):
    # EcoRI site is its own reverse complement
    assert checker.reverse_complement("GAATTC") == "GAATTC"


def test_reverse_complement_with_n(checker):
    assert checker.reverse_complement("NATGC") == "GCATN"


def test_validate_sequence_valid(checker):
    assert checker.validate_sequence("ATGCN") is True


def test_validate_sequence_invalid(checker):
    assert checker.validate_sequence("ATGCZ") is False


def test_validate_sequence_lowercase(checker):
    assert checker.validate_sequence("atgcn") is True


def test_validate_sequence_empty(checker):
    assert checker.validate_sequence("") is True


def test_get_region_normal(checker):
    # 4-base target at positions 101-104 (1-based), 100 bp flank on each side
    seq = "A" * 100 + "GGGG" + "A" * 100
    region = checker.get_region_sequence(seq, 101, 104)
    # 50 bp padding either side = 4 + 50 + 50 = 104
    assert len(region) == 104
    assert "GGGG" in region


def test_get_region_near_start(checker):
    # Target abuts the very start of the sequence — left padding clipped
    seq = "GGGG" + "A" * 100
    region = checker.get_region_sequence(seq, 1, 4)
    assert region.startswith("GGGG")
    assert len(region) == 4 + 50  # no left flank


def test_get_region_near_end(checker):
    # Target abuts the end — right padding clipped
    seq = "A" * 100 + "GGGG"
    region = checker.get_region_sequence(seq, 101, 104)
    assert region.endswith("GGGG")
    assert len(region) == 4 + 50  # no right flank


def test_get_region_padding_zero():
    from multiplex_checker.MultiplexCheckerBlast import BlastChecker
    c = BlastChecker(region_padding=0)
    seq = "ATGCATGCATGC"
    region = c.get_region_sequence(seq, 5, 8)
    # 1-based [5,8] inclusive → 0-based [4,8) → "ATGC"
    assert region == seq[4:8]
