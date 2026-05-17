import pytest


def test_load_fasta_multi_record(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1\nATGCATGC\n>chr2\nGGGGCCCC\n")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert seqs == {"chr1": "ATGCATGC", "chr2": "GGGGCCCC"}


def test_load_fasta_multiline_record(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1\nATGC\nATGC\nATGC\n")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert seqs == {"chr1": "ATGCATGCATGC"}


def test_load_fasta_header_uses_first_token(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1 some extra description here\nATGC\n")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert "chr1" in seqs
    assert seqs["chr1"] == "ATGC"


def test_load_fasta_skips_invalid_sequence(checker, tmp_path):
    # chr2 has an invalid Z and should be skipped, others kept
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1\nATGC\n>chr2\nATGCZ\n>chr3\nGGGG\n")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert "chr1" in seqs
    assert "chr2" not in seqs
    assert "chr3" in seqs


def test_load_fasta_uppercases_input(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1\natgc\n")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert seqs == {"chr1": "ATGC"}


def test_load_fasta_empty_file(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text("")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert seqs == {}


def test_load_fasta_blank_lines_ignored(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1\n\nATGC\n\n>chr2\n\nGGGG\n")
    seqs = checker.load_fasta_sequences(str(fasta))
    assert seqs == {"chr1": "ATGC", "chr2": "GGGG"}


def test_load_fasta_missing_file_raises(checker, tmp_path):
    with pytest.raises(FileNotFoundError):
        checker.load_fasta_sequences(str(tmp_path / "does_not_exist.fa"))


def test_load_fasta_sequence_before_header_raises(checker, tmp_path):
    fasta = tmp_path / "test.fa"
    fasta.write_text("ATGC\n>chr1\nGGGG\n")
    with pytest.raises(RuntimeError, match="Sequence data found before header"):
        checker.load_fasta_sequences(str(fasta))
