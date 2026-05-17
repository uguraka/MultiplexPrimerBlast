import pytest

from multiplex_checker.MultiplexCheckerBlast import BlastChecker
from multiplex_checker.MultiplexSpecifityChecker import PrimerHit


@pytest.fixture
def checker():
    """Default BlastChecker — used for testing base-class methods that don't
    depend on any backend-specific behaviour."""
    return BlastChecker(tm_threshold=30.0, min_amplicon_size=50, max_amplicon_size=1000)


def make_hit(primer_name, chromosome, start, end, strand, tm=60.0):
    """Convenience for building PrimerHit objects in tests."""
    return PrimerHit(
        primer_name=primer_name,
        chromosome=chromosome,
        start=start,
        end=end,
        strand=strand,
        tm=tm,
    )
