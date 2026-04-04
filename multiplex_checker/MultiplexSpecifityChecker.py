import time
import pandas as pd
import sys
import os
import primer3
import csv
import subprocess
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from abc import ABC, abstractmethod

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class NotImplementedException(Exception):
    pass

@dataclass
class PrimerHit:
    """Data class to represent a primer hit on the reference genome"""
    primer_name: str
    chromosome: str
    start: int
    end: int
    strand: str
    tm: float


@dataclass
class Amplicon:
    """Data class to represent a potential amplicon"""
    chromosome: str
    forward_primer: str
    reverse_primer: str
    start: int
    end: int
    size: int
    forward_tm: float
    reverse_tm: float


class PCRSpecificityChecker(ABC):
    def __init__(self, tm_threshold: float = 30.0, max_amplicon_size: int = 1000,
                 min_amplicon_size: int = 50, region_padding: int = 50):
        """
        Initialize the PCR specificity checker.

        Args:
            tm_threshold: Minimum melting temperature for primer binding
            max_amplicon_size: Maximum allowed amplicon size
            min_amplicon_size: Minimum allowed amplicon size
            region_padding: Nucleotides to extract around binding site for Tm calculation
        """

        self.tm_threshold = tm_threshold
        self.max_amplicon_size = max_amplicon_size
        self.min_amplicon_size = min_amplicon_size
        self.region_padding = region_padding
        self.reference_cache = {}  # Cache for reference sequences
        self.primer_cache = {}  # Cache for primer sequences

    def complement(self, seq: str) -> str:
        """Returns the complement of a DNA sequence."""
        base_complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        try:
            return ''.join(base_complement[base.upper()] for base in seq)
        except KeyError as e:
            raise ValueError(f"Invalid base '{e.args[0]}' in sequence: {seq}") from None

    def reverse_complement(self, seq: str) -> str:
        """Returns the reverse complement of a DNA sequence."""
        return self.complement(seq)[::-1]

    def validate_sequence(self, seq: str, seq_name: str = "sequence") -> bool:
        """Validate that sequence contains only valid DNA bases."""
        valid_bases = set('ATGCN')
        invalid_bases = set(seq.upper()) - valid_bases

        if invalid_bases:
            return False
        return True

    def load_fasta_sequences(self, fasta_path: str) -> Dict[str, str]:
        """
        Load all sequences from a FASTA file into memory for efficient access.

        Args:
            fasta_path: Path to the FASTA file

        Returns:
            Dictionary mapping sequence names to sequences
        """
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        sequences = {}
        current_name = None
        current_seq = []

        logger.info(f"Loading sequences from {fasta_path}")

        try:
            with open(fasta_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        # Save previous sequence
                        if current_name is not None:
                            seq = ''.join(current_seq)
                            if self.validate_sequence(seq, current_name):
                                sequences[current_name] = seq.upper()
                            else:
                                logger.warning(f"Skipping invalid sequence: {current_name}")

                        # Start new sequence
                        current_name = line[1:].split()[0]
                        current_seq = []
                    else:
                        if current_name is None:
                            raise ValueError(f"Sequence data found before header at line {line_num}")
                        current_seq.append(line)

                # Save last sequence
                if current_name is not None:
                    seq = ''.join(current_seq)
                    if self.validate_sequence(seq, current_name):
                        sequences[current_name] = seq.upper()

        except Exception as e:
            raise RuntimeError(f"Error reading FASTA file {fasta_path}: {e}")
        logger.info(f"Loaded {len(sequences)} sequences")
        return sequences

    def get_sequence_by_name(self, fasta_path: str, name: str) -> Optional[str]:
        """
        Get a specific sequence by name, using cache for efficiency.

        Args:
            fasta_path: Path to the FASTA file
            name: Name of the sequence to retrieve

        Returns:
            The sequence string or None if not found
        """
        # Use cache if available
        if fasta_path not in self.reference_cache:
            self.reference_cache[fasta_path] = self.load_fasta_sequences(fasta_path)

        return self.reference_cache[fasta_path].get(name)

    def get_region_sequence(self, seq: str, start: int, end: int) -> str:
        """
        Extract a region from a sequence with padding for Tm calculation.

        Args:
            seq: The full sequence
            start: Start position (1-based)
            end: End position (1-based)

        Returns:
            The extracted region with padding
        """
        # Convert to 0-based indexing and add padding
        start_idx = max(0, start - 1 - self.region_padding)
        end_idx = min(len(seq), end + self.region_padding)
        return seq[start_idx:end_idx]

    @abstractmethod
    def run_alignment(self, ref_fasta: str, primer_fasta: str, prefix: str = 'primer_alignment') -> bool:
        """
        Run alignment between reference and primer sequences.

        Args:
            ref_fasta: Path to reference FASTA file
            primer_fasta: Path to primer FASTA file
            prefix: Output file prefix

        Returns:
            True if alignment was successful, False otherwise
        """
        raise NotImplementedError('Implement run_alignment method')


    @abstractmethod
    def parse_alignment_results(self, prefix: str) -> List[Dict]:
        """
        Parse alignment file and extract alignment coordinates.

        Args:
            prefix: Path to alignment file

        Returns:
            List of dictionaries containing raw alignment data


        alignment = {
        'primer_name': current_primer,
        'chromosome': current_chr,
        'ref_start': ref_start,
        'ref_end': ref_end,
        'primer_start': primer_start,
        'primer_end': primer_end,
        'strand': strand,
        'is_forward_strand': is_forward_strand
    }
        """
        raise NotImplementedError('Implement parse_alignment_results method')

    def calculate_primer_thermodynamics(self, alignments: List[Dict], ref_fasta: str, primer_fasta: str) -> List[
        PrimerHit]:
        """
        Calculate binding thermodynamics for parsed alignments and filter by Tm threshold.

        Args:
            alignments: List of alignment dictionaries from parse_alignment_results
            ref_fasta: Path to reference FASTA file
            primer_fasta: Path to primer FASTA file

        Returns:
            List of PrimerHit objects representing valid primer binding sites
        """

        # Load primer sequences
        if primer_fasta not in self.primer_cache:
            self.primer_cache[primer_fasta] = self.load_fasta_sequences(primer_fasta)
        primer_sequences = self.primer_cache[primer_fasta]

        hits = []

        logger.info(f"Calculating thermodynamics for {len(alignments)} alignments")

        for alignment in alignments:
            primer_name = alignment['primer_name']
            chromosome = alignment['chromosome']
            ref_start = alignment['ref_start']
            ref_end = alignment['ref_end']
            strand = alignment['strand']
            is_forward_strand = alignment['is_forward_strand']
            # Get reference sequence region
            chr_seq = self.get_sequence_by_name(ref_fasta, chromosome)
            if chr_seq is None:

                logger.warning(f"Chromosome {chromosome} not found in reference")
                continue

            # Get primer sequence
            primer_seq = primer_sequences.get(primer_name)
            if primer_seq is None:
                logger.warning(f"Primer {primer_name} not found in primer file")
                continue

            # Extract binding region
            region_seq = self.get_region_sequence(chr_seq, ref_start, ref_end)

            # Skip regions with too many N's
            if region_seq.count('N') > 1:
                continue

            # Prepare template sequence for Tm calculation
            if is_forward_strand:
                template_for_tm = self.reverse_complement(region_seq)
            else:
                template_for_tm = region_seq

            # Calculate melting temperature
            try:
                binding_tm = primer3.bindings.calc_end_stability(primer_seq, template_for_tm).tm

                # Only keep hits above threshold
                if binding_tm >= self.tm_threshold:
                    hit = PrimerHit(
                        primer_name=primer_name,
                        chromosome=chromosome,
                        start=ref_start,
                        end=ref_end,
                        strand=strand,
                        tm=binding_tm
                    )
                    hits.append(hit)

            except Exception as e:
                logger.debug(f"Failed to calculate Tm for {primer_name} at {chromosome}:{ref_start}-{ref_end}: {e}")
                continue

        logger.info(f"Found {len(hits)} valid primer binding sites")
        return hits

    def parse_results(self, prefix: str, ref_fasta: str, primer_fasta: str) -> List[PrimerHit]:
        """
        Parse alignment file and calculate binding thermodynamics.

        This is now a high-level orchestrator function that combines parsing and thermodynamics.
        Args:
            prefix: Path to alignment file
            ref_fasta: Path to reference FASTA file
            primer_fasta: Path to primer FASTA file

        Returns:
            List of PrimerHit objects representing valid primer binding sites
        """

        # Parse the delta file to get raw alignments
        st = time.time()
        alignments = self.parse_alignment_results(prefix)
        end = time.time()
        print('parsing:',end-st)

        # Calculate thermodynamics and filter by Tm threshold
        st = time.time()
        hits = self.calculate_primer_thermodynamics(alignments, ref_fasta, primer_fasta)
        end = time.time()
        print('thermodynamic calculations:', end - st)
        return hits

    def find_amplicons(self, hits: List[PrimerHit]) -> List[Amplicon]:
        """
           O(F log F + R log R) for sorting
           O(F + R) for pairing.
           :param hits:
           :return:
           """
        amplicons = []

        # Group hits by chromosome
        chr_hits = {}
        for hit in hits:
            if hit.chromosome not in chr_hits:
                chr_hits[hit.chromosome] = []
            chr_hits[hit.chromosome].append(hit)

        logger.info("Searching for unwanted amplicons...")

        # Look for forward-reverse pairs on same chromosome
        for chr_name, sites in chr_hits.items():
            forward_hits = [h for h in sites if h.strand == 'forward']
            reverse_hits = [h for h in sites if h.strand == 'reverse']

            forward_hits.sort(key=lambda h: h.start)
            reverse_hits.sort(key=lambda h: h.end)

            logger.info(
                f"Sequence: {chr_name},found {len(forward_hits)} forward bindings, {len(reverse_hits)} reverse bindings")

            start_idx = 0  # reverse list window start
            end_idx = 0  # reverse list window end
            len_r = len(reverse_hits)

            for fwd_hit in forward_hits:
                # get the window
                min_pos = fwd_hit.start + self.min_amplicon_size
                max_pos = fwd_hit.start + self.max_amplicon_size
                while start_idx < len_r and reverse_hits[start_idx].end < min_pos:
                    start_idx += 1

                while end_idx < len_r and reverse_hits[end_idx].end <= max_pos:
                    end_idx += 1

                for rev_hit in reverse_hits[start_idx:end_idx]:
                    # Double check size constraints (not necessary)
                    if not fwd_hit.start < rev_hit.end:
                        continue

                    amplicon_size = rev_hit.end - fwd_hit.start + 1

                    # Double check size constraints (not necessary)
                    if self.min_amplicon_size <= amplicon_size <= self.max_amplicon_size:
                        amplicon = Amplicon(
                            chromosome=chr_name,
                            forward_primer=fwd_hit.primer_name,
                            reverse_primer=rev_hit.primer_name,
                            start=fwd_hit.start,
                            end=rev_hit.end,
                            size=amplicon_size,
                            forward_tm=fwd_hit.tm,
                            reverse_tm=rev_hit.tm
                        )
                        amplicons.append(amplicon)
            logger.info(f"Done: Found {len(amplicons)} amplicons in {chr_name}")

            # Sort by size
            # amplicons.sort(key=lambda x: x.size)

        logger.info(f"Found {len(amplicons)} potential unwanted amplicons")
        return amplicons


    def write_amplicons_to_csv(self, amplicons: List[Amplicon], filename: str = 'unwanted_amplicon_results_20_06_25.csv'):
        """Write amplicons to CSV file with detailed information."""
        try:
            with open(filename, mode='w', newline='') as file:
                fieldnames = [
                    'chromosome', 'forward_primer', 'reverse_primer',
                    'start', 'end', 'size', 'forward_tm', 'reverse_tm',
                    'avg_tm', 'priority'
                ]
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()

                for amplicon in amplicons:
                    avg_tm = (amplicon.forward_tm + amplicon.reverse_tm) / 2

                    writer.writerow({
                        'chromosome': amplicon.chromosome,
                        'forward_primer': amplicon.forward_primer,
                        'reverse_primer': amplicon.reverse_primer,
                        'start': amplicon.start,
                        'end': amplicon.end,
                        'size': amplicon.size,
                        'forward_tm': round(amplicon.forward_tm, 2),
                        'reverse_tm': round(amplicon.reverse_tm, 2),
                        'avg_tm': round(avg_tm, 2),
                    })

            logger.info(f"Amplicons written to {filename}")

        except Exception as e:
            logger.error(f"Error writing amplicons to CSV: {e}")
            raise

    def analyze_specificity(self, ref_fasta: str, primer_fasta: str,
                            prefix: str = 'primer_alignment', run_alignment: bool = True) -> bool:
        """
        Complete analysis pipeline for primer specificity.

        Args:
            ref_fasta: Path to reference genome FASTA
            primer_fasta: Path to primer FASTA
            prefix: Output file prefix
            run_alignment: Whether to run alignment (False to use existing alignment file)

        Returns:
            True if analysis completed successfully
        """
        try:
            logger.info('Started alignment')
            # Step 1: Run alignment if requested
            if run_alignment:
                if not self.run_alignment(ref_fasta, primer_fasta, prefix):
                    logger.error("Alignment failed")
                    return False

            # Step 2: Parse results and find binding sites
            logger.info('Finished alignment, starting parsing')

            hits = self.parse_results(prefix, ref_fasta, primer_fasta)
            logger.info(f"Found {len(hits)} potential amplicons")
            if not hits:
                logger.warning("No primer binding sites found above threshold")
                return True

            # Step 3: Find unwanted amplicons
            amplicons = self.find_amplicons(hits)

            # Step 4: Write results
            output_file = f"{prefix}_amplicon_results.csv"
            self.write_amplicons_to_csv(amplicons, output_file)

            # Step 5: Summary
            logger.info(f"Analysis complete:")
            logger.info(f"  - Total primer binding sites: {len(hits)}")
            logger.info(f"  - Potential unwanted amplicons: {len(amplicons)}")
            logger.info(f"  - Results written to: {output_file}")

            if amplicons:
                # Show top 5 most problematic amplicons
                logger.info("Top 5 most problematic amplicons:")
                for i, amp in enumerate(amplicons[:5], 1):
                    avg_tm = (amp.forward_tm + amp.reverse_tm) / 2
                    logger.info(f"  {i}. {amp.forward_primer} + {amp.reverse_primer}: "
                                f"{amp.size}bp, Tm={avg_tm:.1f}°C")

            return True

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            return False