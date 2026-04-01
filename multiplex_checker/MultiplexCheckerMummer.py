from typing import List, Dict, Optional
import subprocess
import os
from multiplex_checker.MultiplexSpecifityChecker import PCRSpecificityChecker
import logging
import sys


# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MummerChecker(PCRSpecificityChecker):
    def __init__(self, tm_threshold=30.0, max_amplicon_size=650, min_amplicon_size=50, region_padding=50):
        super().__init__(tm_threshold=tm_threshold, max_amplicon_size=max_amplicon_size, min_amplicon_size=min_amplicon_size, region_padding=region_padding)
        self._l = 14 # minimum match length
        self._c = 6  # minimum_match length
        self.num_threads = 25

    def parse_alignment_results(self, prefix: str) -> List[Dict]:
        """
             Parse MUMmer delta file and extract alignment coordinates.

             Args:
                 prefix: Path to MUMmer delta file

             Returns:
                 List of dictionaries containing raw alignment data
             """
        mummer_delta_file = f"{prefix}.delta"
        if not os.path.exists(mummer_delta_file):
            raise FileNotFoundError(f"Delta file not found: {mummer_delta_file}")

        alignments = []
        current_chr = None
        current_primer = None

        logger.info(f"Parsing MUMmer delta file: {mummer_delta_file}")

        try:
            with open(mummer_delta_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # Skip empty lines and header
                    if not line or line.startswith('/') or line.startswith('NUCMER'):
                        continue

                    # End of alignment block
                    if line == '0':
                        continue

                    # New alignment block header
                    if line.startswith('>'):
                        try:
                            parts = line[1:].split()
                            if len(parts) >= 2:
                                current_chr = parts[0]
                                current_primer = parts[1]
                            else:
                                logger.warning(f"Malformed header at line {line_num}: {line}")
                                continue
                        except Exception as e:
                            logger.warning(f"Error parsing header at line {line_num}: {e}")
                            continue

                    # Alignment coordinates
                    elif current_chr and current_primer:
                        try:
                            coords = line.split()
                            if len(coords) >= 4:
                                ref_start = int(coords[0])
                                ref_end = int(coords[1])
                                primer_start = int(coords[2])
                                primer_end = int(coords[3])

                                # Determine strand based on primer coordinates
                                is_forward_strand = primer_start < primer_end
                                strand = 'forward' if is_forward_strand else 'reverse'

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
                                alignments.append(alignment)

                        except (ValueError, IndexError) as e:
                            logger.warning(f"Error parsing coordinates at line {line_num}: {e}")
                            continue

        except Exception as e:
            raise RuntimeError(f"Error parsing MUMmer delta file: {e}")

        logger.info(f"Parsed {len(alignments)} alignments from delta file")
        return alignments

    def run_alignment(self, ref_fasta: str, primer_fasta: str, prefix: str = 'primer_alignment',db_prefix: Optional[str] = None) -> bool:

        """
        Run MUMmer alignment between reference and primer sequences.

        Args:
            ref_fasta: Path to reference FASTA file
            primer_fasta: Path to primer FASTA file
            prefix: Output file prefix

        Returns:
            True if alignment was successful, False otherwise
        """
        if not os.path.exists(ref_fasta):
            raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")
        if not os.path.exists(primer_fasta):
            raise FileNotFoundError(f"Primer FASTA not found: {primer_fasta}")

        logger.info(f"Running MUMmer alignment: {ref_fasta} vs {primer_fasta}")

        try:
            # Run nucmer with appropriate parameters for short primers
            cmd = [
                'nucmer',
                '--maxmatch',  # Find all maximal matches
                '-l', str(self._l),  # Minimum match length
                '-c', str(self._c),  # Minimum cluster length
                '--save', 'abx',  # Save additional info
                ref_fasta,
                primer_fasta,
                '-p', prefix, # Output prefix
                '-t', str(self.num_threads)
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Check if delta file was created and is not empty
            delta_file = f"{prefix}.delta"
            if not os.path.exists(delta_file) or os.path.getsize(delta_file) == 0:
                logger.error("MUMmer alignment failed - empty or missing delta file")
                return False

            logger.info(f"MUMmer alignment completed successfully. Output: {delta_file}")
            return True

        except subprocess.CalledProcessError as e:
            logger.error(f"MUMmer alignment failed: {e}")
            logger.error(f"STDERR: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Error running MUMmer: {e}")
            return False



def main():
    """Main function for command line usage."""
    # Configuration
    # ref_fasta = "hg38_primary.fa"
    ref_fasta = "hg38_primary.fa"
    # primer_fasta = "str_primers.txt"
    # primer_fasta = "dys437.fasta"
    primer_fasta = "dys437.fasta"
    prefix = '2025-08-11_dys437_results_mummer'
    run_alignment = True  # Set to True to run alignment, False to use existing delta file

    # Initialize checker with custom parameters
    checker = MummerChecker(
        tm_threshold=40.0,
        max_amplicon_size=650,
        min_amplicon_size=50,
        region_padding=50
    )

    # Run analysis
    success = checker.analyze_specificity(ref_fasta, primer_fasta, prefix, run_alignment)

    if success:
        logger.info("Analysis completed successfully")
    else:
        logger.error("Analysis failed")
        sys.exit(1)


if __name__ == '__main__':
    main()