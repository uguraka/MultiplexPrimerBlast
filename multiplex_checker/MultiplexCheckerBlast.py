from typing import List, Dict
import subprocess
import os
from multiplex_checker.MultiplexSpecifityChecker import PCRSpecificityChecker
import logging
import sys
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)


class BlastChecker(PCRSpecificityChecker):
    def __init__(self, tm_threshold=30.0, max_amplicon_size=650, min_amplicon_size=50, region_padding=50):
        super().__init__(tm_threshold=tm_threshold, max_amplicon_size=max_amplicon_size, min_amplicon_size=min_amplicon_size, region_padding=region_padding)
        self.num_threads = 25
        self.word_size = 12


    def _db_files_exist(self, db_prefix: str) -> bool:
        required_exts = [".nhr", ".nin", ".nsq"]
        return all(os.path.exists(db_prefix + ext) for ext in required_exts)



    def _create_blast_db(self, ref_fasta: str, db_prefix: str):
        # if not shutil.which("makeblastdb"):
        #     raise RuntimeError("makeblastdb not found in PATH")
        cmd = ["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", db_prefix]
        logger.info(f"Creating BLAST DB: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)



    def run_alignment(self, ref_fasta: str, primer_fasta: str, prefix: str = 'primer_alignment', db_prefix: str = None) -> bool:
        if db_prefix is None:
            db_prefix = ref_fasta  # reuse same path for DB prefix
            logger.info(f"No db_prefix given. Using {ref_fasta} as DB prefix.")
            if not self._db_files_exist(db_prefix):
                logger.info(f"BLAST DB not found. Creating...")
                self._create_blast_db(ref_fasta, db_prefix)
            else:
                logger.info(f"BLAST DB already exists at {db_prefix}")
        else:
            if not self._db_files_exist(db_prefix):
                logger.info(f"BLAST DB not found at {db_prefix}. Creating...")
                self._create_blast_db(ref_fasta, db_prefix)
            else:
                logger.info(f"Using existing BLAST DB at {db_prefix}")
        # build blastn command
        cmd = [
            "blastn",
            "-task", "blastn-short",
            "-query", primer_fasta,
            "-db", db_prefix,
            "-out", prefix,
            "-word_size", str(self.word_size),
            "-outfmt", "7 qseqid sseqid qstart qend sstart send length pident evalue sstrand",
            "-num_threads", str(self.num_threads)
        ]

        try:
            logger.info(f"Running blastn command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            logger.info(f"blastn completed successfully.")
            return True

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running blastn: {e}")
            return False

    def parse_alignment_results(self, prefix: str) -> List[Dict]:
        if not os.path.exists(prefix):
            logger.error(f"Alignment result file not found: {prefix}")
            return []
        try:
            df = pd.read_csv(
                prefix,
                sep='\t',
                comment='#',
                header=None,
                names=[
                    "qseqid", "sseqid", "qstart", "qend",
                    "sstart", "send", "length", "pident",
                    "evalue", "sstrand"
                ]
            )

            # todo Specify dtype option on import or set low_memory=False.
            df['sseqid'] = df['sseqid'].astype(str)
        except Exception as e:
            logger.error(f"Failed to read alignment result: {e}")
            return []
        if df.empty:
            logger.warning(f"No alignments found in {prefix}")
            return []

        df = df[df['pident'] >= 80.0]  # todo make this a parameter

        parsed_results = []
        for row in df.itertuples(index=False):
            strand = 'forward' if row.sstrand == 'plus' else 'reverse'

            result = {
                'primer_name': row.qseqid,
                'chromosome': str(row.sseqid),
                'ref_start': int(min(row.sstart, row.send)),
                'ref_end': int(max(row.sstart, row.send)),
                'primer_start': int(min(row.qstart, row.qend)),
                'primer_end': int(max(row.qstart, row.qend)),
                'strand': strand,
                'is_forward_strand': (strand == 'forward')
            }
            parsed_results.append(result)

        logger.info(f"Parsed {len(parsed_results)} alignments from {prefix}")
        return parsed_results


def main(primer_fasta):
    """Main function for command line usage."""
    # Configuration
    # ref_fasta = "hg38_primary.fa"
    ref_fasta = "/home/burakdemiralay/MultiplexSpecifity/hg38_primary.fa"
    # primer_fasta = "str_primers.txt"
    # primer_fasta = "dys437.fasta"
    # primer_fasta = "/home/burakdemiralay/MultiplexSpecifity/str_primers.txt"
    prefix = '2025-12-09_trisomy_results_blast' ##prefix = '2025-08-11_dys437_results_blast'
    run_alignment = True  # Set to True to run alignment, False to use existing delta file wowzie #todo

    # Initialize checker with custom parameters
    checker = BlastChecker(
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

    primer_fasta = "/home/burakdemiralay/MultiplexSpeci fity/2025-12-09_trisomy_primer_probe_mchecker.fasta" #primer_fasta = "/home/burakdemiralay/MultiplexSpecifity/dys437.fasta"

    main(primer_fasta)





