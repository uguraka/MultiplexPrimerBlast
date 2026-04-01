import argparse
import sys
from multiplex_checker.MultiplexCheckerBlast import BlastChecker
from multiplex_checker.MultiplexCheckerMummer import MummerChecker

def main():
    parser = argparse.ArgumentParser(description="Multiplex PCR Primer Specificity Checker")
    parser.add_argument("--ref", required=True, help="Path to reference FASTA file")
    parser.add_argument("--primers", required=True, help="Path to primer FASTA file")
    parser.add_argument("--prefix", default="primer_alignment", help="Output prefix")
    parser.add_argument("--tool", choices=["blast", "mummer"], default="blast", help="Alignment tool to use (blast or mummer)")
    parser.add_argument("--tm-threshold", type=float, default=30.0, help="Minimum melting temperature (Tm) threshold")
    parser.add_argument("--max-amplicon", type=int, default=1000, help="Maximum amplicon size")
    parser.add_argument("--min-amplicon", type=int, default=50, help="Minimum amplicon size")
    parser.add_argument("--skip-alignment", action="store_true", help="Skip alignment and use existing results")

    args = parser.parse_args()

    if args.tool == "blast":
        checker = BlastChecker(
            tm_threshold=args.tm_threshold,
            max_amplicon_size=args.max_amplicon,
            min_amplicon_size=args.min_amplicon,
        )
    elif args.tool == "mummer":
        checker = MummerChecker(
            tm_threshold=args.tm_threshold,
            max_amplicon_size=args.max_amplicon,
            min_amplicon_size=args.min_amplicon,
        )
    else:
        print(f"Unknown tool: {args.tool}")
        sys.exit(1)

    success = checker.analyze_specificity(
        ref_fasta=args.ref,
        primer_fasta=args.primers,
        prefix=args.prefix,
        run_alignment=not args.skip_alignment
    )

    if success:
        print("Analysis completed successfully.")
    else:
        print("Analysis failed.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
