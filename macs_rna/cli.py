"""CLI argument parsing for macs-rna.

Mirrors macs3 callpeak arguments with additional strand-specific options.
"""

import argparse
import sys

from macs_rna import __version__
from macs_rna.callpeak import run_callpeak

# Preset recipes for common experiment types
RECIPES: dict[str, dict] = {
    "m6A-MeRIP": {
        "keep_dup": "all",
    },
    "seCLIP": {
        "keep_dup": 1,
    },
}


def get_callpeak_parser() -> argparse.ArgumentParser:
    """Build argument parser for the callpeak subcommand.

    Returns
    -------
    argparse.ArgumentParser
        Parser with all macs3 callpeak args + strand-specific extras.
    """
    p = argparse.ArgumentParser(
        prog="macs-rna callpeak",
        description="Strand-specific MACS3 peak calling for RNA-seq data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Strand-specific arguments (new) ---
    strand_group = p.add_argument_group("Strand-specific arguments")
    strand_group.add_argument(
        "--libtype",
        required=True,
        choices=("FR", "RF", "F", "R"),
        help="Library strandedness type",
    )
    strand_group.add_argument(
        "--chrom-sizes",
        help="Chromsizes file or shortcut (hg38, mm10) for bdg-to-bigwig conversion",
    )
    strand_group.add_argument(
        "--min-mapq", type=int, default=0,
        help="Minimum mapping quality for strand split",
    )
    strand_group.add_argument(
        "--primary", action="store_true",
        help="Only keep primary alignments during strand split",
    )
    strand_group.add_argument(
        "--n-jobs", type=int, default=8,
        help="Number of parallel threads for samtools",
    )
    strand_group.add_argument(
        "--mode",
        choices=tuple(RECIPES.keys()),
        help="Preset recipe to override defaults",
    )
    strand_group.add_argument(
        "--macs-path", default="macs3",
        help="Path to MACS3 executable",
    )
    strand_group.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing",
    )

    # --- MACS3 callpeak arguments (passthrough) ---
    required_group = p.add_argument_group("Required MACS3 arguments")
    required_group.add_argument(
        "-t", "--treatment", nargs="+", required=True,
        help="Treatment BAM/BED file(s)",
    )
    required_group.add_argument(
        "-c", "--control", nargs="+", default=None,
        help="Control BAM/BED file(s)",
    )
    required_group.add_argument(
        "-n", "--name", required=True,
        help="Name prefix for output files",
    )
    required_group.add_argument(
        "--outdir", required=True,
        help="Output directory",
    )
    required_group.add_argument(
        "-f", "--format",
        choices=("AUTO", "BAM", "BAMPE", "SAM", "BED", "BEDPE"),
        default=None, required=True,
        help="Input file format",
    )
    required_group.add_argument(
        "-g", "--gsize", default="hs",
        help="Effective genome size (hs, mm, ce, dm, or a number)",
    )
    required_group.add_argument(
        "--keep-dup", required=True,
        help="Duplicate handling: 'auto', 'all', or integer",
    )
    required_group.add_argument(
        "--nomodel", action="store_true",
        help="Skip fragment size model building",
    )

    # --- Cutoff (mutually exclusive) ---
    cutoff = p.add_mutually_exclusive_group(required=True)
    cutoff.add_argument(
        "-q", "--qvalue", type=float,
        help="Q-value (FDR) cutoff for peak calling",
    )
    cutoff.add_argument(
        "-p", "--pvalue", type=float,
        help="P-value cutoff for peak calling",
    )

    # --- Optional MACS3 arguments ---
    opt_group = p.add_argument_group("Optional MACS3 arguments")
    opt_group.add_argument(
        "-B", "--bdg", action="store_true",
        help="Save bedGraph files in output",
    )
    opt_group.add_argument(
        "--SPMR", action="store_true",
        help="Save signal per million reads in bedGraph (implies -B)",
    )
    opt_group.add_argument(
        "--scale-to", choices=("small", "large"), default="small",
        help="Scale to small or large library size",
    )
    opt_group.add_argument(
        "--broad", action="store_true",
        help="Call broad peaks",
    )
    opt_group.add_argument(
        "--broad-cutoff", type=float, default=0.1,
        help="Cutoff for broad region linking",
    )
    opt_group.add_argument(
        "--call-summits", action="store_true",
        help="Call summits of peaks (only used in initial macs3 callpeak run)",
    )
    opt_group.add_argument(
        "--shift", type=int, default=None,
        help="Shift size for reads in bp",
    )
    opt_group.add_argument(
        "--extsize", type=int, default=None,
        help="Extension size for reads in bp",
    )
    opt_group.add_argument(
        "--max-gap", type=int, default=None,
        help="Maximum gap between peaks to merge",
    )
    opt_group.add_argument(
        "--min-length", type=int, default=None,
        help="Minimum peak length in bp",
    )
    opt_group.add_argument(
        "--nolambda", action="store_true",
        help="Use fixed background lambda",
    )
    opt_group.add_argument(
        "--slocal", type=int, default=None,
        help="Small local region size for dynamic lambda (bp)",
    )
    opt_group.add_argument(
        "--llocal", type=int, default=None,
        help="Large local region size for dynamic lambda (bp)",
    )
    opt_group.add_argument(
        "--fe-cutoff", type=float, default=None,
        help="Minimum fold enrichment cutoff",
    )
    opt_group.add_argument(
        "--buffer-size", type=int, default=None,
        help="Buffer size for read storage",
    )
    opt_group.add_argument(
        "--verbose", type=int, default=2,
        help="Verbosity level (0-3)",
    )

    return p


def main() -> None:
    """Main entry point for macs-rna CLI."""
    # Top-level parser with subcommands
    parser = argparse.ArgumentParser(
        prog="macs-rna",
        description="Strand-specific MACS3 wrapper for RNA-seq peak calling",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Register callpeak subcommand
    callpeak_parser = get_callpeak_parser()
    subparsers.add_parser(
        "callpeak",
        parents=[callpeak_parser],
        add_help=False,
        help="Strand-specific peak calling",
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "callpeak":
        _apply_recipe(args)
        run_callpeak(args)


def _apply_recipe(args: argparse.Namespace) -> None:
    """Apply preset recipe overrides to parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments. Modified in place.
    """
    if args.mode is None:
        return

    recipe = RECIPES[args.mode]
    for key, value in recipe.items():
        current = getattr(args, key.replace("-", "_"), None)
        if current != value:
            print(f"[recipe:{args.mode}] Overriding --{key}: {current} -> {value}")
            setattr(args, key.replace("-", "_"), value)
