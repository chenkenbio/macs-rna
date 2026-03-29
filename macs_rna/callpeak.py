"""Strand-specific MACS3 callpeak orchestration.

Workflow:
1. Split treatment (and control) BAMs by strand
2. Count total reads across both strands for global scaling
3. For each strand:
   a. Run macs3 callpeak --SPMR -B to generate normalized bedGraphs
   b. Re-score with macs3 bdgcmp using total-reads scaling factor
   c. Convert p-values to q-values if needed (macs3 bdgopt -m p2q)
   d. Call peaks with macs3 bdgpeakcall / bdgbroadcall
4. Set strand in peak files, combine, convert bdg to bigwig
"""

import argparse
import logging
import math
import os
import shutil
import tempfile
from typing import Optional

from macs_rna.strand import split_bam_by_strand
from macs_rna.utils import bdg_to_bigwig, count_reads, resolve_chromsizes, run_cmd

logger = logging.getLogger("macs_rna")


def run_callpeak(args: argparse.Namespace) -> None:
    """Run the full strand-specific callpeak pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments from cli.py.
    """
    from macs_rna.utils import setup_logging

    setup_logging(args.verbose)

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    name = args.name
    dry_run = args.dry_run

    # Resolve chromsizes if provided
    chrom_sizes: Optional[str] = None
    if args.chrom_sizes is not None:
        chrom_sizes = resolve_chromsizes(args.chrom_sizes)
        logger.info("Resolved chromsizes: %s", chrom_sizes)

    # Create a temp directory for intermediate files
    tmpdir = tempfile.mkdtemp(prefix="macs_rna_", dir=outdir)
    logger.info("Temp directory: %s", tmpdir)

    try:
        _run_pipeline(args, outdir, name, tmpdir, chrom_sizes, dry_run)
    finally:
        # Clean up temp directory
        if not dry_run and os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
            logger.info("Cleaned up temp directory: %s", tmpdir)


def _run_pipeline(
    args: argparse.Namespace,
    outdir: str,
    name: str,
    tmpdir: str,
    chrom_sizes: Optional[str],
    dry_run: bool,
) -> None:
    """Internal pipeline implementation.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    outdir : str
        Final output directory.
    name : str
        Experiment name prefix.
    tmpdir : str
        Temporary directory for intermediate files.
    chrom_sizes : str or None
        Resolved chromsizes path, or None.
    dry_run : bool
        If True, print commands without executing.
    """
    # --- Step 1: Split BAMs by strand ---
    logger.info("=== Step 1/5: Splitting BAMs by strand ===")

    treat_prefix = os.path.join(tmpdir, f"{name}_treat")
    treat_fwd, treat_rev = split_bam_by_strand(
        bam=args.treatment[0],
        libtype=args.libtype,
        prefix=treat_prefix,
        min_mapq=args.min_mapq,
        primary=args.primary,
        n_jobs=args.n_jobs,
        dry_run=dry_run,
    )

    ctrl_fwd: Optional[str] = None
    ctrl_rev: Optional[str] = None
    if args.control is not None:
        ctrl_prefix = os.path.join(tmpdir, f"{name}_ctrl")
        ctrl_fwd, ctrl_rev = split_bam_by_strand(
            bam=args.control[0],
            libtype=args.libtype,
            prefix=ctrl_prefix,
            min_mapq=args.min_mapq,
            primary=args.primary,
            n_jobs=args.n_jobs,
            dry_run=dry_run,
        )

    # --- Step 2: Count total reads (before splitting) ---
    logger.info("=== Step 2/5: Counting total reads ===")

    if not dry_run:
        total_treat = count_reads(args.treatment[0])
        total_treat_M = total_treat / 1e6
        logger.info("Total treatment reads: %d (%.3f M)", total_treat, total_treat_M)

        if args.control is not None:
            total_ctrl = count_reads(args.control[0])
            logger.info("Total control reads: %d (%.3f M)", total_ctrl, total_ctrl / 1e6)
    else:
        total_treat_M = 1.0  # placeholder for dry-run
        logger.info("Dry-run: using placeholder total_treat_M=1.0")

    # --- Step 3: Estimate fragment size if needed ---
    extsize = args.extsize
    if not args.nomodel and extsize is None:
        logger.info("=== Step 3/5: Estimating fragment size from full BAM ===")
        extsize = _estimate_fragment_size(args, dry_run)

    # --- Step 4: Per-strand peak calling ---
    logger.info("=== Step 4/5: Per-strand peak calling ===")

    strand_configs = [
        ("fwd", "+", treat_fwd, ctrl_fwd),
        ("rev", "-", treat_rev, ctrl_rev),
    ]

    peak_files = []
    for strand_name, strand_char, t_bam, c_bam in strand_configs:
        logger.info("--- Processing %s strand ---", strand_name)
        strand_dir = os.path.join(tmpdir, strand_name)
        os.makedirs(strand_dir, exist_ok=True)

        strand_prefix = f"{name}_{strand_name}"

        # 3a. Run macs3 callpeak to generate SPMR bedGraphs
        _run_macs3_callpeak(
            args=args,
            treatment=t_bam,
            control=c_bam,
            name=strand_prefix,
            outdir=strand_dir,
            extsize=extsize,
            dry_run=dry_run,
        )

        # 3b-3d. Re-score and call peaks with global scaling
        peak_file = _rescale_and_call_peaks(
            strand_dir=strand_dir,
            strand_prefix=strand_prefix,
            outdir=outdir,
            total_treat_M=total_treat_M,
            args=args,
            dry_run=dry_run,
        )
        peak_files.append((peak_file, strand_char))

        # Convert bdg files to bigwig if chromsizes provided
        _handle_bdg_output(
            strand_dir=strand_dir,
            strand_prefix=strand_prefix,
            outdir=outdir,
            chrom_sizes=chrom_sizes,
            user_requested_bdg=args.bdg or args.SPMR,
            dry_run=dry_run,
        )

    # --- Step 5: Combine strand peaks ---
    logger.info("=== Step 5/5: Combining strand peaks ===")
    _combine_strand_peaks(peak_files, outdir, name, args.broad, dry_run)


def _estimate_fragment_size(args: argparse.Namespace, dry_run: bool) -> Optional[int]:
    """Estimate fragment size using macs3 predictd on the full BAM.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    dry_run : bool
        If True, skip estimation.

    Returns
    -------
    int or None
        Predicted fragment size, or None if dry_run.
    """
    cmd = (
        f"{args.macs_path} predictd"
        f" -i {args.treatment[0]}"
        f" -f {args.format}"
        f" -g {args.gsize}"
    )
    result = run_cmd(cmd, dry_run=dry_run, capture=True, check=False)

    if dry_run or result is None:
        return None

    # Parse predicted d from stderr (MACS3 prints to stderr)
    output = (result.stderr or "") + (result.stdout or "")
    for line in output.splitlines():
        if "predicted fragment length" in line.lower() or "# d =" in line.lower():
            # Extract the number
            for token in line.split():
                try:
                    d = int(token)
                    if 20 < d < 2000:  # sanity check
                        logger.info("Predicted fragment size d = %d", d)
                        return d
                except ValueError:
                    continue

    logger.warning("Could not parse fragment size from predictd output, using default 200")
    return 200


def _run_macs3_callpeak(
    args: argparse.Namespace,
    treatment: str,
    control: Optional[str],
    name: str,
    outdir: str,
    extsize: Optional[int],
    dry_run: bool,
) -> None:
    """Run macs3 callpeak with --SPMR -B to generate normalized bedGraphs.

    The peaks from this run use per-strand scaling and are NOT the final
    output. We only use the bedGraph files for re-scoring.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    treatment : str
        Strand-specific treatment BAM path.
    control : str or None
        Strand-specific control BAM path.
    name : str
        Output name prefix for this strand.
    outdir : str
        Output directory for this strand's macs3 run.
    extsize : int or None
        Fragment extension size (overrides args.extsize if set).
    dry_run : bool
        If True, print command without executing.
    """
    cmd_parts = [
        args.macs_path, "callpeak",
        f"-t {treatment}",
        f"-f {args.format}",
        f"-g {args.gsize}",
        f"-n {name}",
        f"--outdir {outdir}",
        f"--keep-dup {args.keep_dup}",
        f"--scale-to {args.scale_to}",
        "-B --SPMR",  # always generate SPMR bedGraphs
    ]

    if control is not None:
        cmd_parts.append(f"-c {control}")

    # Fragment size handling
    if args.nomodel or extsize is not None:
        cmd_parts.append("--nomodel")
        if extsize is not None:
            cmd_parts.append(f"--extsize {extsize}")
        elif args.extsize is not None:
            cmd_parts.append(f"--extsize {args.extsize}")
    if args.shift is not None:
        cmd_parts.append(f"--shift {args.shift}")

    # Cutoff — macs3 callpeak still needs one for its internal run
    if args.pvalue is not None:
        cmd_parts.append(f"-p {args.pvalue}")
    else:
        qval = args.qvalue if args.qvalue is not None else 0.05
        cmd_parts.append(f"-q {qval}")

    # Optional args
    if args.nolambda:
        cmd_parts.append("--nolambda")
    if args.slocal is not None:
        cmd_parts.append(f"--slocal {args.slocal}")
    if args.llocal is not None:
        cmd_parts.append(f"--llocal {args.llocal}")
    if args.max_gap is not None:
        cmd_parts.append(f"--max-gap {args.max_gap}")
    if args.broad:
        cmd_parts.append("--broad")
    if args.fe_cutoff is not None:
        cmd_parts.append(f"--fe-cutoff {args.fe_cutoff}")
    if args.buffer_size is not None:
        cmd_parts.append(f"--buffer-size {args.buffer_size}")

    cmd = " ".join(cmd_parts)
    run_cmd(cmd, dry_run=dry_run)


def _rescale_and_call_peaks(
    strand_dir: str,
    strand_prefix: str,
    outdir: str,
    total_treat_M: float,
    args: argparse.Namespace,
    dry_run: bool,
) -> str:
    """Re-score bedGraphs with global scaling and call peaks.

    Parameters
    ----------
    strand_dir : str
        Directory containing macs3 callpeak output for this strand.
    strand_prefix : str
        Name prefix (e.g., 'experiment_fwd').
    outdir : str
        Final output directory for peak files.
    total_treat_M : float
        Total treatment reads in millions (the global scaling factor).
    args : argparse.Namespace
        Parsed CLI arguments.
    dry_run : bool
        If True, print commands without executing.

    Returns
    -------
    str
        Path to the final peak file for this strand.
    """
    treat_bdg = os.path.join(strand_dir, f"{strand_prefix}_treat_pileup.bdg")
    ctrl_bdg = os.path.join(strand_dir, f"{strand_prefix}_control_lambda.bdg")

    # 3b. Re-score with global scaling using bdgcmp
    ppois_bdg = os.path.join(strand_dir, f"{strand_prefix}_ppois.bdg")
    run_cmd(
        f"{args.macs_path} bdgcmp"
        f" -t {treat_bdg}"
        f" -c {ctrl_bdg}"
        f" -S {total_treat_M:.6f}"
        f" -m ppois"
        f" -o {ppois_bdg}",
        dry_run=dry_run,
    )

    # 3c. Determine score track and cutoff for peak calling
    if args.pvalue is not None:
        # Use p-value scores directly
        score_bdg = ppois_bdg
        cutoff = -math.log10(args.pvalue)
    else:
        # Convert to q-values for FDR cutoff
        qval = args.qvalue if args.qvalue is not None else 0.05
        qvalue_bdg = os.path.join(strand_dir, f"{strand_prefix}_qvalue.bdg")
        run_cmd(
            f"{args.macs_path} bdgopt"
            f" -i {ppois_bdg}"
            f" -m p2q"
            f" -o {qvalue_bdg}",
            dry_run=dry_run,
        )
        score_bdg = qvalue_bdg
        cutoff = -math.log10(qval)

    # 3d. Call peaks from score track
    if args.broad:
        peak_suffix = "broadPeak"
        peak_file = os.path.join(outdir, f"{strand_prefix}_peaks.{peak_suffix}")
        cmd_parts = [
            f"{args.macs_path} bdgbroadcall",
            f"-i {score_bdg}",
            f"-c {cutoff:.4f}",
            f"-C {-math.log10(args.broad_cutoff):.4f}",
            f"-o {peak_file}",
        ]
        if args.min_length is not None:
            cmd_parts.append(f"-l {args.min_length}")
    else:
        peak_suffix = "narrowPeak"
        peak_file = os.path.join(outdir, f"{strand_prefix}_peaks.{peak_suffix}")
        cmd_parts = [
            f"{args.macs_path} bdgpeakcall",
            f"-i {score_bdg}",
            f"-c {cutoff:.4f}",
            f"-o {peak_file}",
        ]
        if args.min_length is not None:
            cmd_parts.append(f"-l {args.min_length}")
        if args.max_gap is not None:
            cmd_parts.append(f"-g {args.max_gap}")

    run_cmd(" ".join(cmd_parts), dry_run=dry_run)
    return peak_file


def _handle_bdg_output(
    strand_dir: str,
    strand_prefix: str,
    outdir: str,
    chrom_sizes: Optional[str],
    user_requested_bdg: bool,
    dry_run: bool,
) -> None:
    """Convert bedGraph files to bigWig or gzip, then clean up.

    Parameters
    ----------
    strand_dir : str
        Directory containing macs3 callpeak bedGraph output.
    strand_prefix : str
        Name prefix for this strand.
    outdir : str
        Final output directory.
    chrom_sizes : str or None
        Path to chromsizes file, or None to skip bigwig conversion.
    user_requested_bdg : bool
        Whether the user explicitly requested -B/--SPMR output.
    dry_run : bool
        If True, skip file operations.
    """
    if dry_run:
        return

    bdg_files = {
        "treat_pileup": os.path.join(strand_dir, f"{strand_prefix}_treat_pileup.bdg"),
        "control_lambda": os.path.join(strand_dir, f"{strand_prefix}_control_lambda.bdg"),
    }

    for label, bdg_path in bdg_files.items():
        if not os.path.exists(bdg_path):
            continue

        # Rename with .SPMR. to clarify the normalization type
        spmr_bdg = bdg_path.replace(f"_{label}.bdg", f"_{label}.SPMR.bdg")
        os.rename(bdg_path, spmr_bdg)

        if user_requested_bdg:
            # Move to outdir and convert to bigwig if possible
            final_bdg = os.path.join(outdir, os.path.basename(spmr_bdg))
            shutil.move(spmr_bdg, final_bdg)

            if chrom_sizes is not None:
                bw_path = final_bdg.replace(".SPMR.bdg", ".SPMR.bw")
                bdg_to_bigwig(final_bdg, bw_path, chrom_sizes)
        else:
            # User didn't request bdg output; just delete
            os.remove(spmr_bdg)


def _combine_strand_peaks(
    peak_files: list[tuple[str, str]],
    outdir: str,
    name: str,
    broad: bool,
    dry_run: bool,
) -> None:
    """Combine forward and reverse strand peak files.

    Sets column 6 (strand) to '+' or '-' and produces a merged file
    sorted by chromosome and position.

    Parameters
    ----------
    peak_files : list of (path, strand_char)
        List of (peak_file_path, '+' or '-') tuples.
    outdir : str
        Output directory.
    name : str
        Experiment name prefix.
    broad : bool
        Whether broad peaks were called.
    dry_run : bool
        If True, skip file operations.
    """
    if dry_run:
        return

    suffix = "broadPeak" if broad else "narrowPeak"
    combined = os.path.join(outdir, f"{name}_peaks.{suffix}")

    with open(combined, "w") as out_f:
        for peak_path, strand_char in peak_files:
            if not os.path.exists(peak_path):
                logger.warning("Peak file not found: %s", peak_path)
                continue

            with open(peak_path) as in_f:
                for line in in_f:
                    if line.startswith("#") or line.startswith("track"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 6:
                        fields[5] = strand_char
                    out_f.write("\t".join(fields) + "\n")

    # Sort by chr, start
    run_cmd(f"sort -k1,1 -k2,2n -o {combined} {combined}")
    logger.info("Combined peaks written to: %s", combined)
