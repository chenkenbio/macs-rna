"""Strand-specific MACS3 callpeak orchestration.

Workflow:
1. Split treatment (and control) BAMs by strand
2. Count total reads across both strands for global scaling
3. For each strand:
   a. Run macs3 callpeak -B to generate bedGraphs
   b. Re-score with global scaling (bdgopt multiply + bdgcmp, or
      bdgcmp -S when --SPMR)
   c. Convert p-values to q-values if needed (macs3 bdgopt -m p2q)
   d. Call peaks with macs3 bdgpeakcall / bdgbroadcall
4. Set strand in peak files, combine, convert bdg to bigwig
"""

import argparse
import logging
import math
import os
import shlex
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

    # Warn about format pitfalls for RNA-seq
    if args.format == "AUTO":
        logger.warning(
            "Using AUTO format detection. For paired-end RNA-seq, MACS3 may "
            "auto-detect BAMPE which uses fragment coordinates spanning introns. "
            "This produces incorrect pileup for spliced reads. "
            "Please use -f BAM with --nomodel --extsize instead."
        )

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
        # Always clean up temp directory (including dry-run)
        if os.path.exists(tmpdir):
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
    use_spmr = args.SPMR

    # --- Step 1: Split BAMs by strand ---
    logger.info("=== Step 1/5: Splitting BAMs by strand ===")

    treat_prefix = os.path.join(tmpdir, f"{name}_treat")
    treat_fwd, treat_rev = split_bam_by_strand(
        bam=args.treatment,
        libtype=args.libtype,
        prefix=treat_prefix,
        read=args.read,
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
            bam=args.control,
            libtype=args.libtype,
            prefix=ctrl_prefix,
            read=args.read,
            min_mapq=args.min_mapq,
            primary=args.primary,
            n_jobs=args.n_jobs,
            dry_run=dry_run,
        )

    # --- Step 2: Count total and per-strand reads ---
    logger.info("=== Step 2/5: Counting reads ===")

    if not dry_run:
        strand_treat_counts = {
            "fwd": count_reads(treat_fwd),
            "rev": count_reads(treat_rev),
        }
        # Total = sum of strand counts (already filtered by --read)
        total_treat = strand_treat_counts["fwd"] + strand_treat_counts["rev"]
        total_treat_M = total_treat / 1e6
        logger.info("Total treatment reads (after --read %s filter): %d (%.3f M)",
                     args.read, total_treat, total_treat_M)
        for s, c in strand_treat_counts.items():
            logger.info("  %s strand: %d (%.3f M)", s, c, c / 1e6)

        if total_treat == 0:
            raise RuntimeError(
                "Zero treatment reads after strand splitting and --read filter. "
                "Check that the BAM is not empty, the --libtype matches the "
                "library, and --read is appropriate for the data."
            )

        if args.control is not None:
            ctrl_fwd_n = count_reads(ctrl_fwd)
            ctrl_rev_n = count_reads(ctrl_rev)
            total_ctrl = ctrl_fwd_n + ctrl_rev_n
            logger.info("Total control reads (after --read %s filter): %d (%.3f M)",
                         args.read, total_ctrl, total_ctrl / 1e6)
            if total_ctrl == 0:
                raise RuntimeError("Zero control reads after strand splitting and --read filter.")
    else:
        total_treat = 0
        total_treat_M = 1.0  # placeholder for dry-run
        strand_treat_counts = {"fwd": 0, "rev": 0}
        logger.info("Dry-run: using placeholder values")

    # --- Step 3: Estimate fragment size if needed ---
    extsize = args.extsize
    if not args.nomodel and extsize is None:
        logger.info("=== Step 3/5: Estimating fragment size from full BAM ===")
        logger.info(
            "Note: predictd runs on the full BAM (all reads, both strands). "
            "Use --nomodel --extsize to specify fragment size explicitly."
        )
        extsize = _estimate_fragment_size(args, dry_run)

    # --- Step 4: Per-strand peak calling ---
    logger.info("=== Step 4/5: Per-strand peak calling ===")

    strand_configs = [
        ("fwd", "+", treat_fwd, ctrl_fwd),
        ("rev", "-", treat_rev, ctrl_rev),
    ]

    peak_files = []
    for strand_name, strand_char, t_bam, c_bam in strand_configs:
        strand_reads = strand_treat_counts[strand_name]
        if not dry_run and strand_reads == 0:
            logger.warning(
                "Skipping %s strand: 0 treatment reads after filtering", strand_name,
            )
            continue

        logger.info("--- Processing %s strand ---", strand_name)
        strand_dir = os.path.join(tmpdir, strand_name)
        os.makedirs(strand_dir, exist_ok=True)

        strand_prefix = f"{name}_{strand_name}"

        # 4a. Run macs3 callpeak to generate bedGraphs
        _run_macs3_callpeak(
            args=args,
            treatment=t_bam,
            control=c_bam,
            name=strand_prefix,
            outdir=strand_dir,
            extsize=extsize,
            dry_run=dry_run,
        )

        # 4b-4d. Re-score and call peaks with global scaling
        peak_file = _rescale_and_call_peaks(
            strand_dir=strand_dir,
            strand_prefix=strand_prefix,
            outdir=outdir,
            total_treat=total_treat,
            total_treat_M=total_treat_M,
            strand_treat=strand_reads,
            use_spmr=use_spmr,
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
            is_spmr=use_spmr,
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
        f" -i {shlex.quote(args.treatment)}"
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
            for token in line.split():
                try:
                    d = int(token)
                    if 20 < d < 2000:  # sanity check
                        logger.info("Predicted fragment size d = %d", d)
                        return d
                except ValueError:
                    continue

    logger.warning(
        "Could not parse fragment size from predictd output. "
        "Falling back to extsize=200. predictd output:\n%s",
        output.strip() or "(empty)",
    )
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
    """Run macs3 callpeak with -B to generate bedGraphs for re-scoring.

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
    # Convert strand-split BAM to BED — MACS3 >=3.0.3 silently drops PE reads
    # with -f BAM when the BAM has paired flags (bit 1 set).
    treatment_bed = treatment.replace(".bam", ".bed")
    run_cmd(f"bedtools bamtobed -i {shlex.quote(treatment)} > {shlex.quote(treatment_bed)}", dry_run=dry_run)
    control_bed = None
    if control is not None:
        control_bed = control.replace(".bam", ".bed")
        run_cmd(f"bedtools bamtobed -i {shlex.quote(control)} > {shlex.quote(control_bed)}", dry_run=dry_run)

    cmd_parts = [
        args.macs_path, "callpeak",
        f"-t {shlex.quote(treatment_bed)}",
        "-f BED",
        f"-g {args.gsize}",
        f"-n {shlex.quote(name)}",
        f"--outdir {shlex.quote(outdir)}",
        f"--keep-dup {args.keep_dup}",
        f"--scale-to {args.scale_to}",
        "-B",  # always generate bedGraphs for re-scoring
    ]

    # Only add --SPMR if user explicitly requested it
    if args.SPMR:
        cmd_parts.append("--SPMR")

    if control_bed is not None:
        cmd_parts.append(f"-c {shlex.quote(control_bed)}")

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
    total_treat: int,
    total_treat_M: float,
    strand_treat: int,
    use_spmr: bool,
    args: argparse.Namespace,
    dry_run: bool,
) -> str:
    """Re-score bedGraphs with global scaling and call peaks.

    Two modes:
    - Raw (default): scale treat/ctrl bedGraphs by total_treat/strand_treat
      via bdgopt multiply, then bdgcmp with default S=1.0.
    - SPMR (--SPMR): use bdgcmp -S total_treat_M directly on SPMR bedGraphs.

    Parameters
    ----------
    strand_dir : str
        Directory containing macs3 callpeak output for this strand.
    strand_prefix : str
        Name prefix (e.g., 'experiment_fwd').
    outdir : str
        Final output directory for peak files.
    total_treat : int
        Total treatment reads (both strands).
    total_treat_M : float
        Total treatment reads in millions.
    strand_treat : int
        Treatment reads in this strand.
    use_spmr : bool
        Whether --SPMR was passed by the user.
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

    if use_spmr:
        # SPMR mode: bedGraphs are per-million-reads normalized,
        # use bdgcmp -S to set effective library size to total reads
        bdgcmp_treat = treat_bdg
        bdgcmp_ctrl = ctrl_bdg
        scaling_flag = f"-S {total_treat_M:.6f}"
    else:
        # Raw mode: scale both bedGraphs by total/strand to equalize
        # effective library size across strands
        if strand_treat > 0 and total_treat > 0:
            scale_factor = total_treat / strand_treat
        else:
            scale_factor = 1.0
        logger.info(
            "Scaling factor for %s: %.4f (total=%d, strand=%d)",
            strand_prefix, scale_factor, total_treat, strand_treat,
        )

        bdgcmp_treat = os.path.join(strand_dir, f"{strand_prefix}_treat_scaled.bdg")
        bdgcmp_ctrl = os.path.join(strand_dir, f"{strand_prefix}_ctrl_scaled.bdg")

        run_cmd(
            f"{shlex.quote(args.macs_path)} bdgopt"
            f" -i {shlex.quote(treat_bdg)}"
            f" -m multiply"
            f" -p {scale_factor:.6f}"
            f" -o {shlex.quote(bdgcmp_treat)}",
            dry_run=dry_run,
        )
        run_cmd(
            f"{shlex.quote(args.macs_path)} bdgopt"
            f" -i {shlex.quote(ctrl_bdg)}"
            f" -m multiply"
            f" -p {scale_factor:.6f}"
            f" -o {shlex.quote(bdgcmp_ctrl)}",
            dry_run=dry_run,
        )
        scaling_flag = ""  # default S=1.0

    # Re-score with bdgcmp
    ppois_bdg = os.path.join(strand_dir, f"{strand_prefix}_ppois.bdg")
    bdgcmp_cmd = (
        f"{shlex.quote(args.macs_path)} bdgcmp"
        f" -t {shlex.quote(bdgcmp_treat)}"
        f" -c {shlex.quote(bdgcmp_ctrl)}"
    )
    if scaling_flag:
        bdgcmp_cmd += f" {scaling_flag}"
    bdgcmp_cmd += f" -m ppois -o {shlex.quote(ppois_bdg)}"
    run_cmd(bdgcmp_cmd, dry_run=dry_run)

    # Determine score track and cutoff for peak calling
    if args.pvalue is not None:
        # Use p-value scores directly
        score_bdg = ppois_bdg
        cutoff = -math.log10(args.pvalue)
    else:
        # Convert to q-values for FDR cutoff
        qval = args.qvalue if args.qvalue is not None else 0.05
        qvalue_bdg = os.path.join(strand_dir, f"{strand_prefix}_qvalue.bdg")
        run_cmd(
            f"{shlex.quote(args.macs_path)} bdgopt"
            f" -i {shlex.quote(ppois_bdg)}"
            f" -m p2q"
            f" -o {shlex.quote(qvalue_bdg)}",
            dry_run=dry_run,
        )
        score_bdg = qvalue_bdg
        cutoff = -math.log10(qval)

    # Call peaks from score track
    if args.broad:
        peak_suffix = "broadPeak"
        peak_file = os.path.join(outdir, f"{strand_prefix}_peaks.{peak_suffix}")
        cmd_parts = [
            f"{args.macs_path} bdgbroadcall",
            f"-i {shlex.quote(score_bdg)}",
            f"-c {cutoff:.4f}",
            f"-C {-math.log10(args.broad_cutoff):.4f}",
            f"-o {shlex.quote(peak_file)}",
        ]
        if args.min_length is not None:
            cmd_parts.append(f"-l {args.min_length}")
        if args.max_gap is not None:
            cmd_parts.append(f"-g {args.max_gap}")
    else:
        peak_suffix = "narrowPeak"
        peak_file = os.path.join(outdir, f"{strand_prefix}_peaks.{peak_suffix}")
        cmd_parts = [
            f"{args.macs_path} bdgpeakcall",
            f"-i {shlex.quote(score_bdg)}",
            f"-c {cutoff:.4f}",
            f"-o {shlex.quote(peak_file)}",
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
    is_spmr: bool,
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
    is_spmr : bool
        Whether the bedGraphs are SPMR-normalized.
    dry_run : bool
        If True, skip file operations.
    """
    if dry_run:
        return

    bdg_files = {
        "treat_pileup": os.path.join(strand_dir, f"{strand_prefix}_treat_pileup.bdg"),
        "control_lambda": os.path.join(strand_dir, f"{strand_prefix}_control_lambda.bdg"),
    }

    # Tag for filename: explicitly label SPMR vs raw
    tag = ".SPMR" if is_spmr else ""

    for label, bdg_path in bdg_files.items():
        if not os.path.exists(bdg_path):
            continue

        if user_requested_bdg:
            # Rename with normalization tag for clarity
            tagged_name = f"{strand_prefix}_{label}{tag}.bdg"
            final_bdg = os.path.join(outdir, tagged_name)
            shutil.move(bdg_path, final_bdg)

            if chrom_sizes is not None:
                bw_path = final_bdg.removesuffix(".bdg") + ".bw"
                try:
                    bdg_to_bigwig(final_bdg, bw_path, chrom_sizes)
                    os.remove(final_bdg)
                    logger.info("Removed %s (bigwig created)", os.path.basename(final_bdg))
                except Exception as e:
                    logger.warning("bedGraph-to-bigWig conversion failed for %s: %s (keeping bdg)", os.path.basename(final_bdg), e)
        # else: bdg stays in tmpdir and is cleaned up automatically


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
    run_cmd(f"sort -k1,1 -k2,2n -o {shlex.quote(combined)} {shlex.quote(combined)}")
    logger.info("Combined peaks written to: %s", combined)
