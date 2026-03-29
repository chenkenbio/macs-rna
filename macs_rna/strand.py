"""BAM strand splitting using samtools flag filters.

Replicates the flag logic from biock_split_bam_by_strand.py to split
strand-specific reads from RNA-seq BAM files.
"""

import logging
import os
import shlex
from typing import Literal, Optional

from macs_rna.utils import run_cmd

logger = logging.getLogger("macs_rna")

# SAM flag filter expressions for strand splitting.
# Paired-end FR: read1-forward/read2-reverse is the standard Illumina layout.
#   Forward strand: (flag & 99)==99 (read1 fwd, mate rev)
#                   OR (flag & 147)==147 (read2 rev, mate fwd)
#   Reverse strand: (flag & 83)==83 (read1 rev, mate fwd)
#                   OR (flag & 163)==163 (read2 fwd, mate rev)
#
# Flag breakdowns:
#   99  = 1+2+32+64       = paired, proper, mate_rev, read1
#   147 = 1+2+16+128      = paired, proper, rev, read2
#   83  = 1+2+16+64       = paired, proper, rev, read1
#   163 = 1+2+32+128      = paired, proper, mate_rev, read2

_STRAND_FILTERS: dict[str, dict[str, list[str]]] = {
    "FR": {
        "fwd": [
            "-e '(flag & 99) == 99 || (flag & 147) == 147'",
        ],
        "rev": [
            "-e '(flag & 83) == 83 || (flag & 163) == 163'",
        ],
    },
    "RF": {
        "fwd": [
            "-e '(flag & 83) == 83 || (flag & 163) == 163'",
        ],
        "rev": [
            "-e '(flag & 99) == 99 || (flag & 147) == 147'",
        ],
    },
    "F": {
        "fwd": ["-F 16"],
        "rev": ["-f 16"],
    },
    "R": {
        "fwd": ["-f 16"],
        "rev": ["-F 16"],
    },
}


def split_bam_by_strand(
    bam: str,
    libtype: Literal["FR", "RF", "F", "R"],
    prefix: str,
    read: str = "1",
    min_mapq: int = 0,
    primary: bool = False,
    n_jobs: int = 8,
    dry_run: bool = False,
) -> tuple[str, str]:
    """Split a BAM file into forward and reverse strand BAMs.

    Parameters
    ----------
    bam : str
        Input BAM file path.
    libtype : {'FR', 'RF', 'F', 'R'}
        Library strandedness type.
    prefix : str
        Output prefix. Produces {prefix}.fwd.bam and {prefix}.rev.bam.
    read : {'1', '2', 'both'}
        Which read to keep for PE data. '1' = R1 only (flag 64),
        '2' = R2 only (flag 128), 'both' = all reads.
        Ignored for SE libtypes (F/R).
    min_mapq : int
        Minimum mapping quality filter.
    primary : bool
        If True, only keep primary alignments (exclude secondary 0x100
        and supplementary 0x800 alignments).
    n_jobs : int
        Number of threads for samtools.
    dry_run : bool
        If True, print commands without executing.

    Returns
    -------
    tuple[str, str]
        Paths to (forward_strand.bam, reverse_strand.bam).
    """
    if libtype not in _STRAND_FILTERS:
        raise ValueError(f"Unknown libtype: {libtype}. Must be one of {list(_STRAND_FILTERS.keys())}")

    is_pe = libtype in ("FR", "RF")

    if not is_pe and read != "both":
        logger.warning(
            "--read %s is ignored for single-end libtype %s", read, libtype,
        )

    fwd_bam = f"{prefix}.fwd.bam"
    rev_bam = f"{prefix}.rev.bam"

    for strand, out_bam in [("fwd", fwd_bam), ("rev", rev_bam)]:
        filters = _STRAND_FILTERS[libtype][strand]

        # Build samtools view command
        parts = ["samtools", "view", "-b", f"-@ {n_jobs}", "-h"]
        if min_mapq > 0:
            parts.append(f"-q {min_mapq}")
        # Always exclude unmapped (0x4) and supplementary (0x800).
        # Supplementary alignments are chimeric fragments of the same
        # read and inflate counts if included.
        exclude_flags = 0x4 | 0x800  # 2052
        if primary:
            exclude_flags |= 0x100  # also exclude secondary → 2308
        parts.append(f"-F {exclude_flags}")
        # PE read selection: -f 64 (R1) or -f 128 (R2)
        if is_pe and read == "1":
            parts.append("-f 64")
        elif is_pe and read == "2":
            parts.append("-f 128")
        parts.extend(filters)
        parts.append(shlex.quote(bam))
        parts.append(f"-o {shlex.quote(out_bam)}")

        cmd = " ".join(parts)
        run_cmd(cmd, dry_run=dry_run)

        # Index the output BAM
        run_cmd(f"samtools index -@ {n_jobs} {shlex.quote(out_bam)}", dry_run=dry_run)

    read_label = f" (R{read})" if is_pe and read != "both" else ""
    logger.info("Split %s -> %s, %s%s", bam, fwd_bam, rev_bam, read_label)
    return fwd_bam, rev_bam
