"""BAM strand splitting using samtools flag filters.

Replicates the flag logic from biock_split_bam_by_strand.py to split
strand-specific reads from RNA-seq BAM files.
"""

import logging
import os
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
    min_mapq : int
        Minimum mapping quality filter.
    primary : bool
        If True, only keep primary alignments (exclude flag 256).
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

    fwd_bam = f"{prefix}.fwd.bam"
    rev_bam = f"{prefix}.rev.bam"

    for strand, out_bam in [("fwd", fwd_bam), ("rev", rev_bam)]:
        filters = _STRAND_FILTERS[libtype][strand]

        # Build samtools view command
        parts = ["samtools", "view", "-b", f"-@ {n_jobs}", f"-h"]
        if min_mapq > 0:
            parts.append(f"-q {min_mapq}")
        if primary:
            parts.append("-F 256")
        parts.extend(filters)
        parts.append(bam)
        parts.append(f"-o {out_bam}")

        cmd = " ".join(parts)
        run_cmd(cmd, dry_run=dry_run)

        # Index the output BAM
        run_cmd(f"samtools index -@ {n_jobs} {out_bam}", dry_run=dry_run)

    logger.info("Split %s -> %s, %s", bam, fwd_bam, rev_bam)
    return fwd_bam, rev_bam
