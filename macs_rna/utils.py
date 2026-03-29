"""Utility functions for macs-rna."""

import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Union

import yaml

logger = logging.getLogger("macs_rna")

# Path to seedat resources config
_RESOURCES_YAML = Path("~/.config/seedat/resources.yaml").expanduser()


def setup_logging(verbose: int = 2) -> None:
    """Configure logging for macs-rna.

    Parameters
    ----------
    verbose : int
        Verbosity level: 0=critical, 1=warning, 2=info, 3=debug.
    """
    level_map = {0: logging.CRITICAL, 1: logging.WARNING, 2: logging.INFO, 3: logging.DEBUG}
    level = level_map.get(verbose, logging.INFO)
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def run_cmd(
    cmd: Union[str, list[str]],
    dry_run: bool = False,
    check: bool = True,
    capture: bool = False,
) -> Optional[subprocess.CompletedProcess]:
    """Run a shell command via subprocess.

    Parameters
    ----------
    cmd : str or list[str]
        Command to run. If a string, executed via shell.
    dry_run : bool
        If True, log the command without executing.
    check : bool
        If True, raise CalledProcessError on non-zero exit.
    capture : bool
        If True, capture stdout/stderr and return them.

    Returns
    -------
    subprocess.CompletedProcess or None
        Result of the command, or None if dry_run.
    """
    cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
    logger.info("CMD: %s", cmd_str)

    if dry_run:
        return None

    result = subprocess.run(
        cmd_str,
        shell=True,
        check=check,
        capture_output=capture,
        text=True if capture else None,
    )
    return result


def count_reads(bam: str) -> int:
    """Count mapped reads in a BAM file using samtools.

    Parameters
    ----------
    bam : str
        Path to BAM file.

    Returns
    -------
    int
        Number of mapped reads (excludes unmapped, flag 4).
    """
    result = run_cmd(f"samtools view -c -F 4 {bam}", capture=True)
    return int(result.stdout.strip())


def resolve_chromsizes(name_or_path: str) -> str:
    """Resolve a chromsizes file path from a shortcut or direct path.

    Supports shortcuts like 'hg38', 'mm10' via seedat resources.yaml,
    or accepts a direct file path.

    Parameters
    ----------
    name_or_path : str
        Either a genome build shortcut (e.g., 'hg38') or a file path.

    Returns
    -------
    str
        Resolved absolute path to the chromsizes file.

    Raises
    ------
    FileNotFoundError
        If the resolved path does not exist.
    """
    # If it's already a file, return it
    expanded = os.path.expanduser(os.path.expandvars(name_or_path))
    if os.path.isfile(expanded):
        return expanded

    # Try resolving via seedat resources.yaml
    if _RESOURCES_YAML.is_file():
        with open(_RESOURCES_YAML) as f:
            resources = yaml.safe_load(f)
        if name_or_path in resources and "chromsize" in resources[name_or_path]:
            resolved = os.path.expanduser(resources[name_or_path]["chromsize"])
            if os.path.isfile(resolved):
                return resolved
            raise FileNotFoundError(
                f"Chromsizes file from resources.yaml not found: {resolved}"
            )

    raise FileNotFoundError(
        f"Cannot resolve chromsizes: '{name_or_path}' is not a file and "
        f"not found in {_RESOURCES_YAML}"
    )


def bdg_to_bigwig(
    bdg: str,
    bw: str,
    chrom_sizes: str,
) -> None:
    """Convert a bedGraph file to bigWig format.

    Sorts the bedGraph, converts via bedGraphToBigWig, then removes the
    original bdg file.

    Parameters
    ----------
    bdg : str
        Input bedGraph file path.
    bw : str
        Output bigWig file path.
    chrom_sizes : str
        Path to chromosome sizes file.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        sorted_bdg = os.path.join(tmpdir, "sorted.bdg")
        run_cmd(f"sort -k1,1 -k2,2n {bdg} > {sorted_bdg}")
        run_cmd(f"bedGraphToBigWig {sorted_bdg} {chrom_sizes} {bw}")

    # Remove original bdg after successful conversion
    os.remove(bdg)
    logger.info("Converted %s -> %s (bdg removed)", bdg, bw)
