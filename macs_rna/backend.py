"""Backend selection for routing macs-rna subcommands to macs3-rs or stock macs3.

macs-rna shells out to a small set of ``macs3`` subcommands. The bundled Rust
binary ``macs3-rs`` is a bit-exact drop-in for those subcommands, with an
identical CLI flag surface. This module decides, per subcommand + flag set,
whether the Rust accelerator can be used and where to find it, falling back to
stock ``macs3`` whenever a flag is unsupported or the binary is unavailable.
"""

import logging
import os
import shutil
from importlib import resources
from pathlib import Path
from typing import Iterable, Optional

logger = logging.getLogger("macs_rna")

# Flags macs3-rs supports for each subcommand macs-rna can route to it.
#
# Sets are derived from the Rust ``crates/macs3-cli/src/cmd/*.rs`` ``Args``
# structs (clap). Both the short and long spellings of each option are listed so
# a flag matches regardless of which form macs-rna emits. ``--call-summits`` and
# ``--cutoff-analysis`` are deliberately excluded from ``callpeak`` /
# ``bdgpeakcall`` because macs3-rs does not implement them; their presence forces
# a fallback to stock macs3.
SUPPORTED: dict[str, set[str]] = {
    "callpeak": {
        "-t", "--treatment",
        "-c", "--control",
        "-f", "--format",
        "-g", "--gsize",
        "-s", "--tsize",
        "--keep-dup",
        "-n", "--name",
        "--outdir",
        "-B", "--bdg",
        "--SPMR",
        "--trackline",
        "-q", "--qvalue",
        "-p", "--pvalue",
        "--scale-to",
        "--nolambda",
        "--slocal",
        "--llocal",
        "--max-gap",
        "--min-length",
        "--nomodel",
        "--shift",
        "--extsize",
        "--bw",
        "--d-min",
        "-m", "--mfold",
        "--broad",
        "--broad-cutoff",
        "--fe-cutoff",
        "--buffer-size",
        "--verbose",
    },
    "predictd": {
        "-i", "--ifile",
        "-f", "--format",
        "-g", "--gsize",
        "-s", "--tsize",
        "--bw",
        "--d-min",
        "-m", "--mfold",
        "--rfile",
        "--outdir",
        "--buffer-size",
        "--verbose",
    },
    "pileup": {
        "-i", "--ifile",
        "-o", "--ofile",
        "--outdir",
        "-f", "--format",
        "-B", "--both-direction",
        "--extsize",
        "--buffer-size",
        "--verbose",
    },
    "bdgcmp": {
        "-t", "--tfile",
        "-c", "--cfile",
        "-S", "--scaling-factor",
        "-p", "--pseudocount",
        "-m", "--method",
        "--o-prefix",
        "-o", "--ofile",
        "--outdir",
        "--verbose",
    },
    "bdgopt": {
        "-i", "--ifile",
        "-m", "--method",
        "-p", "--extra-param",
        "-o", "--ofile",
        "--outdir",
        "--verbose",
    },
    "bdgpeakcall": {
        "-i", "--ifile",
        "-c", "--cutoff",
        "-l", "--min-length",
        "-g", "--max-gap",
        "--no-trackline",
        "--o-prefix",
        "-o", "--ofile",
        "--outdir",
        "--verbose",
    },
    "bdgbroadcall": {
        "-i", "--ifile",
        "-c", "--cutoff-peak",
        "-C", "--cutoff-link",
        "-l", "--min-length",
        "-g", "--lvl1-max-gap",
        "-G", "--lvl2-max-gap",
        "--no-trackline",
        "--o-prefix",
        "-o", "--ofile",
        "--outdir",
        "--verbose",
    },
    "bdgdiff": {
        "--t1", "--c1", "--t2", "--c2",
        "-C", "--cutoff",
        "-l", "--min-len",
        "-g", "--max-gap",
        "--d1", "--depth1",
        "--d2", "--depth2",
        "--o-prefix",
        "-o", "--ofile",
        "--outdir",
        "--verbose",
    },
}


def find_accel() -> Optional[str]:
    """Locate the bundled ``macs3-rs`` accelerator binary.

    Discovery order:

    1. ``$MACS_RNA_ACCEL`` -- explicit path override.
    2. Packaged binary ``macs_rna/_accel/macs3-rs`` (installed wheel).
    3. Source-tree build ``rust/macs3-rs/target/release/macs3-rs`` relative to
       the repository (convenient during development).
    4. ``macs3-rs`` on ``$PATH``.

    If ``$MACS_RNA_NO_ACCEL`` is set (to any value), the accelerator is disabled
    and ``None`` is returned (forcing stock macs3).

    Returns
    -------
    str or None
        Absolute path to an executable ``macs3-rs`` binary, or ``None`` if no
        usable binary was found (or the accelerator is disabled).
    """
    if os.environ.get("MACS_RNA_NO_ACCEL"):
        logger.debug("MACS_RNA_NO_ACCEL set; accelerator disabled")
        return None

    # 1. Explicit override.
    explicit = os.environ.get("MACS_RNA_ACCEL")
    if explicit:
        path = os.path.expanduser(os.path.expandvars(explicit))
        if _is_executable(path):
            return os.path.abspath(path)
        logger.warning(
            "MACS_RNA_ACCEL=%s is not an executable file; ignoring", explicit
        )

    # 2. Packaged binary inside the installed wheel.
    try:
        packaged = resources.files("macs_rna").joinpath("_accel", "macs3-rs")
        with resources.as_file(packaged) as p:
            if _is_executable(str(p)):
                return os.path.abspath(str(p))
    except (FileNotFoundError, ModuleNotFoundError, AttributeError):
        pass

    # 3. Source-tree build relative to the repo root (macs_rna/ -> repo root).
    repo_root = Path(__file__).resolve().parent.parent
    src_build = repo_root / "rust" / "macs3-rs" / "target" / "release" / "macs3-rs"
    if _is_executable(str(src_build)):
        return str(src_build)

    # 4. On PATH.
    which = shutil.which("macs3-rs")
    if which:
        return os.path.abspath(which)

    return None


def resolve_backend(
    subcmd: str,
    flags: Iterable[str],
    *,
    prefer: str = "auto",
    macs_path: str = "macs3",
) -> str:
    """Choose the executable to run for a given subcommand and flag set.

    Parameters
    ----------
    subcmd : str
        macs3 subcommand name (e.g. ``"callpeak"``, ``"bdgopt"``).
    flags : Iterable[str]
        The option flags that will be passed (e.g. ``["-i", "-m", "-o"]``).
        Only flag tokens matter; values are ignored.
    prefer : {"auto", "rust", "macs3"}, optional
        Backend preference. ``"macs3"`` always returns ``macs_path``;
        ``"rust"`` and ``"auto"`` use the Rust binary where supported and fall
        back to ``macs_path`` otherwise. With ``"rust"`` an explicit warning is
        logged on fallback. Default ``"auto"``.
    macs_path : str, optional
        Path to the stock macs3 executable used for fallback. Default
        ``"macs3"``.

    Returns
    -------
    str
        Executable path to invoke for this subcommand.
    """
    if prefer == "macs3":
        return macs_path

    flag_set = {f for f in flags if f.startswith("-")}
    supported = SUPPORTED.get(subcmd)
    unsupported_flags = (
        flag_set - supported if supported is not None else flag_set
    )

    if supported is None or unsupported_flags:
        if prefer == "rust":
            reason = (
                f"subcommand '{subcmd}' not supported by macs3-rs"
                if supported is None
                else f"flags {sorted(unsupported_flags)} not supported by macs3-rs"
            )
            logger.warning(
                "--backend rust requested but %s; falling back to %s",
                reason, macs_path,
            )
        return macs_path

    accel = find_accel()
    if accel is None:
        if prefer == "rust":
            logger.warning(
                "--backend rust requested but no macs3-rs binary found; "
                "falling back to %s", macs_path,
            )
        return macs_path

    return accel


def _is_executable(path: str) -> bool:
    """Return True if ``path`` is an existing, executable regular file.

    Parameters
    ----------
    path : str
        Filesystem path to check.

    Returns
    -------
    bool
        Whether the path points to an executable file.
    """
    return os.path.isfile(path) and os.access(path, os.X_OK)
