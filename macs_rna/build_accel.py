"""Build the macs3-rs accelerator binary and install it into the package.

Exposes ``macs-rna build-accel``: runs ``cargo build --release -p macs3-cli``
in ``rust/macs3-rs/`` and copies the resulting ``macs3-rs`` binary into
``macs_rna/_accel/macs3-rs`` so the installed package can use it. This is the
after-the-fact counterpart to the wheel build hook (``hatch_build.py``).
"""

import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger("macs_rna")


def _find_rust_dir() -> Optional[Path]:
    """Locate the ``rust/macs3-rs`` source directory.

    Returns
    -------
    pathlib.Path or None
        Path to the Rust workspace, or ``None`` if it cannot be found (e.g. in
        an installed wheel with no bundled sources).
    """
    # macs_rna/build_accel.py -> repo root -> rust/macs3-rs
    repo_root = Path(__file__).resolve().parent.parent
    rust_dir = repo_root / "rust" / "macs3-rs"
    return rust_dir if (rust_dir / "Cargo.toml").is_file() else None


def build_accel() -> int:
    """Build the macs3-rs binary and copy it into ``macs_rna/_accel/``.

    Returns
    -------
    int
        ``0`` on success, non-zero on failure (cargo missing, build failed, or
        sources unavailable).
    """
    if os.environ.get("MACS_RNA_NO_ACCEL"):
        logger.warning("MACS_RNA_NO_ACCEL set; skipping accelerator build")
        return 0

    if shutil.which("cargo") is None:
        logger.error("cargo not found on PATH; cannot build macs3-rs accelerator")
        return 1

    rust_dir = _find_rust_dir()
    if rust_dir is None:
        logger.error(
            "Could not find rust/macs3-rs sources; this only works from a "
            "source checkout, not an installed wheel"
        )
        return 1

    logger.info("Building macs3-rs accelerator in %s ...", rust_dir)
    try:
        subprocess.run(
            ["cargo", "build", "--release", "-p", "macs3-cli"],
            cwd=str(rust_dir),
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        logger.error("cargo build failed (exit %d)", exc.returncode)
        return exc.returncode or 1

    built = rust_dir / "target" / "release" / "macs3-rs"
    if not built.is_file():
        logger.error("Build succeeded but binary not found at %s", built)
        return 1

    accel_dir = Path(__file__).resolve().parent / "_accel"
    accel_dir.mkdir(parents=True, exist_ok=True)
    dest = accel_dir / "macs3-rs"
    shutil.copy2(built, dest)
    os.chmod(dest, 0o755)
    logger.info("Installed accelerator: %s", dest)
    return 0


def main() -> None:
    """Console-script entry point for ``macs-rna build-accel``."""
    from macs_rna.utils import setup_logging

    setup_logging(2)
    sys.exit(build_accel())
