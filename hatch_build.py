"""Hatchling build hook: bundle the macs3-rs accelerator into the wheel.

At wheel-build time this runs ``cargo build --release -p macs3-cli`` in
``rust/macs3-rs/`` and force-includes the resulting ``macs3-rs`` binary at
``macs_rna/_accel/macs3-rs`` inside the wheel. The build is best-effort: if
``cargo`` is missing, the sources are unavailable, or the build fails, a clear
warning is printed and the binary is simply omitted -- the wheel still builds
and works through the stock ``macs3`` fallback. Setting ``$MACS_RNA_NO_ACCEL``
skips the build entirely.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class AccelBuildHook(BuildHookInterface):
    """Build and force-include the macs3-rs accelerator binary."""

    PLUGIN_NAME = "custom"

    def initialize(self, version: str, build_data: dict) -> None:
        """Build macs3-rs (best-effort) and add it to the wheel.

        Parameters
        ----------
        version : str
            Build version string (unused).
        build_data : dict
            Hatchling build-data dict; the binary is registered via
            ``build_data["force_include"]``.
        """
        if os.environ.get("MACS_RNA_NO_ACCEL"):
            self._warn("MACS_RNA_NO_ACCEL set; skipping macs3-rs build")
            return

        root = Path(self.root)
        rust_dir = root / "rust" / "macs3-rs"
        if not (rust_dir / "Cargo.toml").is_file():
            self._warn(
                f"rust/macs3-rs sources not found at {rust_dir}; "
                "wheel will rely on the stock macs3 fallback"
            )
            return

        if shutil.which("cargo") is None:
            self._warn(
                "cargo not found on PATH; skipping macs3-rs build "
                "(wheel will rely on the stock macs3 fallback)"
            )
            return

        try:
            subprocess.run(
                ["cargo", "build", "--release", "-p", "macs3-cli"],
                cwd=str(rust_dir),
                check=True,
            )
        except subprocess.CalledProcessError as exc:
            self._warn(
                f"cargo build failed (exit {exc.returncode}); skipping "
                "macs3-rs accelerator (wheel will rely on the stock macs3 "
                "fallback)"
            )
            return

        built = rust_dir / "target" / "release" / "macs3-rs"
        if not built.is_file():
            self._warn(
                f"cargo build succeeded but binary missing at {built}; "
                "skipping accelerator"
            )
            return

        # Force-include the built binary at macs_rna/_accel/macs3-rs.
        build_data.setdefault("force_include", {})
        build_data["force_include"][str(built)] = "macs_rna/_accel/macs3-rs"
        # Mark the wheel as platform-specific (it ships a native binary).
        build_data["infer_tag"] = True
        build_data["pure_python"] = False
        self._info(f"Bundling macs3-rs accelerator from {built}")

    def _warn(self, msg: str) -> None:
        """Print a build-time warning to stderr."""
        print(f"[macs-rna build] WARNING: {msg}", file=sys.stderr)

    def _info(self, msg: str) -> None:
        """Print a build-time info message to stderr."""
        print(f"[macs-rna build] {msg}", file=sys.stderr)
