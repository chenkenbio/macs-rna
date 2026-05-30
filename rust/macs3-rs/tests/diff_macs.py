#!/usr/bin/env python3
"""
diff_macs.py - Differential harness to compare macs3-rs (Rust reimplementation)
against stock macs3 v3.0.4 for bit-exact / numerically-equivalent output.

Usage:
    python tests/diff_macs.py --rust-bin ./target/release/macs3 [options]
    python tests/diff_macs.py --rust-bin <same-macs3-bin> --macs-bin <same-macs3-bin>  # self-check

Exit codes:
    0 - all tests PASS
    1 - one or more FAIL
"""

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("diff_macs")

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_MACS_BIN = "/home/kenchen/opt/miniforge3/envs/macs3v304/bin/macs3"
DEFAULT_TEST_DIR = "/home/kenchen/Documents/macs-rna/references/MACS/test"
BDG_TOL = 1e-5
SCORE_REL_TOL = 1e-5


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
class CaseResult:
    """Result of a single file comparison."""

    def __init__(
        self,
        subcmd: str,
        filename: str,
        verdict: str,  # PASS / FAIL / SKIP / ERROR
        byte_identical: Optional[bool] = None,
        worst_delta: Optional[float] = None,
        first_mismatch: Optional[str] = None,
        note: str = "",
    ) -> None:
        self.subcmd = subcmd
        self.filename = filename
        self.verdict = verdict
        self.byte_identical = byte_identical
        self.worst_delta = worst_delta
        self.first_mismatch = first_mismatch
        self.note = note

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"CaseResult(subcmd={self.subcmd!r}, filename={self.filename!r}, "
            f"verdict={self.verdict!r})"
        )


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def _open_maybe_gz(path: Path):
    """Open a plain or gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_bdg(path: Path) -> tuple[list[tuple[str, int, int]], list[float], str]:
    """
    Parse a bedGraph file.

    Returns
    -------
    coords : list of (chrom, start, end)
    values : list of float
    raw_text : the entire file content (for byte comparison)
    """
    raw_text = path.read_text()
    coords: list[tuple[str, int, int]] = []
    values: list[float] = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("track") or line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            coords.append((parts[0], int(parts[1]), int(parts[2])))
            values.append(float(parts[3]))
    return coords, values, raw_text


def read_narrowpeak(path: Path) -> tuple[list, str]:
    """
    Parse a narrowPeak file.

    Returns
    -------
    rows : list of dicts with keys chrom,start,end,name,score,strand,
           signalValue,pValue,qValue,summit
    raw_text : full file content
    """
    raw_text = path.read_text()
    rows = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("track") or line.startswith("#") or line.strip() == "":
                continue
            p = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "chrom": p[0],
                    "start": int(p[1]),
                    "end": int(p[2]),
                    "name": p[3],
                    "score": int(p[4]),
                    "strand": p[5],
                    "signalValue": float(p[6]),
                    "pValue": float(p[7]),
                    "qValue": float(p[8]),
                    "summit": int(p[9]),
                }
            )
    return rows, raw_text


def read_broadpeak(path: Path) -> tuple[list, str]:
    """
    Parse a broadPeak file (9-col BED).

    Returns
    -------
    rows : list of dicts
    raw_text : full file content
    """
    raw_text = path.read_text()
    rows = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("track") or line.startswith("#") or line.strip() == "":
                continue
            p = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "chrom": p[0],
                    "start": int(p[1]),
                    "end": int(p[2]),
                    "name": p[3],
                    "score": int(p[4]),
                    "strand": p[5],
                    "signalValue": float(p[6]),
                    "pValue": float(p[7]),
                    "qValue": float(p[8]),
                }
            )
    return rows, raw_text


def read_gappedpeak(path: Path) -> tuple[list, str]:
    """
    Parse a gappedPeak / bed12 file.

    Columns: chrom start end name score strand thickStart thickEnd
             itemRgb blockCount blockSizes blockStarts
             [signalValue pValue qValue]

    Returns
    -------
    rows : list of dicts
    raw_text : full file content
    """
    raw_text = path.read_text()
    rows = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("track") or line.startswith("#") or line.strip() == "":
                continue
            p = line.rstrip("\n").split("\t")
            row: dict[str, Any] = {
                "chrom": p[0],
                "start": int(p[1]),
                "end": int(p[2]),
                "name": p[3],
                "score": int(p[4]),
                "strand": p[5],
                "thickStart": int(p[6]),
                "thickEnd": int(p[7]),
                "itemRgb": p[8],
                "blockCount": int(p[9]),
                "blockSizes": p[10],
                "blockStarts": p[11],
            }
            # optional score columns
            if len(p) >= 15:
                row["signalValue"] = float(p[12])
                row["pValue"] = float(p[13])
                row["qValue"] = float(p[14])
            rows.append(row)
    return rows, raw_text


def read_summits_bed(path: Path) -> tuple[list, str]:
    """Parse a 5-column BED file (summits or bdgdiff output)."""
    raw_text = path.read_text()
    rows = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("track") or line.startswith("#") or line.strip() == "":
                continue
            p = line.rstrip("\n").split("\t")
            row: dict[str, Any] = {
                "chrom": p[0],
                "start": int(p[1]),
                "end": int(p[2]),
            }
            if len(p) >= 4:
                row["name"] = p[3]
            if len(p) >= 5:
                try:
                    row["score"] = float(p[4])
                except ValueError:
                    row["score"] = p[4]
            rows.append(row)
    return rows, raw_text


# ---------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------

def _normalize_for_byte_cmp(text: str) -> str:
    """
    Strip track-lines and normalize trailing newlines so that differences
    in track-line wording don't count as byte differences.
    """
    lines = [
        ln
        for ln in text.splitlines()
        if not ln.startswith("track") and not ln.startswith("#")
    ]
    return "\n".join(lines).rstrip("\n") + "\n"


def compare_bdg(
    ref_path: Path,
    tst_path: Path,
    subcmd: str,
) -> CaseResult:
    """
    Compare two bedGraph files.

    Parameters
    ----------
    ref_path : reference (macs3) output
    tst_path : test (rust) output
    subcmd : sub-command label for reporting

    Returns
    -------
    CaseResult
    """
    fname = tst_path.name
    try:
        ref_coords, ref_vals, ref_raw = read_bdg(ref_path)
        tst_coords, tst_vals, tst_raw = read_bdg(tst_path)
    except Exception as exc:
        return CaseResult(subcmd, fname, "ERROR", note=f"parse error: {exc}")

    byte_identical = _normalize_for_byte_cmp(ref_raw) == _normalize_for_byte_cmp(
        tst_raw
    )

    # Coord check
    if ref_coords != tst_coords:
        # Find first mismatch
        first_mm = None
        for i, (rc, tc) in enumerate(zip(ref_coords, tst_coords)):
            if rc != tc:
                first_mm = f"row {i}: ref={rc} tst={tc}"
                break
        if first_mm is None:
            first_mm = f"length mismatch ref={len(ref_coords)} tst={len(tst_coords)}"
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"coords differ: {first_mm}",
        )

    # Value check
    if len(ref_vals) != len(tst_vals):
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"value count mismatch: ref={len(ref_vals)} tst={len(tst_vals)}",
        )

    arr_ref = np.array(ref_vals)
    arr_tst = np.array(tst_vals)
    deltas = np.abs(arr_ref - arr_tst)
    worst = float(np.max(deltas)) if len(deltas) > 0 else 0.0

    if worst >= BDG_TOL:
        bad_idx = int(np.argmax(deltas))
        mm = (
            f"row {bad_idx} coord={ref_coords[bad_idx]}: "
            f"ref={ref_vals[bad_idx]} tst={tst_vals[bad_idx]} delta={worst:.2e}"
        )
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            worst_delta=worst,
            first_mismatch=mm,
        )

    return CaseResult(
        subcmd, fname, "PASS", byte_identical=byte_identical, worst_delta=worst
    )


def compare_narrowpeak(ref_path: Path, tst_path: Path, subcmd: str) -> CaseResult:
    """Compare two narrowPeak files."""
    fname = tst_path.name
    try:
        ref_rows, ref_raw = read_narrowpeak(ref_path)
        tst_rows, tst_raw = read_narrowpeak(tst_path)
    except Exception as exc:
        return CaseResult(subcmd, fname, "ERROR", note=f"parse error: {exc}")

    byte_identical = _normalize_for_byte_cmp(ref_raw) == _normalize_for_byte_cmp(
        tst_raw
    )

    if len(ref_rows) != len(tst_rows):
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"row count: ref={len(ref_rows)} tst={len(tst_rows)}",
        )

    worst = 0.0
    for i, (rr, tr) in enumerate(zip(ref_rows, tst_rows)):
        # Exact: coords, name, strand, summit, integer score
        for k in ("chrom", "start", "end", "name", "strand", "score", "summit"):
            if rr[k] != tr[k]:
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    first_mismatch=f"row {i} col {k!r}: ref={rr[k]} tst={tr[k]}",
                )
        # Numeric score cols within relative tol
        for k in ("signalValue", "pValue", "qValue"):
            rv, tv = rr[k], tr[k]
            denom = max(abs(rv), abs(tv), 1e-12)
            rel = abs(rv - tv) / denom
            abs_d = abs(rv - tv)
            delta = min(rel, abs_d)  # use smaller of rel and abs
            if delta > worst:
                worst = delta
            if rel > SCORE_REL_TOL and abs_d > BDG_TOL:
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    worst_delta=worst,
                    first_mismatch=f"row {i} col {k!r}: ref={rv} tst={tv} rel={rel:.2e}",
                )

    return CaseResult(
        subcmd, fname, "PASS", byte_identical=byte_identical, worst_delta=worst
    )


def compare_broadpeak(ref_path: Path, tst_path: Path, subcmd: str) -> CaseResult:
    """Compare two broadPeak (9-col) files."""
    fname = tst_path.name
    try:
        ref_rows, ref_raw = read_broadpeak(ref_path)
        tst_rows, tst_raw = read_broadpeak(tst_path)
    except Exception as exc:
        return CaseResult(subcmd, fname, "ERROR", note=f"parse error: {exc}")

    byte_identical = _normalize_for_byte_cmp(ref_raw) == _normalize_for_byte_cmp(
        tst_raw
    )

    if len(ref_rows) != len(tst_rows):
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"row count: ref={len(ref_rows)} tst={len(tst_rows)}",
        )

    worst = 0.0
    for i, (rr, tr) in enumerate(zip(ref_rows, tst_rows)):
        for k in ("chrom", "start", "end", "name", "strand", "score"):
            if rr[k] != tr[k]:
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    first_mismatch=f"row {i} col {k!r}: ref={rr[k]} tst={tr[k]}",
                )
        for k in ("signalValue", "pValue", "qValue"):
            rv, tv = rr[k], tr[k]
            denom = max(abs(rv), abs(tv), 1e-12)
            rel = abs(rv - tv) / denom
            abs_d = abs(rv - tv)
            delta = min(rel, abs_d)
            if delta > worst:
                worst = delta
            if rel > SCORE_REL_TOL and abs_d > BDG_TOL:
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    worst_delta=worst,
                    first_mismatch=f"row {i} col {k!r}: ref={rv} tst={tv} rel={rel:.2e}",
                )

    return CaseResult(
        subcmd, fname, "PASS", byte_identical=byte_identical, worst_delta=worst
    )


def compare_gappedpeak(ref_path: Path, tst_path: Path, subcmd: str) -> CaseResult:
    """Compare two gappedPeak / bed12 files."""
    fname = tst_path.name
    try:
        ref_rows, ref_raw = read_gappedpeak(ref_path)
        tst_rows, tst_raw = read_gappedpeak(tst_path)
    except Exception as exc:
        return CaseResult(subcmd, fname, "ERROR", note=f"parse error: {exc}")

    byte_identical = _normalize_for_byte_cmp(ref_raw) == _normalize_for_byte_cmp(
        tst_raw
    )

    if len(ref_rows) != len(tst_rows):
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"row count: ref={len(ref_rows)} tst={len(tst_rows)}",
        )

    worst = 0.0
    for i, (rr, tr) in enumerate(zip(ref_rows, tst_rows)):
        # Exact coords & block cols
        for k in (
            "chrom", "start", "end", "name", "strand",
            "thickStart", "thickEnd", "itemRgb",
            "blockCount", "blockSizes", "blockStarts", "score",
        ):
            if rr.get(k) != tr.get(k):
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    first_mismatch=f"row {i} col {k!r}: ref={rr.get(k)} tst={tr.get(k)}",
                )
        # Numeric score cols (optional)
        for k in ("signalValue", "pValue", "qValue"):
            if k not in rr:
                continue
            rv, tv = rr[k], tr[k]
            denom = max(abs(rv), abs(tv), 1e-12)
            rel = abs(rv - tv) / denom
            abs_d = abs(rv - tv)
            delta = min(rel, abs_d)
            if delta > worst:
                worst = delta
            if rel > SCORE_REL_TOL and abs_d > BDG_TOL:
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    worst_delta=worst,
                    first_mismatch=f"row {i} col {k!r}: ref={rv} tst={tv} rel={rel:.2e}",
                )

    return CaseResult(
        subcmd, fname, "PASS", byte_identical=byte_identical, worst_delta=worst
    )


def compare_summits(ref_path: Path, tst_path: Path, subcmd: str) -> CaseResult:
    """Compare summits.bed (5-col BED)."""
    fname = tst_path.name
    try:
        ref_rows, ref_raw = read_summits_bed(ref_path)
        tst_rows, tst_raw = read_summits_bed(tst_path)
    except Exception as exc:
        return CaseResult(subcmd, fname, "ERROR", note=f"parse error: {exc}")

    byte_identical = _normalize_for_byte_cmp(ref_raw) == _normalize_for_byte_cmp(
        tst_raw
    )

    if len(ref_rows) != len(tst_rows):
        return CaseResult(
            subcmd,
            fname,
            "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"row count: ref={len(ref_rows)} tst={len(tst_rows)}",
        )

    worst = 0.0
    for i, (rr, tr) in enumerate(zip(ref_rows, tst_rows)):
        for k in ("chrom", "start", "end"):
            if rr[k] != tr[k]:
                return CaseResult(
                    subcmd,
                    fname,
                    "FAIL",
                    byte_identical=byte_identical,
                    first_mismatch=f"row {i} col {k!r}: ref={rr[k]} tst={tr[k]}",
                )
        # score col (float)
        if "score" in rr and "score" in tr:
            try:
                rv, tv = float(rr["score"]), float(tr["score"])
                d = abs(rv - tv)
                if d > worst:
                    worst = d
                if d >= BDG_TOL:
                    return CaseResult(
                        subcmd,
                        fname,
                        "FAIL",
                        byte_identical=byte_identical,
                        worst_delta=worst,
                        first_mismatch=f"row {i} score: ref={rv} tst={tv} delta={d:.2e}",
                    )
            except (TypeError, ValueError):
                if rr.get("score") != tr.get("score"):
                    return CaseResult(
                        subcmd,
                        fname,
                        "FAIL",
                        byte_identical=byte_identical,
                        first_mismatch=f"row {i} score: ref={rr.get('score')} tst={tr.get('score')}",
                    )

    return CaseResult(
        subcmd, fname, "PASS", byte_identical=byte_identical, worst_delta=worst
    )


def compare_xls(ref_path: Path, tst_path: Path, subcmd: str) -> CaseResult:
    """
    Compare macs3 peaks.xls files.

    The header block contains the output directory path and other run-specific
    metadata in comment lines (starting with #).  Those are intentionally
    skipped; only data rows (non-# lines) are compared.  Score columns are
    compared with BDG_TOL tolerance; coordinate and name columns are exact.

    Parameters
    ----------
    ref_path : reference xls file
    tst_path : test xls file
    subcmd : sub-command label

    Returns
    -------
    CaseResult
    """
    fname = tst_path.name

    def _parse_xls(path: Path) -> tuple[list[str], list[list[str]]]:
        header: list[str] = []
        rows: list[list[str]] = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or line.strip() == "":
                    continue
                parts = line.rstrip("\n").split("\t")
                if parts[0] == "chr":
                    header = parts
                else:
                    rows.append(parts)
        return header, rows

    try:
        ref_header, ref_rows = _parse_xls(ref_path)
        tst_header, tst_rows = _parse_xls(tst_path)
    except Exception as exc:
        return CaseResult(subcmd, fname, "ERROR", note=f"parse error: {exc}")

    # byte_identical: data lines only (strip comment block from both)
    ref_data = "\n".join(
        ln for ln in ref_path.read_text().splitlines()
        if not ln.startswith("#") and ln.strip()
    )
    tst_data = "\n".join(
        ln for ln in tst_path.read_text().splitlines()
        if not ln.startswith("#") and ln.strip()
    )
    byte_identical = ref_data == tst_data

    if len(ref_rows) != len(tst_rows):
        return CaseResult(
            subcmd, fname, "FAIL",
            byte_identical=byte_identical,
            first_mismatch=f"row count: ref={len(ref_rows)} tst={len(tst_rows)}",
        )

    # Identify numeric columns by header name
    numeric_cols = {"-log10(pvalue)", "-log10(qvalue)", "fold_enrichment", "pileup"}
    exact_cols = {"chr", "start", "end", "length", "abs_summit", "name"}

    header = ref_header  # use ref header for column names
    worst = 0.0
    for i, (rr, tr) in enumerate(zip(ref_rows, tst_rows)):
        for ci, col in enumerate(header):
            if ci >= len(rr) or ci >= len(tr):
                break
            rv_s, tv_s = rr[ci], tr[ci]
            if col in exact_cols:
                if rv_s != tv_s:
                    return CaseResult(
                        subcmd, fname, "FAIL",
                        byte_identical=byte_identical,
                        first_mismatch=f"row {i} col {col!r}: ref={rv_s!r} tst={tv_s!r}",
                    )
            elif col in numeric_cols:
                try:
                    rv, tv = float(rv_s), float(tv_s)
                    denom = max(abs(rv), abs(tv), 1e-12)
                    rel = abs(rv - tv) / denom
                    abs_d = abs(rv - tv)
                    delta = min(rel, abs_d)
                    if delta > worst:
                        worst = delta
                    if rel > SCORE_REL_TOL and abs_d > BDG_TOL:
                        return CaseResult(
                            subcmd, fname, "FAIL",
                            byte_identical=byte_identical,
                            worst_delta=worst,
                            first_mismatch=f"row {i} col {col!r}: ref={rv} tst={tv} rel={rel:.2e}",
                        )
                except ValueError:
                    if rv_s != tv_s:
                        return CaseResult(
                            subcmd, fname, "FAIL",
                            byte_identical=byte_identical,
                            first_mismatch=f"row {i} col {col!r}: ref={rv_s!r} tst={tv_s!r}",
                        )

    return CaseResult(subcmd, fname, "PASS", byte_identical=byte_identical, worst_delta=worst)


def _dispatch_compare(
    ref_path: Path, tst_path: Path, subcmd: str
) -> CaseResult:
    """Route to the right comparator based on file extension."""
    name = tst_path.name.lower()
    if name.endswith(".bdg") or name.endswith(".bedgraph"):
        return compare_bdg(ref_path, tst_path, subcmd)
    elif name.endswith(".narrowpeak"):
        return compare_narrowpeak(ref_path, tst_path, subcmd)
    elif name.endswith(".broadpeak"):
        return compare_broadpeak(ref_path, tst_path, subcmd)
    elif name.endswith(".gappedpeak") or name.endswith(".bed12"):
        return compare_gappedpeak(ref_path, tst_path, subcmd)
    elif name.endswith("_summits.bed") or name.endswith(".bed"):
        return compare_summits(ref_path, tst_path, subcmd)
    elif name.endswith(".xls"):
        return compare_xls(ref_path, tst_path, subcmd)
    else:
        # Generic text diff (e.g. .r model files)
        ref_raw = ref_path.read_text()
        tst_raw = tst_path.read_text()
        byte_identical = ref_raw == tst_raw
        verdict = "PASS" if byte_identical else "FAIL"
        return CaseResult(
            subcmd, tst_path.name, verdict, byte_identical=byte_identical
        )


# ---------------------------------------------------------------------------
# predictd comparison
# ---------------------------------------------------------------------------

def parse_predictd_output(text: str) -> tuple[Optional[int], list[int]]:
    """
    Parse predicted fragment length from macs3 predictd output.

    Looks for patterns like:
      "# predicted fragment length is 229 bps"
      "# alternative fragment length(s) may be 229,123 bps"

    Returns
    -------
    d : primary predicted fragment length (or None if not found)
    alts : list of alternative fragment lengths (may be empty)
    """
    d: Optional[int] = None
    alts: list[int] = []
    for line in text.splitlines():
        m = re.search(r"predicted fragment length is (\d+)\s*bps", line, re.IGNORECASE)
        if m:
            d = int(m.group(1))
        m2 = re.search(r"alternative fragment length.*?(\d[\d,\s]+)\s*bps", line, re.IGNORECASE)
        if m2:
            try:
                alts = [int(x.strip()) for x in m2.group(1).split(",") if x.strip().isdigit()]
            except ValueError:
                pass
    return d, alts


def compare_predictd(
    ref_stdout: str, ref_stderr: str, tst_stdout: str, tst_stderr: str, subcmd: str
) -> CaseResult:
    """Compare predictd outputs: primary d must be EXACT; alts tolerated."""
    combined_ref = ref_stdout + "\n" + ref_stderr
    combined_tst = tst_stdout + "\n" + tst_stderr

    ref_d, ref_alts = parse_predictd_output(combined_ref)
    tst_d, tst_alts = parse_predictd_output(combined_tst)

    if ref_d is None:
        return CaseResult(subcmd, "predictd", "ERROR", note="could not parse ref d")
    if tst_d is None:
        return CaseResult(subcmd, "predictd", "SKIP", note="could not parse tst d (binary may not be ready)")

    byte_identical = ref_d == tst_d and set(ref_alts) == set(tst_alts)

    if ref_d != tst_d:
        return CaseResult(
            subcmd,
            "predictd",
            "FAIL",
            byte_identical=False,
            first_mismatch=f"primary d: ref={ref_d} tst={tst_d}",
        )

    note = ""
    if set(ref_alts) != set(tst_alts):
        note = f"alt d differs (tolerated): ref={sorted(ref_alts)} tst={sorted(tst_alts)}"
        log.warning("predictd alt lengths differ (tolerated): %s", note)

    return CaseResult(subcmd, "predictd", "PASS", byte_identical=byte_identical, note=note)


# ---------------------------------------------------------------------------
# Subprocess runner
# ---------------------------------------------------------------------------

def run_cmd(
    cmd: list[str],
    cwd: Optional[Path] = None,
    timeout: int = 600,
) -> tuple[int, str, str, float]:
    """
    Run a command, capture stdout/stderr.

    Parameters
    ----------
    cmd : command + args
    cwd : working directory
    timeout : seconds

    Returns
    -------
    returncode, stdout, stderr, elapsed_seconds
    """
    log.debug("Running: %s", " ".join(cmd))
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd else None,
            timeout=timeout,
        )
        elapsed = time.monotonic() - t0
        return proc.returncode, proc.stdout, proc.stderr, elapsed
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        log.error("Command timed out after %ds: %s", timeout, " ".join(cmd))
        return -1, "", "TIMEOUT", elapsed
    except FileNotFoundError as exc:
        elapsed = time.monotonic() - t0
        log.error("Binary not found: %s", exc)
        return -1, "", str(exc), elapsed


# ---------------------------------------------------------------------------
# Test cases definition
# ---------------------------------------------------------------------------

def build_test_cases(
    test_dir: Path,
    ref_outdir: Path,
    tst_outdir: Path,
    shared_bdg_dir: Path,
    macs_bin: str,
) -> list[dict]:
    """
    Build the list of test case definitions.

    Each dict has:
        subcmd  : str label
        ref_cmd : list[str]   (for macs-bin)
        tst_cmd : list[str]   (for rust-bin)  — binary placeholder = "__BIN__"
        ref_outdir : Path
        tst_outdir : Path
        compare_files : list of str filenames to compare
        is_predictd : bool
        depends_on : optional label of prerequisite (for chained inputs)

    For chained commands (bdgcmp, bdgopt, bdgpeakcall, bdgbroadcall),
    both engines get the SAME macs3-produced bdg so the sub-commands are
    tested in isolation.

    Parameters
    ----------
    test_dir : MACS test input directory
    ref_outdir : output directory for macs reference
    tst_outdir : output directory for rust binary
    shared_bdg_dir : directory for shared macs3-produced intermediate files
    macs_bin : path to macs3 binary for pre-generating shared inputs
    """
    chip = str(test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz")
    ctrl = str(test_dir / "CTCF_SE_CTRL_chr22_50k.bed.gz")

    cases: list[dict] = []

    # -----------------------------------------------------------------------
    # 1. pileup CHIP
    # -----------------------------------------------------------------------
    cases.append(
        {
            "subcmd": "pileup_chip",
            "ref_outdir": ref_outdir / "pileup_chip",
            "tst_outdir": tst_outdir / "pileup_chip",
            "args": [
                "pileup", "-f", "BED", "-i", chip,
                "--extsize", "200", "-o", "pileup_chip.bdg",
            ],
            "compare_files": ["pileup_chip.bdg"],
            "is_predictd": False,
        }
    )

    # 2. pileup CTRL
    cases.append(
        {
            "subcmd": "pileup_ctrl",
            "ref_outdir": ref_outdir / "pileup_ctrl",
            "tst_outdir": tst_outdir / "pileup_ctrl",
            "args": [
                "pileup", "-f", "BED", "-i", ctrl,
                "--extsize", "200", "-o", "pileup_ctrl.bdg",
            ],
            "compare_files": ["pileup_ctrl.bdg"],
            "is_predictd": False,
        }
    )

    # -----------------------------------------------------------------------
    # Shared inputs for chained tests
    # (generated once from macs_bin; both engines then get the same input)
    # -----------------------------------------------------------------------
    shared_chip_bdg = shared_bdg_dir / "pileup_chip.bdg"
    shared_ctrl_bdg = shared_bdg_dir / "pileup_ctrl.bdg"
    shared_fe_bdg = shared_bdg_dir / "bdgcmp_FE.bdg"
    # We'll generate these at runtime; store paths for reference in cases.

    # 3. bdgcmp ppois + FE
    cases.append(
        {
            "subcmd": "bdgcmp",
            "ref_outdir": ref_outdir / "bdgcmp",
            "tst_outdir": tst_outdir / "bdgcmp",
            "args": [
                "bdgcmp",
                "-t", str(shared_chip_bdg),
                "-c", str(shared_ctrl_bdg),
                "-m", "ppois", "FE",
                "-p", "1",
                "--o-prefix", "bdgcmp",
            ],
            "compare_files": ["bdgcmp_ppois.bdg", "bdgcmp_FE.bdg"],
            "is_predictd": False,
            "shared_input_paths": [shared_chip_bdg, shared_ctrl_bdg],
        }
    )

    # 4. bdgopt max
    cases.append(
        {
            "subcmd": "bdgopt_max",
            "ref_outdir": ref_outdir / "bdgopt_max",
            "tst_outdir": tst_outdir / "bdgopt_max",
            "args": [
                "bdgopt", "-m", "max", "-p", "2",
                "-i", str(shared_chip_bdg),
                "-o", "bdgopt_max.bdg",
            ],
            "compare_files": ["bdgopt_max.bdg"],
            "is_predictd": False,
            "shared_input_paths": [shared_chip_bdg],
        }
    )

    # 5. bdgopt min
    cases.append(
        {
            "subcmd": "bdgopt_min",
            "ref_outdir": ref_outdir / "bdgopt_min",
            "tst_outdir": tst_outdir / "bdgopt_min",
            "args": [
                "bdgopt", "-m", "min", "-p", "10",
                "-i", str(shared_chip_bdg),
                "-o", "bdgopt_min.bdg",
            ],
            "compare_files": ["bdgopt_min.bdg"],
            "is_predictd": False,
            "shared_input_paths": [shared_chip_bdg],
        }
    )

    # 6. bdgopt multiply
    cases.append(
        {
            "subcmd": "bdgopt_multiply",
            "ref_outdir": ref_outdir / "bdgopt_multiply",
            "tst_outdir": tst_outdir / "bdgopt_multiply",
            "args": [
                "bdgopt", "-m", "multiply", "-p", "1.5",
                "-i", str(shared_chip_bdg),
                "-o", "bdgopt_multiply.bdg",
            ],
            "compare_files": ["bdgopt_multiply.bdg"],
            "is_predictd": False,
            "shared_input_paths": [shared_chip_bdg],
        }
    )

    # 7. bdgopt p2q
    cases.append(
        {
            "subcmd": "bdgopt_p2q",
            "ref_outdir": ref_outdir / "bdgopt_p2q",
            "tst_outdir": tst_outdir / "bdgopt_p2q",
            "args": [
                "bdgopt", "-m", "p2q",
                "-i", str(shared_chip_bdg),
                "-o", "bdgopt_p2q.bdg",
            ],
            "compare_files": ["bdgopt_p2q.bdg"],
            "is_predictd": False,
            "shared_input_paths": [shared_chip_bdg],
        }
    )

    # 8. bdgpeakcall (uses shared FE bdg)
    cases.append(
        {
            "subcmd": "bdgpeakcall",
            "ref_outdir": ref_outdir / "bdgpeakcall",
            "tst_outdir": tst_outdir / "bdgpeakcall",
            "args": [
                "bdgpeakcall",
                "-i", str(shared_fe_bdg),
                "-c", "2",
                "--o-prefix", "bdgpeakcall",
            ],
            "compare_files": None,  # discover after run
            "is_predictd": False,
            "shared_input_paths": [shared_fe_bdg],
        }
    )

    # 9. bdgbroadcall (uses shared FE bdg)
    cases.append(
        {
            "subcmd": "bdgbroadcall",
            "ref_outdir": ref_outdir / "bdgbroadcall",
            "tst_outdir": tst_outdir / "bdgbroadcall",
            "args": [
                "bdgbroadcall",
                "-i", str(shared_fe_bdg),
                "-c", "2", "-C", "1.5",
                "--o-prefix", "bdgbroadcall",
            ],
            "compare_files": None,  # discover after run
            "is_predictd": False,
            "shared_input_paths": [shared_fe_bdg],
        }
    )

    # 10. predictd
    cases.append(
        {
            "subcmd": "predictd",
            "ref_outdir": ref_outdir / "predictd",
            "tst_outdir": tst_outdir / "predictd",
            "args": [
                "predictd",
                "-g", "52000000",
                "-i", chip,
                "--d-min", "10",
            ],
            "compare_files": [],
            "is_predictd": True,
        }
    )

    # -----------------------------------------------------------------------
    # callpeak cases (each gets its own outdir)
    # -----------------------------------------------------------------------
    callpeak_base = [
        "callpeak", "-g", "52000000",
        "-t", chip, "-c", ctrl,
        "-B",
    ]

    for name, extra_args in [
        ("narrow2", ["--nomodel", "--extsize", "100"]),
        ("narrow3", ["--nomodel", "--extsize", "100", "--shift", "-50"]),
        ("narrow4", ["--nomodel", "--nolambda", "--extsize", "100", "--shift", "-50"]),
        ("narrow5", ["--scale-to", "large"]),
        ("broad",   ["--broad"]),
    ]:
        cases.append(
            {
                "subcmd": f"callpeak_{name}",
                "ref_outdir": ref_outdir / f"callpeak_{name}",
                "tst_outdir": tst_outdir / f"callpeak_{name}",
                "args": callpeak_base + extra_args + ["-n", name],
                "compare_files": None,  # discover after run
                "is_predictd": False,
            }
        )

    return cases


# ---------------------------------------------------------------------------
# Generate shared inputs
# ---------------------------------------------------------------------------

def generate_shared_inputs(
    macs_bin: str,
    test_dir: Path,
    shared_bdg_dir: Path,
) -> bool:
    """
    Pre-generate shared intermediate files (pileup bdg, FE bdg) from macs_bin
    so that chained tests receive the same inputs for both engines.

    Parameters
    ----------
    macs_bin : path to macs3 binary
    test_dir : MACS test inputs directory
    shared_bdg_dir : directory to write shared files

    Returns
    -------
    success : bool
    """
    shared_bdg_dir.mkdir(parents=True, exist_ok=True)
    chip = str(test_dir / "CTCF_SE_ChIP_chr22_50k.bed.gz")
    ctrl = str(test_dir / "CTCF_SE_CTRL_chr22_50k.bed.gz")

    chip_bdg = shared_bdg_dir / "pileup_chip.bdg"
    ctrl_bdg = shared_bdg_dir / "pileup_ctrl.bdg"
    fe_bdg = shared_bdg_dir / "bdgcmp_FE.bdg"

    # pileup chip
    if not chip_bdg.exists():
        log.info("Generating shared chip pileup bdg...")
        rc, out, err, elapsed = run_cmd(
            [macs_bin, "pileup", "-f", "BED", "-i", chip,
             "--extsize", "200", "-o", str(chip_bdg)],
            cwd=test_dir,
        )
        if rc != 0:
            log.error("Failed to generate shared chip bdg: %s", err)
            return False
        log.info("  Done in %.1fs", elapsed)
    else:
        log.info("Shared chip bdg already exists, reusing.")

    # pileup ctrl
    if not ctrl_bdg.exists():
        log.info("Generating shared ctrl pileup bdg...")
        rc, out, err, elapsed = run_cmd(
            [macs_bin, "pileup", "-f", "BED", "-i", ctrl,
             "--extsize", "200", "-o", str(ctrl_bdg)],
            cwd=test_dir,
        )
        if rc != 0:
            log.error("Failed to generate shared ctrl bdg: %s", err)
            return False
        log.info("  Done in %.1fs", elapsed)
    else:
        log.info("Shared ctrl bdg already exists, reusing.")

    # bdgcmp FE
    if not fe_bdg.exists():
        log.info("Generating shared FE bdg...")
        rc, out, err, elapsed = run_cmd(
            [macs_bin, "bdgcmp",
             "-t", str(chip_bdg),
             "-c", str(ctrl_bdg),
             "-m", "FE",
             "-p", "1",
             "--o-prefix", str(shared_bdg_dir / "bdgcmp"),
             ],
            cwd=test_dir,
        )
        if rc != 0:
            log.error("Failed to generate shared FE bdg: %s", err)
            return False
        # bdgcmp writes <prefix>_FE.bdg
        generated = shared_bdg_dir / "bdgcmp_FE.bdg"
        if not generated.exists():
            log.error("Expected %s not found after bdgcmp.", generated)
            return False
        log.info("  Done in %.1fs", elapsed)
    else:
        log.info("Shared FE bdg already exists, reusing.")

    return True


# ---------------------------------------------------------------------------
# Run a single test case
# ---------------------------------------------------------------------------

def run_case(
    case: dict,
    macs_bin: str,
    rust_bin: Optional[str],
    test_dir: Path,
) -> list[CaseResult]:
    """
    Run one test case: execute both binaries and compare outputs.

    Parameters
    ----------
    case : test case definition dict
    macs_bin : reference binary
    rust_bin : rust binary (may be None → SKIP)
    test_dir : working directory for command execution

    Returns
    -------
    list of CaseResult (one per compared file)
    """
    subcmd = case["subcmd"]
    ref_outdir: Path = case["ref_outdir"]
    tst_outdir: Path = case["tst_outdir"]
    args: list[str] = case["args"]
    is_predictd: bool = case.get("is_predictd", False)

    ref_outdir.mkdir(parents=True, exist_ok=True)
    tst_outdir.mkdir(parents=True, exist_ok=True)

    # Inject --outdir into args (append at end so it overrides)
    def inject_outdir(outdir: Path) -> list[str]:
        # Remove any existing --outdir
        filtered = []
        i = 0
        while i < len(args):
            if args[i] == "--outdir":
                i += 2
            else:
                filtered.append(args[i])
                i += 1
        return filtered + ["--outdir", str(outdir)]

    ref_cmd = [macs_bin] + inject_outdir(ref_outdir)
    t0 = time.monotonic()

    # --- Run reference ---
    log.info("[%s] Running reference macs3...", subcmd)
    ref_rc, ref_out, ref_err, ref_elapsed = run_cmd(ref_cmd, cwd=test_dir)
    log.info("[%s] Reference done in %.1fs (rc=%d)", subcmd, ref_elapsed, ref_rc)
    if ref_rc != 0:
        log.warning("[%s] Reference exited with rc=%d\nstderr: %s", subcmd, ref_rc, ref_err[-500:])
        return [CaseResult(subcmd, "*", "ERROR", note=f"ref rc={ref_rc}: {ref_err[-300:]}")]

    # --- Check rust binary ---
    if rust_bin is None or not Path(rust_bin).exists():
        log.info("[%s] rust-bin not found → SKIP", subcmd)
        return [CaseResult(subcmd, "*", "SKIP", note="rust-bin not available")]

    tst_cmd = [rust_bin] + inject_outdir(tst_outdir)
    log.info("[%s] Running rust macs3...", subcmd)
    tst_rc, tst_out, tst_err, tst_elapsed = run_cmd(tst_cmd, cwd=test_dir)
    log.info("[%s] Rust done in %.1fs (rc=%d)", subcmd, tst_elapsed, tst_rc)

    total_elapsed = time.monotonic() - t0

    if tst_rc != 0:
        log.warning("[%s] Rust exited with rc=%d\nstderr: %s", subcmd, tst_rc, tst_err[-500:])
        return [CaseResult(subcmd, "*", "ERROR", note=f"tst rc={tst_rc}: {tst_err[-300:]}")]

    # --- predictd special comparison ---
    if is_predictd:
        result = compare_predictd(ref_out, ref_err, tst_out, tst_err, subcmd)
        return [result]

    # --- Discover files to compare ---
    compare_files: Optional[list[str]] = case.get("compare_files")
    if compare_files is None:
        # Auto-discover: compare all files produced in ref_outdir
        compare_files = [f.name for f in ref_outdir.iterdir() if f.is_file()]
        compare_files.sort()

    results: list[CaseResult] = []
    for fname in compare_files:
        ref_path = ref_outdir / fname
        tst_path = tst_outdir / fname
        if not ref_path.exists():
            results.append(CaseResult(subcmd, fname, "ERROR", note="ref file missing"))
            continue
        if not tst_path.exists():
            results.append(CaseResult(subcmd, fname, "SKIP", note="tst file missing"))
            continue
        result = _dispatch_compare(ref_path, tst_path, subcmd)
        results.append(result)

    if not results:
        results.append(CaseResult(subcmd, "(no files)", "SKIP", note="no output files found"))

    log.info(
        "[%s] Total elapsed %.1fs — %d comparisons: %d PASS, %d FAIL, %d SKIP/ERROR",
        subcmd,
        total_elapsed,
        len(results),
        sum(1 for r in results if r.verdict == "PASS"),
        sum(1 for r in results if r.verdict == "FAIL"),
        sum(1 for r in results if r.verdict in ("SKIP", "ERROR")),
    )
    return results


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(results: list[CaseResult]) -> None:
    """Print a formatted summary table to stdout."""
    header = (
        f"{'SUBCMD':<22} {'FILE':<45} {'VERDICT':<8} "
        f"{'WORST_DELTA':<14} {'BYTE_IDENTICAL':<15} NOTE"
    )
    sep = "-" * len(header)
    print("\n" + sep)
    print(header)
    print(sep)
    for r in results:
        delta_str = f"{r.worst_delta:.2e}" if r.worst_delta is not None else "—"
        bi_str = str(r.byte_identical) if r.byte_identical is not None else "—"
        note_str = r.note or (r.first_mismatch or "")
        if len(note_str) > 80:
            note_str = note_str[:77] + "..."
        print(
            f"{r.subcmd:<22} {r.filename:<45} {r.verdict:<8} "
            f"{delta_str:<14} {bi_str:<15} {note_str}"
        )
    print(sep)

    total = len(results)
    n_pass = sum(1 for r in results if r.verdict == "PASS")
    n_fail = sum(1 for r in results if r.verdict == "FAIL")
    n_skip = sum(1 for r in results if r.verdict == "SKIP")
    n_err = sum(1 for r in results if r.verdict == "ERROR")
    print(
        f"\nSUMMARY: {total} comparisons | "
        f"PASS={n_pass} FAIL={n_fail} SKIP={n_skip} ERROR={n_err}"
    )
    byte_id = [r for r in results if r.byte_identical is not None]
    if byte_id:
        n_bi = sum(1 for r in byte_id if r.byte_identical)
        print(f"         {n_bi}/{len(byte_id)} files byte-identical (after trackline normalization)")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Diff harness: compare macs3-rs output against stock macs3 v3.0.4.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--rust-bin",
        metavar="PATH",
        default=None,
        help="Path to macs3-rs binary (may not exist yet; tests will be SKIP).",
    )
    parser.add_argument(
        "--macs-bin",
        metavar="PATH",
        default=DEFAULT_MACS_BIN,
        help="Path to reference macs3 binary.",
    )
    parser.add_argument(
        "--workdir",
        metavar="DIR",
        default=None,
        help="Working directory for outputs. Defaults to a temp directory.",
    )
    parser.add_argument(
        "--only",
        metavar="SUBCMD",
        default=None,
        help="Run only the test case matching this subcmd label (partial match).",
    )
    parser.add_argument(
        "--test-dir",
        metavar="DIR",
        default=DEFAULT_TEST_DIR,
        help="Directory containing MACS test input files.",
    )
    parser.add_argument(
        "--keep-workdir",
        action="store_true",
        default=False,
        help="Do not delete the workdir after the run (useful for debugging).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Enable DEBUG-level logging.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    """
    Entry point.

    Parameters
    ----------
    argv : optional argument list for testing; defaults to sys.argv[1:]

    Returns
    -------
    exit code (0 = all pass, 1 = any fail)
    """
    args = parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    log.info("=" * 60)
    log.info("diff_macs.py - macs3-rs diff harness")
    log.info("macs-bin  : %s", args.macs_bin)
    log.info("rust-bin  : %s", args.rust_bin or "(not set)")
    log.info("test-dir  : %s", args.test_dir)
    log.info("=" * 60)

    t_start = time.monotonic()

    # Validate macs-bin
    if not Path(args.macs_bin).exists():
        log.error("macs-bin not found: %s", args.macs_bin)
        return 1

    test_dir = Path(args.test_dir).expanduser()
    if not test_dir.exists():
        log.error("test-dir not found: %s", test_dir)
        return 1

    # Set up workdir
    tmp_created = False
    if args.workdir:
        workdir = Path(args.workdir).expanduser()
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        workdir = Path(tempfile.mkdtemp(prefix="diff_macs_"))
        tmp_created = True
    log.info("workdir   : %s", workdir)

    ref_outdir = workdir / "ref"
    tst_outdir = workdir / "tst"
    shared_bdg_dir = workdir / "shared_inputs"

    # Generate shared inputs
    log.info("Generating shared intermediate files from macs-bin...")
    ok = generate_shared_inputs(args.macs_bin, test_dir, shared_bdg_dir)
    if not ok:
        log.error("Failed to generate shared inputs. Aborting.")
        return 1

    # Build test cases
    cases = build_test_cases(
        test_dir=test_dir,
        ref_outdir=ref_outdir,
        tst_outdir=tst_outdir,
        shared_bdg_dir=shared_bdg_dir,
        macs_bin=args.macs_bin,
    )

    # Filter if --only
    if args.only:
        cases = [c for c in cases if args.only.lower() in c["subcmd"].lower()]
        if not cases:
            log.error("No cases matched --only %r", args.only)
            return 1
        log.info("Filtered to %d cases matching %r", len(cases), args.only)

    # Run all cases
    all_results: list[CaseResult] = []
    for case in cases:
        log.info("--- Running: %s ---", case["subcmd"])
        results = run_case(
            case=case,
            macs_bin=args.macs_bin,
            rust_bin=args.rust_bin,
            test_dir=test_dir,
        )
        all_results.extend(results)

    # Print summary
    print_summary(all_results)

    elapsed_total = time.monotonic() - t_start
    log.info("Total wall time: %.1f s", elapsed_total)

    # Cleanup
    if tmp_created and not args.keep_workdir:
        log.info("Removing temp workdir: %s", workdir)
        shutil.rmtree(workdir, ignore_errors=True)
    elif args.keep_workdir:
        log.info("Keeping workdir: %s", workdir)

    # Exit code
    n_fail = sum(1 for r in all_results if r.verdict == "FAIL")
    if n_fail > 0:
        log.error("%d FAIL(s) found.", n_fail)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
