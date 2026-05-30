//! bedGraph reading / writing, ported from `MACS3/IO/BedGraphIO.py` and
//! `MACS3/Signal/cPosValCalculation.c::write_pv_array_to_bedGraph`.
//!
//! This module is **frozen after Phase 0** (read-only for later phases). It
//! covers two distinct on-disk writers that MACS uses, which differ in how they
//! decide where one bedGraph interval ends and the next begins:
//!
//!   * [`write_pv_exact`] — the C writer
//!     (`cPosValCalculation.c::write_pv_array_to_bedGraph`). It coalesces
//!     adjacent `PosVal` entries only when their `value` is **exactly** equal
//!     (`!=` on `float`), starts at `pre_s = 0`, and emits the leading interval
//!     `[0, pv[0].pos)`.
//!   * [`write_bedgraph_predicate`] — the Python writer family
//!     (`BedGraphIO.py::write_bedGraph`, `ScoreTrack.py::write_bedGraph`,
//!     `TwoConditionScores.write_bedGraph`). It walks `(pos, value)` arrays and
//!     breaks on a caller-supplied predicate over consecutive values:
//!       - BedGraphIO: break on every interval (no coalescing).
//!       - ScoreTrackII / CallPeakUnit: break when `abs(pre_v - v) > 1e-5`.
//!       - TwoConditionScores: break when `abs(pre_v - v) >= 1e-6`.
//!
//! All values are formatted with `%.5f` via [`crate::fmt::write_f5`].

use std::io::{self, Read, Write};

use crate::bedgraph::BedGraphTrack;
use crate::fmt::write_f5;
use crate::pileup::PosVal;

/// gzip magic bytes (`1f 8b`).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Parse a leading base-10 integer the way C `atoi` does: skip leading ASCII
/// whitespace, take an optional sign, then consume digits until the first
/// non-digit (returning `0` if no digits are present). MACS calls `atoi` on the
/// bedGraph start/end columns.
fn atoi(s: &[u8]) -> i32 {
    let mut i = 0;
    while i < s.len() && (s[i] as char).is_ascii_whitespace() {
        i += 1;
    }
    let mut sign: i64 = 1;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        if s[i] == b'-' {
            sign = -1;
        }
        i += 1;
    }
    let mut acc: i64 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        acc = acc * 10 + (s[i] - b'0') as i64;
        i += 1;
    }
    (sign * acc) as i32
}

/// Parse the value column to `f64` then narrow to `f32`, matching
/// `BedGraphIO.read_bedGraph`'s `atof(fs[3])` followed by storage in a C
/// `float` track. We parse the whole token as a Rust `f64` (which performs the
/// same correctly-rounded decimal-to-binary conversion as glibc `strtod`/`atof`
/// for valid numbers) and then apply the single `f64 -> f32` rounding that the
/// Cython `add_loc(..., value: cython.float)` coercion performs.
fn atof_f32(s: &[u8]) -> f32 {
    let txt = std::str::from_utf8(s).unwrap_or("");
    let trimmed = txt.trim();
    trimmed.parse::<f64>().unwrap_or(0.0) as f32
}

/// Open `path` and return a reader that transparently decompresses gzip input
/// (detected by the leading `1f 8b` magic bytes), matching MACS's reliance on a
/// gzip-aware open. Plain files are returned as-is.
fn open_maybe_gzip(path: &str) -> io::Result<Box<dyn Read>> {
    let mut f = std::fs::File::open(path)?;
    let mut magic = [0u8; 2];
    let n = read_up_to(&mut f, &mut magic)?;
    // Re-open so the consumed magic bytes are seen again by the chosen reader.
    let f2 = std::fs::File::open(path)?;
    if n == 2 && magic == GZIP_MAGIC {
        Ok(Box::new(flate2::read::MultiGzDecoder::new(f2)))
    } else {
        Ok(Box::new(f2))
    }
}

/// Read into `buf`, retrying short reads until it is full or EOF. Returns the
/// number of bytes read.
fn read_up_to(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buf.len() {
        match r.read(&mut buf[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(filled)
}

/// Read a bedGraph file into a [`BedGraphTrack`], mirroring
/// `BedGraphIO.read_bedGraph` (with `baseline_value = 0`).
///
/// Lines beginning with `track`, `#`, or `browse` are skipped. Remaining lines
/// are whitespace-split into `chrom start end value`; `start`/`end` are parsed
/// with [`atoi`] and `value` with [`atof_f32`], then fed to
/// [`BedGraphTrack::add_loc`] (which merges equal-valued neighbours, exactly as
/// the Cython `add_loc`).
pub fn read_bedgraph(path: &str) -> io::Result<BedGraphTrack> {
    let reader = open_maybe_gzip(path)?;
    let mut buf = io::BufReader::new(reader);
    let mut track = BedGraphTrack::new(0.0);

    let mut line: Vec<u8> = Vec::with_capacity(128);
    loop {
        line.clear();
        let n = read_line(&mut buf, &mut line)?;
        if n == 0 {
            break;
        }
        // Strip trailing newline / carriage return for clean tokenisation.
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.starts_with(b"track") || line.starts_with(b"#") || line.starts_with(b"browse") {
            continue;
        }
        if line.is_empty() {
            continue;
        }
        // Python `bytes.split()` splits on runs of ASCII whitespace.
        let mut fields = line
            .split(|&c| (c as char).is_ascii_whitespace())
            .filter(|f| !f.is_empty());
        let chrom = match fields.next() {
            Some(c) => c,
            None => continue,
        };
        let start = match fields.next() {
            Some(s) => atoi(s),
            None => continue,
        };
        let end = match fields.next() {
            Some(s) => atoi(s),
            None => continue,
        };
        let value = match fields.next() {
            Some(s) => atof_f32(s),
            None => continue,
        };
        track.add_loc(chrom, start, end, value);
    }
    Ok(track)
}

/// Read a single line (including the terminating `\n`) into `out`. Returns the
/// number of bytes read (0 at EOF).
fn read_line(r: &mut impl io::BufRead, out: &mut Vec<u8>) -> io::Result<usize> {
    r.read_until(b'\n', out)
}

/// Write a coalesced bedGraph for one chromosome, replicating the C writer
/// `cPosValCalculation.c::write_pv_array_to_bedGraph`.
///
/// `pv` is a `PosVal` array in bedGraph order (each entry's `pos` is the
/// right-open end of its interval). Intervals are coalesced only when their
/// `value` is **exactly** equal (`!=` comparison on `f32`). The first interval
/// spans `[0, pv[0].pos)`; subsequent breaks set `pre_s = pre_e`. Values are
/// emitted as `%.5f`. When `append` is true the bytes extend an existing file
/// region (the caller controls the underlying handle's position); this matches
/// the C `fopen("a")` vs `fopen("w")` choice at the call site.
///
/// Does nothing for an empty `pv` (the C code dereferences `pv_array[0]`
/// assuming a non-empty array, so MACS never calls it with an empty one).
pub fn write_pv_exact<W: Write>(
    chrom: &[u8],
    pv: &[PosVal],
    out: &mut W,
    _append: bool,
) -> io::Result<()> {
    if pv.is_empty() {
        return Ok(());
    }
    let mut pre_s: i32 = 0;
    let mut pre_e: i32 = pv[0].pos;
    let mut pre_v: f32 = pv[0].val;

    let mut buf = String::with_capacity(64);
    for entry in &pv[1..] {
        if entry.val != pre_v {
            write_bdg_line(out, &mut buf, chrom, pre_s, pre_e, pre_v)?;
            pre_s = pre_e;
            pre_e = entry.pos;
            pre_v = entry.val;
        } else {
            pre_e = entry.pos;
        }
    }
    // last piece
    write_bdg_line(out, &mut buf, chrom, pre_s, pre_e, pre_v)?;
    Ok(())
}

/// Format and write one `chrom\tstart\tend\t%.5f\n` line, reusing `buf`.
fn write_bdg_line<W: Write>(
    out: &mut W,
    buf: &mut String,
    chrom: &[u8],
    start: i32,
    end: i32,
    value: f32,
) -> io::Result<()> {
    buf.clear();
    write_f5(buf, value);
    out.write_all(chrom)?;
    write!(out, "\t{start}\t{end}\t")?;
    out.write_all(buf.as_bytes())?;
    out.write_all(b"\n")?;
    Ok(())
}

/// Break-predicate over two consecutive interval values, deciding whether a new
/// bedGraph line should be emitted. Implementations correspond to the three
/// Python writers (see module docs).
pub trait BreakPredicate {
    /// Return `true` when `cur` is "different enough" from `pre` that a line
    /// boundary must be written before `cur`.
    fn breaks(&self, pre: f32, cur: f32) -> bool;
}

/// Break on every interval — no coalescing (`BedGraphIO.write_bedGraph`).
pub struct AlwaysBreak;
impl BreakPredicate for AlwaysBreak {
    fn breaks(&self, _pre: f32, _cur: f32) -> bool {
        true
    }
}

/// Break when `abs(pre - cur) > 1e-5` (`ScoreTrackII` / `CallPeakUnit`).
pub struct Gt1e5;
impl BreakPredicate for Gt1e5 {
    fn breaks(&self, pre: f32, cur: f32) -> bool {
        (pre - cur).abs() > 1e-5
    }
}

/// Break when `abs(pre - cur) >= 1e-6` (`TwoConditionScores.write_bedGraph`).
pub struct Ge1e6;
impl BreakPredicate for Ge1e6 {
    fn breaks(&self, pre: f32, cur: f32) -> bool {
        (pre - cur).abs() >= 1e-6
    }
}

/// Generic Python-style bedGraph writer for one chromosome, parameterised by a
/// break predicate. Mirrors the common loop in `BedGraphIO.write_bedGraph` /
/// `ScoreTrackII.write_bedGraph` / `TwoConditionScores.write_bedGraph`:
///
/// ```text
/// pre = 0; pre_v = value[0]
/// for i in 1..ln:
///     v = value[i]; p = pos[i-1]
///     if predicate(pre_v, v):
///         write(chrom, pre, p, pre_v); pre_v = v; pre = p
/// write(chrom, pre, pos[ln-1], pre_v)   # last one
/// ```
///
/// `pos` and `value` must have equal length; an empty slice writes nothing.
/// For [`AlwaysBreak`] this reproduces the no-coalesce BedGraphIO output (every
/// interval emitted with `pre = pos[i-1]`).
pub fn write_bedgraph_predicate<W: Write, P: BreakPredicate>(
    chrom: &[u8],
    pos: &[i32],
    value: &[f32],
    pred: &P,
    out: &mut W,
) -> io::Result<()> {
    debug_assert_eq!(pos.len(), value.len());
    let ln = pos.len();
    if ln == 0 {
        return Ok(());
    }
    let mut buf = String::with_capacity(64);
    let mut pre: i32 = 0;
    let mut pre_v: f32 = value[0];
    for i in 1..ln {
        let v = value[i];
        let p = pos[i - 1];
        if pred.breaks(pre_v, v) {
            write_bdg_line(out, &mut buf, chrom, pre, p, pre_v)?;
            pre_v = v;
            pre = p;
        }
    }
    // last one
    write_bdg_line(out, &mut buf, chrom, pre, pos[ln - 1], pre_v)?;
    Ok(())
}

/// Optional UCSC bedGraph track line, as written by
/// `BedGraphIO.write_bedGraph` when `trackline` is set. Quotes in `name` /
/// `description` are backslash-escaped, matching the Python `.replace`.
pub fn bedgraph_trackline(name: &str, description: &str) -> String {
    let n = name.replace('"', "\\\"");
    let d = description.replace('"', "\\\"");
    format!("track type=bedGraph name=\"{n}\" description=\"{d}\" visibility=2 alwaysZero=on\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoi_matches_c() {
        assert_eq!(atoi(b"123"), 123);
        assert_eq!(atoi(b"  -45xyz"), -45);
        assert_eq!(atoi(b"+7"), 7);
        assert_eq!(atoi(b"abc"), 0);
        assert_eq!(atoi(b"0"), 0);
    }

    #[test]
    fn atof_f32_narrows_via_f64() {
        assert_eq!(atof_f32(b"1.5"), 1.5_f32);
        assert_eq!(atof_f32(b"0.10000"), 0.1_f32);
        // round-trip through f64 then narrow: parse 0.1 as f64 then `as f32`
        let want = 0.1_f64 as f32;
        assert_eq!(atof_f32(b"0.1"), want);
    }

    /// The golden case from `cPosValCalculation.c::main()`: the `max_pileup`
    /// array is
    ///   [(10,0.0),(11,0.5),(12,1.5),(13,2.0),(14,2.0),(15,1.5)]
    /// and the reference binary writes (with (13,2.0) coalesced into (14,2.0)):
    ///   chr1  0  10 0.00000
    ///   chr1 10  11 0.50000
    ///   chr1 11  12 1.50000
    ///   chr1 12  14 2.00000
    ///   chr1 14  15 1.50000
    #[test]
    fn write_pv_exact_matches_c_main_example() {
        let pv = vec![
            PosVal::new(10, 0.0),
            PosVal::new(11, 0.5),
            PosVal::new(12, 1.5),
            PosVal::new(13, 2.0),
            PosVal::new(14, 2.0),
            PosVal::new(15, 1.5),
        ];
        let mut out: Vec<u8> = Vec::new();
        write_pv_exact(b"chr1", &pv, &mut out, false).unwrap();
        let want = "chr1\t0\t10\t0.00000\n\
                    chr1\t10\t11\t0.50000\n\
                    chr1\t11\t12\t1.50000\n\
                    chr1\t12\t14\t2.00000\n\
                    chr1\t14\t15\t1.50000\n";
        assert_eq!(String::from_utf8(out).unwrap(), want);
    }

    #[test]
    fn write_pv_exact_single_interval() {
        let pv = vec![PosVal::new(100, 3.0)];
        let mut out: Vec<u8> = Vec::new();
        write_pv_exact(b"chrX", &pv, &mut out, false).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "chrX\t0\t100\t3.00000\n");
    }

    #[test]
    fn predicate_alwaysbreak_emits_every_interval() {
        // No coalescing: BedGraphIO-style output.
        let pos = vec![10, 20, 30];
        let val = vec![1.0_f32, 1.0, 2.0];
        let mut out: Vec<u8> = Vec::new();
        write_bedgraph_predicate(b"chr1", &pos, &val, &AlwaysBreak, &mut out).unwrap();
        // pre starts 0; i=1 p=pos[0]=10 break -> 0..10 v=1; i=2 p=pos[1]=20 break -> 10..20 v=1
        // last -> 20..30 v=2
        let want = "chr1\t0\t10\t1.00000\n\
                    chr1\t10\t20\t1.00000\n\
                    chr1\t20\t30\t2.00000\n";
        assert_eq!(String::from_utf8(out).unwrap(), want);
    }

    #[test]
    fn predicate_gt1e5_coalesces_near_equal() {
        let pos = vec![10, 20, 30];
        let val = vec![1.0_f32, 1.0, 2.0];
        let mut out: Vec<u8> = Vec::new();
        write_bedgraph_predicate(b"chr1", &pos, &val, &Gt1e5, &mut out).unwrap();
        // 1.0 and 1.0 don't break -> merged 0..20; then last 20..30 v=2
        let want = "chr1\t0\t20\t1.00000\n\
                    chr1\t20\t30\t2.00000\n";
        assert_eq!(String::from_utf8(out).unwrap(), want);
    }

    #[test]
    fn read_bedgraph_roundtrip_plain() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("macs3rs_bdg_rt_{}.bdg", std::process::id()));
        let p = path.to_str().unwrap();
        std::fs::write(
            p,
            "track type=bedGraph\n\
             # a comment\n\
             chr1\t0\t100\t3.00000\n\
             chr1\t100\t200\t3.00000\n\
             chr1\t200\t300\t5.00000\n",
        )
        .unwrap();
        let t = read_bedgraph(p).unwrap();
        let c = t.get_data_by_chr(b"chr1").unwrap();
        // add_loc merges the two equal 3.0 intervals.
        assert_eq!(c.pos, vec![200, 300]);
        assert_eq!(c.val, vec![3.0_f32, 5.0]);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn read_bedgraph_roundtrip_gzip() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let dir = std::env::temp_dir();
        let path = dir.join(format!("macs3rs_bdg_rt_{}.bdg.gz", std::process::id()));
        let p = path.to_str().unwrap();
        let f = std::fs::File::create(p).unwrap();
        let mut enc = GzEncoder::new(f, Compression::default());
        enc.write_all(b"chr1\t0\t100\t1.50000\nchr2\t0\t50\t2.00000\n")
            .unwrap();
        enc.finish().unwrap();
        let t = read_bedgraph(p).unwrap();
        assert_eq!(t.get_data_by_chr(b"chr1").unwrap().val, vec![1.5_f32]);
        assert_eq!(t.get_data_by_chr(b"chr2").unwrap().val, vec![2.0_f32]);
        std::fs::remove_file(p).ok();
    }
}
