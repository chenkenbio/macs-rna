//! Input parsers: alignment/region files -> tracks.
//!
//! Ports `references/MACS/MACS3/IO/Parser.py` (`BEDParser`, `SAMParser`,
//! `BAMParser`, `BEDPEParser`, `BAMPEParser`). Single-end formats build an
//! [`FwTrack`] of 5'-end positions split by strand; paired-end formats build a
//! [`PeTrack`] of fragment intervals.
//!
//! Owner: Phase 1 agent A. These signatures are the contract that `cmd/pileup`,
//! `cmd/predictd`, and Phase 2 `cmd/callpeak` rely on.
//!
//! Bit-exactness: each parser matches `Parser.py` exactly for which coordinate
//! becomes the 5' end per strand, how malformed/zero-length records are skipped,
//! the `(+ -> 0, - -> 1)` strand encoding, and the read-mapping flag filters.
//! gzip input is detected by magic bytes; BAM is decoded from the raw
//! BGZF-decompressed byte stream (via `noodles-bgzf`) with the same binary
//! field offsets MACS uses in `bam_fw_binary_parse`/`get_references`.
//!
//! ## Final public API (consumed by `cmd/pileup`, `predictd`, callpeak)
//! ```ignore
//! pub fn parse_bed(path: &str, buffer_size: i64) -> io::Result<FwTrack>;
//! pub fn parse_sam(path: &str, buffer_size: i64) -> io::Result<FwTrack>;
//! pub fn parse_bam(path: &str, buffer_size: i64) -> io::Result<FwTrack>;
//! ```
//! None of these call `finalize()`; the caller finalizes (and `filter_dup`s)
//! after appending all input files, matching `load_tag_files_options`. `parse_bam`
//! additionally calls `FwTrack::set_rlengths` with the BAM-header lengths (as
//! `BAMParser.build_fwtrack` does), so its chromosome iteration order for
//! `pileup` becomes `sorted(valid) ++ sorted(missed)`.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufReader, Read};

use flate2::read::MultiGzDecoder;

use crate::constants::READ_BUFFER_SIZE;
use crate::track_fw::FwTrack;
use crate::track_pe::PeTrack;

/// gzip magic bytes (`1f 8b`).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Return `true` if the file begins with the gzip magic bytes. MACS tries
/// `gzip.open(...).read(10)` and falls back to plain on `IOError`.
fn is_gzipped(path: &str) -> io::Result<bool> {
    let mut f = File::open(path)?;
    let mut magic = [0u8; 2];
    match f.read_exact(&mut magic) {
        Ok(()) => Ok(magic == GZIP_MAGIC),
        Err(_) => Ok(false), // file shorter than 2 bytes -> not gzipped
    }
}

/// Open a (possibly gzipped) file and slurp its whole content into a `Vec<u8>`.
/// MACS streams in `READ_BUFFER_SIZE` chunks, but the line-splitting semantics
/// are identical to reading the entire decompressed stream and splitting on
/// `b"\n"`, so we read it all and split once.
fn read_all_maybe_gz(path: &str) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    if is_gzipped(path)? {
        let f = File::open(path)?;
        let mut dec = MultiGzDecoder::new(BufReader::with_capacity(READ_BUFFER_SIZE, f));
        dec.read_to_end(&mut buf)?;
    } else {
        let mut f = BufReader::with_capacity(READ_BUFFER_SIZE, File::open(path)?);
        f.read_to_end(&mut buf)?;
    }
    Ok(buf)
}

/// C `atoi` on a byte slice: skip leading ASCII whitespace, optional sign, then
/// consume digits until the first non-digit (0 if none). Matches MACS, which
/// calls libc `atoi` on BED/SAM integer columns.
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

/// Strip trailing ASCII whitespace, mirroring Python `bytes.rstrip()`.
fn rstrip(s: &[u8]) -> &[u8] {
    let mut end = s.len();
    while end > 0 && (s[end - 1] as char).is_ascii_whitespace() {
        end -= 1;
    }
    &s[..end]
}

/// Strip a trailing `.fa` suffix from a reference name, matching MACS's
/// `thisref[:thisref.rindex(b".fa")]` (only when `.fa` is present).
fn strip_fa(name: &[u8]) -> &[u8] {
    // rindex of b".fa"
    if name.len() < 3 {
        return name;
    }
    for start in (0..=name.len() - 3).rev() {
        if &name[start..start + 3] == b".fa" {
            return &name[..start];
        }
    }
    name
}

// ------------------------------------------------------------------
// BED
// ------------------------------------------------------------------

/// Parse a BED file (optionally gzipped) into an [`FwTrack`] of 5' ends.
///
/// Ports `Parser.py::BEDParser`. Leading `track`/`browser`/`#` lines are
/// skipped. For `+` strand the 5' end is the leftmost coordinate
/// (`atoi(col2)`, strand 0); for `-` it is the rightmost (`atoi(col3)`, strand
/// 1). Missing strand column defaults to `+`. Lines whose parsed position is
/// `< 0` or whose chromosome is empty are dropped (matching the `fpos < 0 or
/// not chromosome` filter in `build_fwtrack`). The track is **not** finalized.
pub fn parse_bed(path: &str, buffer_size: i64) -> io::Result<FwTrack> {
    let data = read_all_maybe_gz(path)?;
    let mut fwtrack = FwTrack::new(0, String::new(), buffer_size);

    let mut first_data_seen = false;
    for raw in data.split(|&b| b == b'\n') {
        // skip_first_commentlines: drop leading track/browser/# lines.
        if !first_data_seen {
            if is_bed_comment(raw) {
                continue;
            }
            first_data_seen = true;
        }
        if let Some((chrom, fpos, strand)) = bed_parse_line(raw) {
            if fpos < 0 || chrom.is_empty() {
                continue;
            }
            fwtrack.add_loc(chrom, fpos, strand);
        }
    }
    Ok(fwtrack)
}

/// `True` for BED comment/header lines (`track…`, `browser…`, `#…`) and for
/// empty lines, which the comment-skip loop in MACS also slides past.
fn is_bed_comment(line: &[u8]) -> bool {
    if line.is_empty() {
        return true;
    }
    line.starts_with(b"track") || line.starts_with(b"browser") || line[0] == b'#'
}

/// Parse one BED line into `(chrom, fivepos, strand)`. Returns `None` for an
/// empty line after rstrip (its split would yield `[b""]`, which the position
/// filter rejects anyway, but skipping avoids a spurious empty chromosome).
fn bed_parse_line(line: &[u8]) -> Option<(&[u8], i32, i32)> {
    let line = rstrip(line);
    if line.is_empty() {
        return None;
    }
    let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
    let chromname = fields[0];
    // MACS reads thisfields[5]; IndexError -> default plus strand at col1.
    match fields.get(5) {
        Some(&b"+") => Some((chromname, atoi(fields[1]), 0)),
        Some(&b"-") => Some((chromname, atoi(fields[2]), 1)),
        Some(_) => {
            // StrandFormatError in MACS; treat as a skipped malformed line.
            None
        }
        None => Some((chromname, atoi(fields[1]), 0)),
    }
}

// ------------------------------------------------------------------
// SAM
// ------------------------------------------------------------------

/// Parse a SAM file into an [`FwTrack`]. Ports `Parser.py::SAMParser`.
///
/// Header lines (`@…`) are skipped. The mapping filter is MACS's
/// `(bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2))`. For
/// the minus strand (`flag & 16`) the 5' end is `pos - 1 + sum(M/D/N/X/= ops)`
/// from the CIGAR; otherwise it is `pos - 1`. The track is **not** finalized.
pub fn parse_sam(path: &str, buffer_size: i64) -> io::Result<FwTrack> {
    let data = read_all_maybe_gz(path)?;
    let mut fwtrack = FwTrack::new(0, String::new(), buffer_size);

    let mut first_data_seen = false;
    for raw in data.split(|&b| b == b'\n') {
        if !first_data_seen {
            // skip_first_commentlines: drop leading @-prefixed header lines.
            if raw.is_empty() || raw[0] == b'@' {
                continue;
            }
            first_data_seen = true;
        }
        if let Some((chrom, fpos, strand)) = sam_parse_line(raw) {
            if fpos < 0 || chrom.is_empty() {
                continue;
            }
            fwtrack.add_loc(chrom, fpos, strand);
        }
    }
    Ok(fwtrack)
}

/// Parse one SAM alignment line. Returns `None` (mapped to a skipped line) for
/// empty lines or reads failing the mapping filter.
fn sam_parse_line(line: &[u8]) -> Option<(&[u8], i32, i32)> {
    let line = rstrip(line);
    if line.is_empty() {
        return None;
    }
    let fields: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
    if fields.len() < 6 {
        return None;
    }
    let thisref = fields[2];
    let bwflag = atoi(fields[1]);
    let cigar = fields[5];

    // (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2))
    if (bwflag & 2820) != 0
        || ((bwflag & 1) != 0 && ((bwflag & 136) != 0 || (bwflag & 2) == 0))
    {
        return None;
    }

    let (thisstrand, thisstart) = if (bwflag & 16) != 0 {
        // minus strand: shift by the summed alignment length from CIGAR.
        (1, atoi(fields[3]) - 1 + cigar_ref_len(cigar))
    } else {
        (0, atoi(fields[3]) - 1)
    };

    Some((strip_fa(thisref), thisstart, thisstrand))
}

/// Sum the lengths of CIGAR operations that consume the reference, matching
/// MACS's `findall(b"(\d+)[MDNX=]", CIGAR)` summed to `int`. Operations `S`,
/// `I`, `H`, `P` are ignored. A `*` (unavailable CIGAR) contributes 0.
fn cigar_ref_len(cigar: &[u8]) -> i32 {
    let mut total: i32 = 0;
    let mut num: i32 = 0;
    let mut have_num = false;
    for &c in cigar {
        if c.is_ascii_digit() {
            num = num * 10 + (c - b'0') as i32;
            have_num = true;
        } else {
            if have_num && matches!(c, b'M' | b'D' | b'N' | b'X' | b'=') {
                total += num;
            }
            num = 0;
            have_num = false;
        }
    }
    total
}

// ------------------------------------------------------------------
// BAM
// ------------------------------------------------------------------

/// Parse a BAM file (BGZF) into an [`FwTrack`]. Ports `Parser.py::BAMParser`.
///
/// The BGZF stream is decompressed with `noodles-bgzf`; the resulting raw BAM
/// byte stream is parsed with the exact field offsets MACS uses. References and
/// their lengths come from the BAM header (`get_references`), reads are filtered
/// with the same flag mask as `bam_fw_binary_parse`, and the minus-strand 5' end
/// is recovered by adding the reference-consuming CIGAR lengths. The track is
/// **not** finalized, but `set_rlengths` is applied (as in `build_fwtrack`).
pub fn parse_bam(path: &str, buffer_size: i64) -> io::Result<FwTrack> {
    // Decompress the whole BGZF stream into memory and parse the raw BAM bytes.
    let f = File::open(path)?;
    let mut reader = noodles_bgzf::io::Reader::new(BufReader::with_capacity(READ_BUFFER_SIZE, f));
    let mut data = Vec::new();
    reader.read_to_end(&mut data)?;

    let mut fwtrack = FwTrack::new(0, String::new(), buffer_size);

    // magic "BAM\x01"
    if data.len() < 4 || &data[0..4] != b"BAM\x01" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "File is not of a valid BAM format!",
        ));
    }
    let mut cur = 4usize;

    // l_text (i32) then text
    let l_text = read_i32(&data, cur)? as usize;
    cur += 4 + l_text;

    // n_ref (i32)
    let n_ref = read_i32(&data, cur)? as usize;
    cur += 4;

    let mut references: Vec<Vec<u8>> = Vec::with_capacity(n_ref);
    let mut rlengths: BTreeMap<Vec<u8>, i32> = BTreeMap::new();
    for _ in 0..n_ref {
        let l_name = read_i32(&data, cur)? as usize;
        cur += 4;
        // name is NUL-terminated; MACS uses fread(nlength)[:-1] to drop the NUL.
        let name = data[cur..cur + l_name - 1].to_vec();
        cur += l_name;
        let l_ref = read_i32(&data, cur)?;
        cur += 4;
        rlengths.insert(name.clone(), l_ref);
        references.push(name);
    }

    // alignment records
    while cur + 4 <= data.len() {
        let block_size = read_i32(&data, cur)? as usize;
        cur += 4;
        if cur + block_size > data.len() {
            break; // truncated final record; MACS stops on struct.error
        }
        let rec = &data[cur..cur + block_size];
        cur += block_size;

        if let Some((chrid, fpos, strand)) = bam_fw_binary_parse(rec) {
            // chrid == -1 means "skip"; otherwise index into references.
            if chrid >= 0 && (chrid as usize) < references.len() {
                fwtrack.add_loc(&references[chrid as usize], fpos, strand);
            }
        }
    }

    fwtrack.set_rlengths(&rlengths);
    Ok(fwtrack)
}

/// Read a little-endian `i32` at `off`, erroring if out of range.
fn read_i32(data: &[u8], off: usize) -> io::Result<i32> {
    if off + 4 > data.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "BAM: truncated integer field",
        ));
    }
    Ok(i32::from_le_bytes([
        data[off],
        data[off + 1],
        data[off + 2],
        data[off + 3],
    ]))
}

/// Parse one BAM SE alignment record body, mirroring `bam_fw_binary_parse`.
///
/// `rec` is the record bytes *after* the 4-byte block_size prefix. Returns
/// `(refID, fivepos, strand)` or `None` (== MACS's `(-1, -1, -1)` skip) when the
/// read fails the mapping filter. For the minus strand (`flag & 16`) the start
/// is advanced by the reference-consuming CIGAR ops (`M/D/N/=/X`).
fn bam_fw_binary_parse(rec: &[u8]) -> Option<(i32, i32, i32)> {
    if rec.len() < 16 {
        return None;
    }
    let n_cigar_op = u16::from_le_bytes([rec[12], rec[13]]) as usize;
    let bwflag = u16::from_le_bytes([rec[14], rec[15]]) as i32;

    // (bwflag & 2820) or (bwflag & 1 and (bwflag & 136 or not bwflag & 2))
    if (bwflag & 2820) != 0
        || ((bwflag & 1) != 0 && ((bwflag & 136) != 0 || (bwflag & 2) == 0))
    {
        return None;
    }

    let thisref = i32::from_le_bytes([rec[0], rec[1], rec[2], rec[3]]);
    let mut thisstart = i32::from_le_bytes([rec[4], rec[5], rec[6], rec[7]]);

    let thisstrand = if (bwflag & 16) != 0 {
        let l_read_name = rec[8] as usize;
        // CIGAR begins at byte offset 32 + l_read_name.
        let cigar_off = 32 + l_read_name;
        for k in 0..n_cigar_op {
            let off = cigar_off + k * 4;
            if off + 4 > rec.len() {
                break;
            }
            let cigar_code =
                u32::from_le_bytes([rec[off], rec[off + 1], rec[off + 2], rec[off + 3]]);
            let op = cigar_code & 15;
            // CIGAR op M/D/N/=/X consume the reference (op codes 0,2,3,7,8).
            if matches!(op, 0 | 2 | 3 | 7 | 8) {
                thisstart += (cigar_code >> 4) as i32;
            }
        }
        1
    } else {
        0
    };

    Some((thisref, thisstart, thisstrand))
}

// ------------------------------------------------------------------
// Paired-end (low priority — stubbed)
// ------------------------------------------------------------------

/// Parse a BEDPE file into a [`PeTrack`]. Ports `Parser.py::BEDPEParser`. STUB.
pub fn parse_bedpe(_path: &str, _buffer_size: i64) -> io::Result<PeTrack> {
    todo!("Phase 1A (low priority): port Parser.py::BEDPEParser")
}

/// Parse a BAMPE file into a [`PeTrack`]. Ports `Parser.py::BAMPEParser`. STUB.
pub fn parse_bampe(_path: &str, _buffer_size: i64) -> io::Result<PeTrack> {
    todo!("Phase 1A (low priority): port Parser.py::BAMPEParser")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoi_matches_c() {
        assert_eq!(atoi(b"123"), 123);
        assert_eq!(atoi(b"  -45abc"), -45);
        assert_eq!(atoi(b"+7"), 7);
        assert_eq!(atoi(b""), 0);
        assert_eq!(atoi(b"x"), 0);
    }

    #[test]
    fn strip_fa_suffix() {
        assert_eq!(strip_fa(b"chr1.fa"), b"chr1");
        assert_eq!(strip_fa(b"chr1"), b"chr1");
        assert_eq!(strip_fa(b"chr1.fasta"), b"chr1"); // rindex of ".fa" hits
    }

    #[test]
    fn cigar_ref_len_only_ref_consuming() {
        // 10M5I20M -> 30; soft clips / insertions ignored.
        assert_eq!(cigar_ref_len(b"10M5I20M"), 30);
        assert_eq!(cigar_ref_len(b"50M"), 50);
        assert_eq!(cigar_ref_len(b"5S30M2D3N"), 35);
        assert_eq!(cigar_ref_len(b"*"), 0);
    }

    #[test]
    fn bed_plus_minus_five_ends() {
        // + strand -> leftmost (col1, strand 0); - strand -> rightmost (col2, 1).
        let (c, p, s) = bed_parse_line(b"chr1\t100\t200\t.\t.\t+").unwrap();
        assert_eq!((c, p, s), (&b"chr1"[..], 100, 0));
        let (c, p, s) = bed_parse_line(b"chr1\t100\t200\t.\t.\t-").unwrap();
        assert_eq!((c, p, s), (&b"chr1"[..], 200, 1));
        // missing strand -> default plus
        let (c, p, s) = bed_parse_line(b"chr1\t100\t200").unwrap();
        assert_eq!((c, p, s), (&b"chr1"[..], 100, 0));
    }

    #[test]
    fn sam_minus_strand_uses_cigar() {
        // flag=0 -> plus, start = pos-1
        let (_c, p, s) = sam_parse_line(
            b"r1\t0\tchr1\t101\t60\t50M\t*\t0\t0\tACGT\tIIII",
        )
        .unwrap();
        assert_eq!((p, s), (100, 0));
        // flag=16 -> minus, start = pos-1 + cigar_ref_len(50M) = 100+50 = 150
        let (_c, p, s) = sam_parse_line(
            b"r2\t16\tchr1\t101\t60\t50M\t*\t0\t0\tACGT\tIIII",
        )
        .unwrap();
        assert_eq!((p, s), (150, 1));
        // unmapped (flag & 4) -> skipped
        assert!(sam_parse_line(b"r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").is_none());
    }
}
