//! `pileup` arguments, mirroring `bin/macs3::add_pileup_parser`.
//!
//! Field names are the MACS argparse `dest` names. Phase 1A fills `run`.
//!
//! Note: MACS's `pileup` exposes `-B/--both-direction` (dest `bothdirection`)
//! rather than a `--bdg` flag, and `--extsize`; there is no `--shift` /
//! `--scale-to` on this subcommand in 3.0.4. We keep the real flag set.

use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, Context};
use clap::Args as ClapArgs;

use macs3_core::bedgraph_io::write_pv_exact;
use macs3_core::parser::{parse_bam, parse_bed, parse_sam};
use macs3_core::pileup::single_end_pileup;
use macs3_core::track_fw::FwTrack;

/// `macs3 pileup` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-i/--ifile` (dest `ifile`): input file(s), 1+.
    #[arg(short = 'i', long = "ifile", required = true, num_args = 1..)]
    pub ifile: Vec<String>,

    /// `-o/--ofile` (dest `outputfile`): output bedGraph. Required.
    #[arg(short = 'o', long = "ofile", required = true)]
    pub outputfile: String,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `-f/--format` (dest `format`).
    #[arg(short = 'f', long = "format", default_value = "AUTO")]
    pub format: String,

    /// `-B/--both-direction` (dest `bothdirection`): extend in both directions.
    #[arg(short = 'B', long = "both-direction", default_value_t = false)]
    pub bothdirection: bool,

    /// `--extsize` (dest `extsize`).
    #[arg(long = "extsize", default_value_t = 200)]
    pub extsize: i32,

    /// `--buffer-size` (dest `buffer_size`).
    #[arg(long = "buffer-size", default_value_t = 100000)]
    pub buffer_size: i64,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

/// Run `pileup`, porting `Commands/pileup_cmd.py::run` (single-end branch) plus
/// `Signal/Pileup.py::pileup_and_write_se`.
///
/// Pipeline:
/// 1. Pick a parser from `-f` (BED / SAM / BAM; PE formats are unsupported).
/// 2. Build an `FwTrack` from the first input, append the rest, then finalize.
/// 3. Derive 5'/3' shifts (`halfextension=False`): default extends only to the
///    3' side (`five=0, three=extsize`); `-B` extends both sides with
///    `d = extsize*2`, `five=d//2, three=d-d//2`.
/// 4. For each chromosome in `get_rlengths()` insertion order, pile up with
///    `single_end_pileup(plus, minus, five, three, 0, rlength, 1.0, 0.0)` and
///    append via `write_pv_exact`.
pub fn run(a: &Args) -> anyhow::Result<()> {
    let fmt = a.format.to_uppercase();

    if matches!(fmt.as_str(), "BAMPE" | "BEDPE" | "FRAG") {
        bail!(
            "pileup format {fmt:?} (paired-end) is not supported by macs3-rs yet"
        );
    }

    // Resolve output path = outdir/outputfile and remove any existing file
    // (matches `os.path.join(outdir, outputfile)` + `os.unlink`).
    let outfile = if a.outdir.is_empty() {
        a.outputfile.clone()
    } else {
        Path::new(&a.outdir)
            .join(&a.outputfile)
            .to_string_lossy()
            .into_owned()
    };
    if Path::new(&outfile).is_file() {
        let _ = std::fs::remove_file(&outfile);
    }

    if a.extsize <= 0 {
        bail!("--extsize must > 0!");
    }

    // 1+2: parse the first file, append the rest, finalize.
    let mut treat = parse_one(&fmt, &a.ifile[0], a.buffer_size)?;
    for tfile in &a.ifile[1..] {
        let extra = parse_one(&fmt, tfile, a.buffer_size)?;
        append_fwtrack(&mut treat, extra);
    }
    treat.finalize();

    // 3: shift values (halfextension=False).
    let (five_shift, three_shift) = if a.bothdirection {
        let d = a.extsize * 2;
        (d / 2, d - d / 2)
    } else {
        (0, a.extsize)
    };

    // 4: pile up each chromosome and write. Chromosome iteration order is
    // `list(get_rlengths().keys())` (insertion order), NOT bytewise-sorted.
    let chrom_order: Vec<Vec<u8>> = treat
        .rlengths_iter_order()
        .into_iter()
        .map(|c| c.to_vec())
        .collect();
    let rlengths = treat.get_rlengths().clone();

    let file = File::create(&outfile)
        .with_context(|| format!("failed to create output file {outfile}"))?;
    let mut out = BufWriter::new(file);

    for chrom in &chrom_order {
        let Some(loc) = treat.get_locations_by_chr(chrom) else {
            continue;
        };
        let rlength = rlengths.get(chrom).copied().unwrap_or(macs3_core::track_fw::INT_MAX);
        let pv = single_end_pileup(
            &loc.plus,
            &loc.minus,
            five_shift,
            three_shift,
            0,
            rlength,
            1.0,
            0.0,
        );
        // `append` is irrelevant here: we hold a single handle and write each
        // chromosome's bytes sequentially, equivalent to the C append-per-chrom.
        write_pv_exact(chrom, &pv, &mut out, true)?;
    }

    out.flush()?;
    Ok(())
}

/// Build an `FwTrack` from a single input file according to `fmt`.
///
/// For an explicit format we dispatch directly. For `AUTO` (not exercised by
/// the diff harness, which always passes `-f BED`) we sniff: a raw `BAM\x01`
/// magic -> BAM; a leading `@` -> SAM; a gzip magic -> try BAM (BGZF), falling
/// back to BED if the decompressed stream is not a valid BAM; otherwise BED.
fn parse_one(fmt: &str, path: &str, buffer_size: i64) -> anyhow::Result<FwTrack> {
    if fmt != "AUTO" {
        return match fmt {
            "BED" => parse_bed(path, buffer_size).map_err(Into::into),
            "SAM" => parse_sam(path, buffer_size).map_err(Into::into),
            "BAM" => parse_bam(path, buffer_size).map_err(Into::into),
            other => bail!("pileup format {other:?} is not supported by macs3-rs"),
        };
    }

    // AUTO sniff via leading bytes.
    let mut f = File::open(path)?;
    let mut head = [0u8; 4];
    let n = f.read(&mut head).unwrap_or(0);
    if n >= 4 && &head == b"BAM\x01" {
        return parse_bam(path, buffer_size).map_err(Into::into);
    }
    if n >= 1 && head[0] == b'@' {
        return parse_sam(path, buffer_size).map_err(Into::into);
    }
    if n >= 2 && head[0] == 0x1f && head[1] == 0x8b {
        // gzip: try BAM (BGZF) first; fall back to BED on invalid-BAM.
        match parse_bam(path, buffer_size) {
            Ok(t) => return Ok(t),
            Err(e) if e.kind() == std::io::ErrorKind::InvalidData => {
                return parse_bed(path, buffer_size).map_err(Into::into);
            }
            Err(e) => return Err(e.into()),
        }
    }
    parse_bed(path, buffer_size).map_err(Into::into)
}

/// Merge `extra`'s locations into `dst`, mirroring `append_fwtrack` (positions
/// are appended without finalizing). Used for multiple `-i` inputs.
fn append_fwtrack(dst: &mut FwTrack, extra: FwTrack) {
    for chrom in extra.get_chr_names() {
        if let Some(loc) = extra.get_locations_by_chr(chrom) {
            for &p in &loc.plus {
                dst.add_loc(chrom, p, 0);
            }
            for &p in &loc.minus {
                dst.add_loc(chrom, p, 1);
            }
        }
    }
}
