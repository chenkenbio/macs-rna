//! `bdgdiff` arguments and runner, mirroring `bin/macs3::add_bdgdiff_parser`
//! and `Commands/bdgdiff_cmd.py::run`.
//!
//! Field names are the MACS argparse `dest` names. The three BED5 outputs are
//! written directly here (matching `PeakIO._to_bed` with `score_column="score"`
//! and `trackline=False`) so this command stays independent of agent B's
//! `PeakIO::write_to_bed`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use clap::Args as ClapArgs;

use macs3_core::bedgraph_io::read_bedgraph;
use macs3_core::fmt::format_g6;
use macs3_core::peak_io::PeakIO;
use macs3_core::twocond::TwoConditionScores;

/// `macs3 bdgdiff` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `--t1` (dest `t1bdg`): condition-1 treatment bedGraph. Required.
    #[arg(long = "t1", required = true)]
    pub t1bdg: String,

    /// `--c1` (dest `c1bdg`): condition-1 control bedGraph. Required.
    #[arg(long = "c1", required = true)]
    pub c1bdg: String,

    /// `--t2` (dest `t2bdg`): condition-2 treatment bedGraph. Required.
    #[arg(long = "t2", required = true)]
    pub t2bdg: String,

    /// `--c2` (dest `c2bdg`): condition-2 control bedGraph. Required.
    #[arg(long = "c2", required = true)]
    pub c2bdg: String,

    /// `-C/--cutoff` (dest `cutoff`): log10 likelihood-ratio cutoff.
    #[arg(short = 'C', long = "cutoff", default_value_t = 3.0)]
    pub cutoff: f64,

    /// `-l/--min-len` (dest `minlen`).
    #[arg(short = 'l', long = "min-len", default_value_t = 200)]
    pub minlen: i32,

    /// `-g/--max-gap` (dest `maxgap`).
    #[arg(short = 'g', long = "max-gap", default_value_t = 100)]
    pub maxgap: i32,

    /// `--d1/--depth1` (dest `depth1`).
    #[arg(long = "d1", visible_alias = "depth1", default_value_t = 1.0)]
    pub depth1: f64,

    /// `--d2/--depth2` (dest `depth2`).
    #[arg(long = "d2", visible_alias = "depth2", default_value_t = 1.0)]
    pub depth2: f64,

    /// `--o-prefix` (dest `oprefix`).
    #[arg(long = "o-prefix")]
    pub oprefix: Option<String>,

    /// `-o/--ofile` (dest `ofile`): exactly three output filenames.
    #[arg(short = 'o', long = "ofile", num_args = 3)]
    pub ofile: Vec<String>,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

/// Join `outdir` and a filename like `os.path.join` (empty `outdir` -> bare).
fn join_outdir(outdir: &str, name: &str) -> String {
    if outdir.is_empty() {
        name.to_string()
    } else {
        Path::new(outdir).join(name).to_string_lossy().into_owned()
    }
}

/// Run `bdgdiff` (`Commands/bdgdiff_cmd.py::run`).
pub fn run(a: &Args) -> Result<()> {
    if a.maxgap >= a.minlen {
        anyhow::bail!(
            "MAXGAP should be smaller than MINLEN! Your input is MAXGAP = {} and MINLEN = {}",
            a.maxgap,
            a.minlen
        );
    }

    let t1 = read_bedgraph(&a.t1bdg).with_context(|| format!("reading {}", a.t1bdg))?;
    let c1 = read_bedgraph(&a.c1bdg).with_context(|| format!("reading {}", a.c1bdg))?;
    let t2 = read_bedgraph(&a.t2bdg).with_context(|| format!("reading {}", a.t2bdg))?;
    let c2 = read_bedgraph(&a.c2bdg).with_context(|| format!("reading {}", a.c2bdg))?;

    // Down-scale the larger condition to the smaller one (bdgdiff_cmd.py).
    let (mut depth1, mut depth2) = (a.depth1, a.depth2);
    if depth1 > depth2 {
        depth1 = depth2 / depth1;
        depth2 = 1.0;
    } else if depth1 < depth2 {
        depth2 = depth1 / depth2;
        depth1 = 1.0;
    } else {
        depth1 = 1.0;
        depth2 = 1.0;
    }

    // Build differential scores (pseudocount default 0.01).
    let mut tcs = TwoConditionScores::new(
        t1,
        c1,
        t2,
        c2,
        depth1 as f32,
        depth2 as f32,
        0.01,
    );
    tcs.build();
    let (cat1, cat2, cat3) =
        tcs.call_peaks(a.cutoff as f32, a.minlen, a.maxgap);

    // Resolve the three output paths and per-file peak-name prefixes.
    let (ofiles, name_prefixes): ([String; 3], [String; 3]) = if !a.ofile.is_empty() {
        (
            [
                join_outdir(&a.outdir, &a.ofile[0]),
                join_outdir(&a.outdir, &a.ofile[1]),
                join_outdir(&a.outdir, &a.ofile[2]),
            ],
            // _to_bed does `name_prefix % name`; the raw ofile string carries no
            // `%s`, so the `%` fails and the prefix is the ofile string as-is.
            [a.ofile[0].clone(), a.ofile[1].clone(), a.ofile[2].clone()],
        )
    } else {
        let prefix = a.oprefix.as_deref().unwrap_or("diffpeak");
        (
            [
                join_outdir(&a.outdir, &format!("{prefix}_c{:.1}_cond1.bed", a.cutoff)),
                join_outdir(&a.outdir, &format!("{prefix}_c{:.1}_cond2.bed", a.cutoff)),
                join_outdir(&a.outdir, &format!("{prefix}_c{:.1}_common.bed", a.cutoff)),
            ],
            [
                format!("{prefix}_cond1_"),
                format!("{prefix}_cond2_"),
                format!("{prefix}_common_"),
            ],
        )
    };

    // Track-line metadata per category, matching bdgdiff_cmd.py's
    // write_to_bed(name=..., description=...) calls.
    write_bed5(
        &ofiles[0],
        &cat1,
        &name_prefixes[0],
        "condition 1",
        "unique regions in condition 1",
    )?;
    write_bed5(
        &ofiles[1],
        &cat2,
        &name_prefixes[1],
        "condition 2",
        "unique regions in condition 2",
    )?;
    write_bed5(
        &ofiles[2],
        &cat3,
        &name_prefixes[2],
        "common",
        "common regions in both conditions",
    )?;

    Ok(())
}

/// Alphabetical subpeak label for a zero-based index (`PeakIO.subpeak_letters`):
/// `a`, `b`, ..., `z`, `aa`, `ab`, ... .
fn subpeak_letters(i: usize) -> String {
    if i < 26 {
        ((b'a' + i as u8) as char).to_string()
    } else {
        let mut s = subpeak_letters(i / 26);
        s.push((b'a' + (i % 26) as u8) as char);
        s
    }
}

/// Write a [`PeakIO`] as a BED5 file, replicating `PeakIO._to_bed` with
/// `score_column="score"` and `trackline=True` (bdgdiff's default).
///
/// A UCSC `track` header is written first
/// (`track name="<name> (peaks)" description="<description>" visibility=1`),
/// with quotes in `name`/`description` backslash-escaped. Chromosomes are
/// emitted bytewise-sorted; within a chromosome, peaks sharing the same `end`
/// are grouped (a global `n_peak` counter increments once per group, and grouped
/// peaks get an alphabetical suffix). The score is the peak score formatted with
/// C `%.6g`.
fn write_bed5(
    path: &str,
    peaks: &PeakIO,
    name_prefix: &str,
    track_name: &str,
    track_desc: &str,
) -> Result<()> {
    let f = File::create(path).with_context(|| format!("creating {path}"))?;
    let mut w = BufWriter::new(f);

    // _to_bed trackline: name escaped, description used as-is (the `description
    // % name` interpolation fails for these literal descriptions, so the raw
    // string is kept).
    let name_esc = track_name.replace('"', "\\\"");
    let desc_esc = track_desc.replace('"', "\\\"");
    writeln!(
        w,
        "track name=\"{name_esc} (peaks)\" description=\"{desc_esc}\" visibility=1"
    )?;

    let mut n_peak: i64 = 0;
    for chrom in peaks.get_chr_names() {
        let chrom_peaks = match peaks.peaks_by_chr(chrom) {
            Some(p) => p,
            None => continue,
        };
        // Group consecutive peaks by equal `end` (itertools.groupby semantics:
        // grouping is over *consecutive* equal keys in insertion order).
        let mut i = 0usize;
        while i < chrom_peaks.len() {
            let end = chrom_peaks[i].end;
            let mut j = i;
            while j < chrom_peaks.len() && chrom_peaks[j].end == end {
                j += 1;
            }
            n_peak += 1;
            let group = &chrom_peaks[i..j];
            let chrom_str = String::from_utf8_lossy(chrom);
            if group.len() > 1 {
                for (k, peak) in group.iter().enumerate() {
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}{}{}\t{}",
                        chrom_str,
                        peak.start,
                        peak.end,
                        name_prefix,
                        n_peak,
                        subpeak_letters(k),
                        format_g6(peak.score),
                    )?;
                }
            } else {
                let peak = &group[0];
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}{}\t{}",
                    chrom_str,
                    peak.start,
                    peak.end,
                    name_prefix,
                    n_peak,
                    format_g6(peak.score),
                )?;
            }
            i = j;
        }
    }
    w.flush()?;
    Ok(())
}
