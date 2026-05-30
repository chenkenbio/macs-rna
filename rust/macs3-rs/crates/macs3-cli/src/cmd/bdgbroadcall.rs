//! `bdgbroadcall` arguments, mirroring `bin/macs3::add_bdgbroadcall_parser`.
//!
//! Field names are the MACS argparse `dest` names.
//!
//! Ports `Commands/bdgbroadcall_cmd.py`: read the score bedGraph, call
//! `call_broadpeaks(cutoffpeak, cutofflink, minlen, lvl1maxgap, lvl2maxgap)`,
//! and write the result as gappedPeak (BED12+3).

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use clap::Args as ClapArgs;

use macs3_core::bedgraph_io::read_bedgraph;

/// `macs3 bdgbroadcall` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-i/--ifile` (dest `ifile`): input score bedGraph. Required.
    #[arg(short = 'i', long = "ifile", required = true)]
    pub ifile: String,

    /// `-c/--cutoff-peak` (dest `cutoffpeak`).
    #[arg(short = 'c', long = "cutoff-peak", default_value_t = 2.0)]
    pub cutoffpeak: f64,

    /// `-C/--cutoff-link` (dest `cutofflink`).
    #[arg(short = 'C', long = "cutoff-link", default_value_t = 1.0)]
    pub cutofflink: f64,

    /// `-l/--min-length` (dest `minlen`).
    #[arg(short = 'l', long = "min-length", default_value_t = 200)]
    pub minlen: i32,

    /// `-g/--lvl1-max-gap` (dest `lvl1maxgap`).
    #[arg(short = 'g', long = "lvl1-max-gap", default_value_t = 30)]
    pub lvl1maxgap: i32,

    /// `-G/--lvl2-max-gap` (dest `lvl2maxgap`).
    #[arg(short = 'G', long = "lvl2-max-gap", default_value_t = 800)]
    pub lvl2maxgap: i32,

    /// `--no-trackline` (dest `trackline`, `store_false`, default True).
    #[arg(long = "no-trackline", default_value_t = false)]
    pub no_trackline: bool,

    /// `--o-prefix` (dest `oprefix`).
    #[arg(long = "o-prefix")]
    pub oprefix: Option<String>,

    /// `-o/--ofile` (dest `ofile`).
    #[arg(short = 'o', long = "ofile")]
    pub ofile: Option<String>,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

impl Args {
    /// Effective `trackline` value (`store_false` semantics: true unless
    /// `--no-trackline` was given).
    pub fn trackline(&self) -> bool {
        !self.no_trackline
    }
}

/// Run `bdgbroadcall` (`bdgbroadcall_cmd.py::run`).
pub fn run(a: &Args) -> anyhow::Result<()> {
    let btrack = read_bedgraph(&a.ifile)?;

    let bpeaks = btrack.call_broadpeaks(
        a.cutoffpeak as f32,
        a.cutofflink as f32,
        a.minlen,
        a.lvl1maxgap,
        a.lvl2maxgap,
    );

    // Output filename / prefix selection (bdgbroadcall_cmd.py).
    let (oprefix, fname) = if let Some(ofile) = &a.ofile {
        (ofile.clone(), ofile.clone())
    } else {
        let oprefix = a
            .oprefix
            .clone()
            .ok_or_else(|| anyhow::anyhow!("either -o/--ofile or --o-prefix is required"))?;
        let fname = format!(
            "{}_c{:.1}_C{:.2}_l{}_g{}_G{}_broad.bed12",
            oprefix, a.cutoffpeak, a.cutofflink, a.minlen, a.lvl1maxgap, a.lvl2maxgap
        );
        (oprefix, fname)
    };

    let out_path = Path::new(&a.outdir).join(&fname);
    let f = File::create(&out_path)?;
    let mut w = BufWriter::new(f);

    // write_to_gappedPeak defaults: name=b"peak", description=b"%s" (which
    // formats to b"peak"). Only name_prefix and score_column/trackline differ.
    let name_prefix = format!("{oprefix}_broadRegion");
    bpeaks.write_to_gapped_peak(
        &mut w,
        name_prefix.as_bytes(),
        b"peak",
        b"%s",
        "score",
        a.trackline(),
    )?;
    use std::io::Write;
    w.flush()?;
    Ok(())
}
