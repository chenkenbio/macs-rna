//! `bdgpeakcall` arguments, mirroring `bin/macs3::add_bdgpeakcall_parser`.
//!
//! Field names are the MACS argparse `dest` names.
//!
//! Note `--no-trackline` is `store_false` with default `True` in MACS, i.e. the
//! effective `trackline` defaults to true and the flag turns it off. We model
//! that with a `--no-trackline` flag plus a `trackline()` accessor.
//!
//! Ports `Commands/bdgpeakcall_cmd.py`: read the score bedGraph, call
//! `call_peaks(cutoff, minlen, maxgap)`, and write narrowPeak.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use clap::Args as ClapArgs;

use macs3_core::bedgraph_io::read_bedgraph;

/// `macs3 bdgpeakcall` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-i/--ifile` (dest `ifile`): input score bedGraph. Required.
    #[arg(short = 'i', long = "ifile", required = true)]
    pub ifile: String,

    /// `-c/--cutoff` (dest `cutoff`).
    #[arg(short = 'c', long = "cutoff", default_value_t = 5.0)]
    pub cutoff: f64,

    /// `-l/--min-length` (dest `minlen`).
    #[arg(short = 'l', long = "min-length", default_value_t = 200)]
    pub minlen: i32,

    /// `-g/--max-gap` (dest `maxgap`).
    #[arg(short = 'g', long = "max-gap", default_value_t = 30)]
    pub maxgap: i32,

    /// `--call-summits` (dest `call_summits`).
    #[arg(long = "call-summits", default_value_t = false)]
    pub call_summits: bool,

    /// `--cutoff-analysis` (dest `cutoff_analysis`).
    #[arg(long = "cutoff-analysis", default_value_t = false)]
    pub cutoff_analysis: bool,

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

/// Run `bdgpeakcall` (`bdgpeakcall_cmd.py::run`).
pub fn run(a: &Args) -> anyhow::Result<()> {
    let btrack = read_bedgraph(&a.ifile)?;

    // cutoff_analysis path is not exercised by the macs-rna harness; only the
    // call-peaks path is ported.
    if a.cutoff_analysis {
        return Err(anyhow::anyhow!(
            "bdgpeakcall --cutoff-analysis not implemented in macs3-rs"
        ));
    }

    let peaks = btrack.call_peaks(
        a.cutoff as f32,
        a.minlen,
        a.maxgap,
        a.call_summits,
    );

    // Output filename / prefix selection (bdgpeakcall_cmd.py).
    let (oprefix, fname) = if let Some(ofile) = &a.ofile {
        (ofile.clone(), ofile.clone())
    } else {
        let oprefix = a
            .oprefix
            .clone()
            .ok_or_else(|| anyhow::anyhow!("either -o/--ofile or --o-prefix is required"))?;
        let fname = format!(
            "{}_c{:.1}_l{}_g{}_peaks.narrowPeak",
            oprefix, a.cutoff, a.minlen, a.maxgap
        );
        (oprefix, fname)
    };

    let out_path = Path::new(&a.outdir).join(&fname);
    let f = File::create(&out_path)?;
    let mut w = BufWriter::new(f);

    let name_prefix = format!("{oprefix}_narrowPeak");
    peaks.write_to_narrow_peak(
        &mut w,
        name_prefix.as_bytes(),
        oprefix.as_bytes(),
        "score",
        a.trackline(),
    )?;
    use std::io::Write;
    w.flush()?;
    Ok(())
}
