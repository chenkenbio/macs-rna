//! `bdgcmp` arguments and runner, mirroring `bin/macs3::add_bdgcmp_parser`
//! and `Commands/bdgcmp_cmd.py::run`.
//!
//! Field names are the MACS argparse `dest` names.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};
use clap::Args as ClapArgs;

use macs3_core::bedgraph_io::read_bedgraph;
use macs3_core::scoretrack::ScoreTrack2;

/// `macs3 bdgcmp` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-t/--tfile` (dest `tfile`): treatment bedGraph. Required.
    #[arg(short = 't', long = "tfile", required = true)]
    pub tfile: String,

    /// `-c/--cfile` (dest `cfile`): control bedGraph. Required.
    #[arg(short = 'c', long = "cfile", required = true)]
    pub cfile: String,

    /// `-S/--scaling-factor` (dest `sfactor`).
    #[arg(short = 'S', long = "scaling-factor", default_value_t = 1.0)]
    pub sfactor: f64,

    /// `-p/--pseudocount` (dest `pseudocount`).
    #[arg(short = 'p', long = "pseudocount", default_value_t = 0.0)]
    pub pseudocount: f64,

    /// `-m/--method` (dest `method`): one or more scoring methods.
    #[arg(
        short = 'm',
        long = "method",
        num_args = 1..,
        value_parser = ["ppois", "qpois", "subtract", "logFE", "FE", "logLR", "slogLR", "max"],
        default_values_t = [String::from("ppois")]
    )]
    pub method: Vec<String>,

    /// `--o-prefix` (dest `oprefix`): output filename prefix.
    #[arg(long = "o-prefix")]
    pub oprefix: Option<String>,

    /// `-o/--ofile` (dest `ofile`): explicit output filenames, one per method.
    #[arg(short = 'o', long = "ofile", num_args = 1..)]
    pub ofile: Vec<String>,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

/// Join `outdir` and a filename the way `os.path.join` does: an empty `outdir`
/// yields the bare filename (matching MACS, which defaults `outdir` to "").
fn join_outdir(outdir: &str, name: &str) -> String {
    if outdir.is_empty() {
        name.to_string()
    } else {
        Path::new(outdir).join(name).to_string_lossy().into_owned()
    }
}

/// Map a method name to the `change_score_method` byte
/// (`Commands/bdgcmp_cmd.py`).
fn method_byte(method: &str) -> u8 {
    match method {
        "ppois" => b'p',
        "qpois" => b'q',
        "subtract" => b'd',
        "logFE" => b'f',
        "FE" => b'F',
        "logLR" => b'l',
        "slogLR" => b's',
        "max" => b'M',
        other => unreachable!("clap restricts method values; got {other}"),
    }
}

/// Run `bdgcmp` (`Commands/bdgcmp_cmd.py::run`).
pub fn run(a: &Args) -> Result<()> {
    // opt_validate_bdgcmp: # of --ofile must == # of -m.
    if !a.ofile.is_empty() && a.ofile.len() != a.method.len() {
        bail!("The number and the order of arguments for --ofile must be the same as for -m.");
    }

    // pseudo_depth = 1 / scaling_factor (reciprocal trick to override SPMR).
    let scaling_factor = a.sfactor;
    let pseudo_depth = (1.0 / scaling_factor) as f32;

    // Read and build treatment / control bedGraphs.
    let tbtrack = read_bedgraph(&a.tfile)
        .with_context(|| format!("reading treatment bedGraph {}", a.tfile))?;
    let cbtrack = read_bedgraph(&a.cfile)
        .with_context(|| format!("reading control bedGraph {}", a.cfile))?;

    // Build ScoreTrackII via the make_ScoreTrackII_for_macs merge-walk.
    let mut sbtrack =
        ScoreTrack2::from_two_bedgraphs(&tbtrack, &cbtrack, pseudo_depth, pseudo_depth);

    // -S != 1: override SPMR by multiplying via normalization method 'M'.
    if (scaling_factor - 1.0).abs() > 1e-6 {
        sbtrack.change_normalization_method(b'M');
    }
    sbtrack.set_pseudocount(a.pseudocount as f32);

    // For each requested method (deduplicated, preserving first-seen order),
    // compute scores and write the bedGraph.
    let mut already: Vec<&str> = Vec::new();
    for (i, method) in a.method.iter().enumerate() {
        if already.contains(&method.as_str()) {
            continue;
        }
        already.push(method.as_str());

        let ofile = if !a.ofile.is_empty() {
            join_outdir(&a.outdir, &a.ofile[i])
        } else {
            let prefix = a
                .oprefix
                .as_deref()
                .unwrap_or("bdgcmp_output");
            join_outdir(&a.outdir, &format!("{prefix}_{method}.bdg"))
        };

        sbtrack.change_score_method(method_byte(method));

        let f = File::create(&ofile)
            .with_context(|| format!("creating output {ofile}"))?;
        let mut w = BufWriter::new(f);
        let name = format!("{}_Scores", method.to_uppercase());
        let desc = format!("Scores calculated by {}", method.to_uppercase());
        sbtrack
            .write_bedgraph(&mut w, &name, &desc, 3)
            .with_context(|| format!("writing {ofile}"))?;
        w.flush()?;
    }

    Ok(())
}
