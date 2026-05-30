//! `bdgopt` arguments, mirroring `bin/macs3::add_bdgopt_parser`.
//!
//! Field names are the MACS argparse `dest` names.
//!
//! Ports `Commands/bdgopt_cmd.py`: read the input bedGraph, apply the chosen
//! score-column modification (`multiply`/`add`/`max`/`min` via `apply_func`, or
//! `p2q`), and write the modified track via the BedGraphIO writer (every
//! interval emitted, `%.5f`, with a UCSC trackline). Chromosomes are written
//! bytewise-sorted (`get_chr_names`).

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use clap::Args as ClapArgs;

use macs3_core::bedgraph_io::{
    bedgraph_trackline, read_bedgraph, write_bedgraph_predicate, AlwaysBreak,
};

/// `macs3 bdgopt` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-i/--ifile` (dest `ifile`): input bedGraph. Required.
    #[arg(short = 'i', long = "ifile", required = true)]
    pub ifile: String,

    /// `-m/--method` (dest `method`): score-column modification method.
    #[arg(
        short = 'm',
        long = "method",
        value_parser = ["multiply", "add", "p2q", "max", "min"],
        default_value = "p2q"
    )]
    pub method: String,

    /// `-p/--extra-param` (dest `extraparam`): extra parameter(s) for `method`.
    #[arg(short = 'p', long = "extra-param", num_args = 0..)]
    pub extraparam: Vec<f64>,

    /// `-o/--ofile` (dest `ofile`): output bedGraph. Required.
    #[arg(short = 'o', long = "ofile", required = true)]
    pub ofile: String,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

/// Run `bdgopt` (`bdgopt_cmd.py::run`).
pub fn run(a: &Args) -> anyhow::Result<()> {
    let mut btrack = read_bedgraph(&a.ifile)?;

    let method = a.method.to_lowercase();
    if method == "p2q" {
        btrack.p2q();
    } else {
        // multiply/add require an extra parameter; max/min use it too.
        // opt_validate_bdgopt enforces extraparam presence for multiply/add.
        let extraparam = *a
            .extraparam
            .first()
            .ok_or_else(|| anyhow::anyhow!("Need EXTRAPARAM for method {}!", method))?
            as f32;
        match method.as_str() {
            "multiply" => btrack.apply_func(|x| x * extraparam),
            "add" => btrack.apply_func(|x| x + extraparam),
            "max" => btrack.apply_func(|x| if x > extraparam { x } else { extraparam }),
            "min" => btrack.apply_func(|x| if x < extraparam { x } else { extraparam }),
            other => return Err(anyhow::anyhow!("Invalid method: {other}")),
        }
    }

    let out_path = Path::new(&a.outdir).join(&a.ofile);
    let f = File::create(&out_path)?;
    let mut w = BufWriter::new(f);

    // BedGraphIO.write_bedGraph: trackline on by default, name/description as
    // built by bdgopt_cmd.py.
    let upper = a.method.to_uppercase();
    let name = format!("{upper}_modified_scores");
    let description = format!("Scores calculated by {upper}");
    let trackline = bedgraph_trackline(&name, &description);
    use std::io::Write;
    w.write_all(trackline.as_bytes())?;

    for chrom in btrack.get_chr_names() {
        let c = btrack.get_data_by_chr(chrom).unwrap();
        write_bedgraph_predicate(chrom, &c.pos, &c.val, &AlwaysBreak, &mut w)?;
    }
    w.flush()?;
    Ok(())
}
