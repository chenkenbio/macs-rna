//! `predictd` arguments + driver, mirroring `bin/macs3::add_predictd_parser`
//! and `MACS3/Commands/predictd_cmd.py`.
//!
//! Field names are the MACS argparse `dest` names. `run` ports `predictd_cmd.run`:
//! parse the input by `-f`, finalize the [`FwTrack`], build a [`PeakModel`], and
//! print `# predicted fragment length is N bps` to stderr (the line macs-rna
//! greps). The R model file is written best-effort.

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use clap::Args as ClapArgs;

use macs3_core::constants::{effective_gsize, MAX_PAIRNUM};
use macs3_core::parser::{parse_bam, parse_bed, parse_sam};
use macs3_core::peakmodel::{PeakModel, PeakModelOptions};
use macs3_core::track_fw::FwTrack;

/// `macs3 predictd` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-i/--ifile` (dest `ifile`): input file(s), 1+.
    #[arg(short = 'i', long = "ifile", required = true, num_args = 1..)]
    pub ifile: Vec<String>,

    /// `-f/--format` (dest `format`).
    #[arg(short = 'f', long = "format", default_value = "AUTO")]
    pub format: String,

    /// `-g/--gsize` (dest `gsize`).
    #[arg(short = 'g', long = "gsize", default_value = "hs")]
    pub gsize: String,

    /// `-s/--tsize` (dest `tsize`).
    #[arg(short = 's', long = "tsize")]
    pub tsize: Option<i32>,

    /// `--bw` (dest `bw`).
    #[arg(long = "bw", default_value_t = 300)]
    pub bw: i32,

    /// `--d-min` (dest `d_min`).
    #[arg(long = "d-min", default_value_t = 20)]
    pub d_min: i32,

    /// `-m/--mfold` (dest `mfold`): exactly two values (lower upper).
    #[arg(short = 'm', long = "mfold", num_args = 2, default_values_t = [5, 50])]
    pub mfold: Vec<i32>,

    /// `--rfile` (dest `rfile`): R model output filename.
    #[arg(long = "rfile", default_value = "predictd_model.R")]
    pub rfile: String,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `--buffer-size` (dest `buffer_size`).
    #[arg(long = "buffer-size", default_value_t = 100000)]
    pub buffer_size: i64,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

/// Run `predictd` (`predictd_cmd.run`).
pub fn run(a: &Args) -> anyhow::Result<()> {
    // bin/macs3 creates --outdir before dispatch if it does not exist.
    if !a.outdir.is_empty() && !Path::new(&a.outdir).exists() {
        std::fs::create_dir_all(&a.outdir)?;
    }

    // --- gsize: efgsize shortcut, else float (OptValidator.opt_validate_predictd).
    let gsize = match effective_gsize(&a.gsize) {
        Some(v) => v,
        None => a.gsize.parse::<f64>().map_err(|_| {
            anyhow::anyhow!("Error when interpreting --gsize option: {}", a.gsize)
        })?,
    };

    // --- format -> parser. AUTO sniffs BAM then BED (guess_parser order).
    let format = a.format.to_uppercase();

    // PE_MODE: BAMPE / BEDPE. predictd only prints the average insert length.
    if format == "BAMPE" || format == "BEDPE" {
        // PE paths are low priority (parser is stubbed); we cannot compute the
        // average template length yet. Match the MACS message shape on stderr.
        anyhow::bail!("predictd: paired-end formats (BAMPE/BEDPE) are not supported yet");
    }

    // --- mfold (lower, upper).
    if a.mfold.len() != 2 {
        anyhow::bail!("--mfold requires exactly two integers");
    }
    let lmfold = a.mfold[0];
    let umfold = a.mfold[1];
    if lmfold > umfold {
        anyhow::bail!("Upper limit of mfold should be greater than lower limit!");
    }
    if a.d_min < 0 {
        anyhow::bail!("Minimum fragment size shouldn't be negative!");
    }

    // --- 1. read tag files (load_tag_files_options): parse each input, append,
    //        then finalize ONCE. predictd does NOT call filter_dup.
    eprintln!("# read alignment files...");
    let mut treat = load_tag_files(a, &format)?;
    treat.finalize();
    let t0 = treat.total;
    eprintln!("# total tags in alignment file: {}", t0);

    // --- 2. build model.
    eprintln!("# Build Peak Model...");
    let opt = PeakModelOptions {
        gsize,
        umfold,
        lmfold,
        bw: a.bw,
        d_min: a.d_min,
    };
    let mut peakmodel = PeakModel::new(opt, MAX_PAIRNUM);

    match peakmodel.build(&treat) {
        Ok(()) => {
            // The must-have line macs-rna greps (wording identical to MACS).
            eprintln!("# predicted fragment length is {} bps", peakmodel.d);
            let alts: Vec<String> =
                peakmodel.alternative_d.iter().map(|d| d.to_string()).collect();
            eprintln!(
                "# alternative fragment length(s) may be {} bps",
                alts.join(",")
            );

            // Best-effort R model file (OptValidator: modelR = outdir/rfile).
            let model_r = if a.outdir.is_empty() {
                a.rfile.clone()
            } else {
                Path::new(&a.outdir)
                    .join(&a.rfile)
                    .to_string_lossy()
                    .into_owned()
            };
            if let Err(e) = write_model_r(&peakmodel, &model_r) {
                eprintln!("# warning: failed to write R model {}: {}", model_r, e);
            }
        }
        Err(_) => {
            eprintln!("# Can't find enough pairs of symmetric peaks to build model!");
        }
    }

    Ok(())
}

/// Port of `load_tag_files_options`: parse `ifile[0]`, append the rest, return
/// the (not-yet-finalized) track. AUTO format sniffs BAM then BED per input.
fn load_tag_files(a: &Args, format: &str) -> anyhow::Result<FwTrack> {
    let mut treat = parse_one(&a.ifile[0], format, a.buffer_size)?;
    for f in &a.ifile[1..] {
        let extra = parse_one(f, format, a.buffer_size)?;
        append_fwtrack(&mut treat, extra);
    }
    Ok(treat)
}

/// Dispatch one input file to the right parser. For AUTO we replicate
/// `guess_parser`'s BAM-before-BED ordering: if the file is gzip-magic'd we try
/// BAM first (it cleanly errors on non-BAM BGZF/gzip content) and fall back to
/// BED; a plain (non-gzip) file goes straight to BED. Only the formats macs-rna
/// feeds predictd (BED/SAM/BAM/AUTO) are handled.
fn parse_one(path: &str, format: &str, buffer_size: i64) -> anyhow::Result<FwTrack> {
    match format {
        "BED" => Ok(parse_bed(path, buffer_size)?),
        "SAM" => Ok(parse_sam(path, buffer_size)?),
        "BAM" => Ok(parse_bam(path, buffer_size)?),
        "AUTO" => {
            if is_gzip_magic(path) {
                // gzip member: could be BAM (BGZF) or a gzipped BED. Try BAM,
                // fall back to BED on any decode/format error.
                match parse_bam(path, buffer_size) {
                    Ok(t) => Ok(t),
                    Err(_) => Ok(parse_bed(path, buffer_size)?),
                }
            } else {
                Ok(parse_bed(path, buffer_size)?)
            }
        }
        other => anyhow::bail!("Format \"{}\" cannot be recognized for predictd!", other),
    }
}

/// Return `true` if `path` begins with the gzip magic bytes (`1f 8b`).
fn is_gzip_magic(path: &str) -> bool {
    let mut head = [0u8; 2];
    match File::open(path) {
        Ok(mut f) => f.read_exact(&mut head).is_ok() && head == [0x1f, 0x8b],
        Err(_) => false,
    }
}

/// Merge `extra`'s locations into `treat` (port of `append_fwtrack`). Both are
/// pre-finalize; we just re-add every 5' position so the combined track behaves
/// like a single parse, then the caller finalizes once.
fn append_fwtrack(treat: &mut FwTrack, extra: FwTrack) {
    for chrom in extra.get_chr_names() {
        if let Some(loc) = extra.get_locations_by_chr(chrom) {
            for &p in &loc.plus {
                treat.add_loc(chrom, p, 0);
            }
            for &p in &loc.minus {
                treat.add_loc(chrom, p, 1);
            }
        }
    }
}

/// Best-effort R model script (port of `OutputWriter.model2r_script`). The
/// stderr `d` line is the contract; this file is informational only.
fn write_model_r(pm: &PeakModel, path: &str) -> std::io::Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "# R script for MACS3 model (macs3-rs)")?;
    writeln!(f, "p <- c({})", join_f32(&pm.plus_line))?;
    writeln!(f, "m <- c({})", join_f32(&pm.minus_line))?;
    writeln!(f, "ycorr <- c({})", join_f32(&pm.ycorr))?;
    writeln!(f, "xcorr <- c({})", join_f32(&pm.xcorr))?;
    writeln!(f, "altd <- c({})", join_i32(&pm.alternative_d))?;
    writeln!(f, "d <- {}", pm.d)?;
    Ok(())
}

fn join_f32(v: &[f32]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn join_i32(v: &[i32]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}
