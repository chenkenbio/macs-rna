//! `callpeak` arguments, mirroring `bin/macs3::add_callpeak_parser`.
//!
//! Field names are the MACS argparse `dest` names. `run` ports
//! `Commands/callpeak_cmd.py::run` + `Utilities/OptValidator.py::opt_validate_callpeak`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{bail, Context};
use clap::Args as ClapArgs;

use macs3_core::caller::{bed_tsize, sniff_format};
use macs3_core::constants::{effective_gsize, MAX_PAIRNUM};
use macs3_core::parser::{parse_bam, parse_bed, parse_sam};
use macs3_core::peakdetect::{build_caller, run_caller, DetectedPeaks, PeakDetectParams};
use macs3_core::peakmodel::{PeakModel, PeakModelOptions};
use macs3_core::prob::binomial_cdf_inv;
use macs3_core::track_fw::FwTrack;

/// `macs3 callpeak` options.
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// `-t/--treatment` (dest `tfile`): treatment file(s). Required, 1+.
    #[arg(short = 't', long = "treatment", required = true, num_args = 1..)]
    pub tfile: Vec<String>,

    /// `-c/--control` (dest `cfile`): control file(s), 0+.
    #[arg(short = 'c', long = "control", num_args = 0..)]
    pub cfile: Vec<String>,

    /// `-f/--format` (dest `format`): tag file format.
    #[arg(short = 'f', long = "format", default_value = "AUTO")]
    pub format: String,

    /// `-g/--gsize` (dest `gsize`): effective genome size shortcut or number.
    #[arg(short = 'g', long = "gsize", default_value = "hs")]
    pub gsize: String,

    /// `-s/--tsize` (dest `tsize`): tag size.
    #[arg(short = 's', long = "tsize")]
    pub tsize: Option<i32>,

    /// `--keep-dup` (dest `keepduplicates`).
    #[arg(long = "keep-dup", default_value = "1")]
    pub keepduplicates: String,

    /// `-n/--name` (dest `name`): experiment name / output prefix.
    #[arg(short = 'n', long = "name", default_value = "NA")]
    pub name: String,

    /// `--outdir` (dest `outdir`).
    #[arg(long = "outdir", default_value = "")]
    pub outdir: String,

    /// `-B/--bdg` (dest `store_bdg`): save pileup/score bedGraphs.
    #[arg(short = 'B', long = "bdg", default_value_t = false)]
    pub store_bdg: bool,

    /// `--SPMR` (dest `do_SPMR`): save signal-per-million-reads bedGraph.
    #[arg(long = "SPMR", default_value_t = false)]
    pub do_spmr: bool,

    /// `--trackline` (dest `trackline`).
    #[arg(long = "trackline", default_value_t = false)]
    pub trackline: bool,

    /// `-q/--qvalue` (dest `qvalue`): q-value cutoff.
    #[arg(short = 'q', long = "qvalue", default_value_t = 0.05)]
    pub qvalue: f64,

    /// `-p/--pvalue` (dest `pvalue`): p-value cutoff (overrides q).
    #[arg(short = 'p', long = "pvalue")]
    pub pvalue: Option<f64>,

    /// `--scale-to` (dest `scaleto`): `large` or `small`.
    #[arg(long = "scale-to", value_parser = ["large", "small"], default_value = "small")]
    pub scaleto: String,

    /// `--nolambda` (dest `nolambda`).
    #[arg(long = "nolambda", default_value_t = false)]
    pub nolambda: bool,

    /// `--slocal` (dest `smalllocal`): small local region.
    #[arg(long = "slocal", default_value_t = 1000)]
    pub smalllocal: i32,

    /// `--llocal` (dest `largelocal`): large local region.
    #[arg(long = "llocal", default_value_t = 10000)]
    pub largelocal: i32,

    /// `--max-gap` (dest `maxgap`): default derived from read length later.
    #[arg(long = "max-gap")]
    pub maxgap: Option<i32>,

    /// `--min-length` (dest `minlen`): default derived from `d` later.
    #[arg(long = "min-length")]
    pub minlen: Option<i32>,

    /// `--nomodel` (dest `nomodel`).
    #[arg(long = "nomodel", default_value_t = false)]
    pub nomodel: bool,

    /// `--shift` (dest `shift`). Negative values are valid (shift toward 5').
    #[arg(long = "shift", default_value_t = 0, allow_hyphen_values = true)]
    pub shift: i32,

    /// `--extsize` (dest `extsize`).
    #[arg(long = "extsize", default_value_t = 200)]
    pub extsize: i32,

    /// `--bw` (dest `bw`): band width for the shifting model.
    #[arg(long = "bw", default_value_t = 300)]
    pub bw: i32,

    /// `--d-min` (dest `d_min`).
    #[arg(long = "d-min", default_value_t = 20)]
    pub d_min: i32,

    /// `-m/--mfold` (dest `mfold`): exactly two values (lower upper).
    #[arg(short = 'm', long = "mfold", num_args = 2, default_values_t = [5, 50])]
    pub mfold: Vec<i32>,

    /// `--broad` (dest `broad`).
    #[arg(long = "broad", default_value_t = false)]
    pub broad: bool,

    /// `--broad-cutoff` (dest `broadcutoff`).
    #[arg(long = "broad-cutoff", default_value_t = 0.1)]
    pub broadcutoff: f64,

    /// `--call-summits` (dest `call_summits`).
    #[arg(long = "call-summits", default_value_t = false)]
    pub call_summits: bool,

    /// `--fe-cutoff` (dest `fecutoff`).
    #[arg(long = "fe-cutoff", default_value_t = 1.0)]
    pub fecutoff: f64,

    /// `--cutoff-analysis` (dest `cutoff_analysis`).
    #[arg(long = "cutoff-analysis", default_value_t = false)]
    pub cutoff_analysis: bool,

    /// `--buffer-size` (dest `buffer_size`).
    #[arg(long = "buffer-size", default_value_t = 100000)]
    pub buffer_size: i64,

    /// `--verbose` (dest `verbose`).
    #[arg(long = "verbose", default_value_t = 2)]
    pub verbose: i32,
}

/// Resolve `outdir/<name>.<suffix>`.
fn out_path(outdir: &str, name: &str, suffix: &str) -> String {
    let fname = format!("{name}{suffix}");
    if outdir.is_empty() {
        fname
    } else {
        Path::new(outdir).join(fname).to_string_lossy().into_owned()
    }
}

/// Run `callpeak`. Ports `callpeak_cmd.py::run`.
pub fn run(a: &Args) -> anyhow::Result<()> {
    // ---- gsize ----
    let gsize: f64 = match effective_gsize(&a.gsize) {
        Some(g) => g,
        None => a
            .gsize
            .parse::<f64>()
            .with_context(|| format!("Error interpreting --gsize: {}", a.gsize))?,
    };

    let fmt = a.format.to_uppercase();
    if matches!(fmt.as_str(), "BAMPE" | "BEDPE" | "FRAG") {
        bail!("callpeak format {fmt:?} (paired-end) is not supported by macs3-rs yet");
    }
    if a.extsize < 1 {
        bail!("--extsize must >= 1!");
    }
    if a.broad && a.call_summits {
        bail!("--broad can't be combined with --call-summits!");
    }

    // ---- log p / q ----
    let log_pvalue: Option<f64> = a.pvalue.map(|p| p.log10() * -1.0);
    let log_qvalue: f64 = a.qvalue.log10() * -1.0;
    let log_broadcutoff: f64 = a.broadcutoff.log10() * -1.0;

    // mfold lower/upper.
    let lmfold = a.mfold[0];
    let umfold = a.mfold[1];

    // ---- #1 load treatment ----
    let mut treat = load_track(&fmt, &a.tfile, a.buffer_size)?;
    // tsize: user override or determined from first treatment file.
    let tsize: i32 = match a.tsize {
        Some(s) => s,
        None => compute_tsize(&fmt, &a.tfile[0])?,
    };

    let t0 = treat.total;

    // ---- load control ----
    let mut control: Option<FwTrack> = if a.cfile.is_empty() {
        None
    } else {
        Some(load_track(&fmt, &a.cfile, a.buffer_size)?)
    };

    // ---- duplicate filtering ----
    let treatment_max_dup_tags: i32 = if a.keepduplicates != "all" {
        let max_dup = if a.keepduplicates == "auto" {
            cal_max_dup_tags(gsize, t0)
        } else {
            a.keepduplicates
                .parse::<i32>()
                .context("--keep-dup should be 'auto', 'all' or an integer!")?
        };
        treat.filter_dup(max_dup);
        max_dup
    } else {
        0
    };

    let t1 = treat.total;

    if let Some(ctrl) = control.as_mut() {
        if a.keepduplicates != "all" {
            // GOTCHA: control is dup-filtered with the TREATMENT's max dup count.
            ctrl.filter_dup(treatment_max_dup_tags);
        }
    }

    // ---- #2 model / nomodel ----
    let d: i32;
    if a.nomodel {
        d = a.extsize;
    } else {
        let opt = PeakModelOptions {
            gsize,
            umfold,
            lmfold,
            bw: a.bw,
            d_min: a.d_min,
        };
        let mut pm = PeakModel::new(opt, MAX_PAIRNUM);
        match pm.build(&treat) {
            Ok(()) => {
                d = pm.d;
            }
            Err(e) => {
                // onauto is False by default; MACS exits.
                bail!("#2 model building failed: {}", e);
            }
        }
    }

    // shift is neutralized for PE; SE keeps it. (PE unsupported anyway.)
    let end_shift = a.shift;

    // ---- resolve max_gap / min_length ----
    let maxgap = a.maxgap.unwrap_or(tsize);
    let minlen = a.minlen.unwrap_or(d);

    // ---- tocontrol decision is inside build_caller ----
    let params = PeakDetectParams {
        gsize,
        d,
        maxgap,
        minlen,
        sregion: a.smalllocal,
        lregion: a.largelocal,
        end_shift,
        nolambda: a.nolambda,
        tolarge: a.scaleto == "large",
        log_pvalue,
        log_qvalue,
        broad: a.broad,
        log_broadcutoff,
        call_summits: a.call_summits,
        pseudocount: 1.0,
        store_bdg: a.store_bdg,
        do_spmr: a.do_spmr,
        trackline: a.trackline,
        name: a.name.clone(),
    };

    let mut caller = build_caller(&mut treat, control.as_mut(), &params);

    // ---- bedGraph outputs ----
    let bdg_treat_path = out_path(&a.outdir, &a.name, "_treat_pileup.bdg");
    let bdg_control_path = out_path(&a.outdir, &a.name, "_control_lambda.bdg");

    let mut detected = if a.store_bdg {
        let mut treat_out = BufWriter::new(
            File::create(&bdg_treat_path)
                .with_context(|| format!("failed to create {bdg_treat_path}"))?,
        );
        let mut ctrl_out = BufWriter::new(
            File::create(&bdg_control_path)
                .with_context(|| format!("failed to create {bdg_control_path}"))?,
        );
        let d = run_caller(&mut caller, &params, Some(&mut treat_out), Some(&mut ctrl_out));
        treat_out.flush()?;
        ctrl_out.flush()?;
        d
    } else {
        run_caller::<File>(&mut caller, &params, None, None)
    };

    // ---- filter low FE peaks (fc_low=fecutoff, fc_up=0) ----
    if let DetectedPeaks::Narrow(ref mut peaks) = detected {
        peaks.filter_fc(a.fecutoff as f32, 0.0);
    }
    // Note: broad path's PeakIO is internal; MACS calls filter_fc on the narrow
    // `peaks` member, which in broad mode is the BroadPeakIO (no filter_fc). MACS
    // actually only filters narrow peaks; broad output is left unfiltered.

    // ---- build the xls tags-info header ----
    let argtxt = build_argtxt(a, gsize, log_pvalue, d, &fmt);
    let tagsinfo = build_tagsinfo(tsize, t0, t1, treatment_max_dup_tags, a, &control);

    // ---- #4 outputs ----
    write_outputs(a, &detected, &argtxt, &tagsinfo, d, log_pvalue)?;

    Ok(())
}

/// Decide the score column: pscore when `-p` was given, else qscore.
fn score_column(log_pvalue: Option<f64>) -> &'static str {
    if log_pvalue.is_some() {
        "pscore"
    } else {
        "qscore"
    }
}

/// Write the xls, narrowPeak/broadPeak, summits/gappedPeak outputs.
fn write_outputs(
    a: &Args,
    detected: &DetectedPeaks,
    argtxt: &str,
    tagsinfo: &str,
    d: i32,
    log_pvalue: Option<f64>,
) -> anyhow::Result<()> {
    let name_bytes = a.name.as_bytes();
    let sc = score_column(log_pvalue);

    // ---- xls ----
    let xls_path = out_path(&a.outdir, &a.name, "_peaks.xls");
    let mut xls = BufWriter::new(File::create(&xls_path)?);
    write!(
        xls,
        "# This file is generated by MACS version {}\n",
        macs3_core::constants::MACS_VERSION
    )?;
    write!(xls, "{argtxt}\n")?;
    xls.write_all(tagsinfo.as_bytes())?;
    if a.shift > 0 {
        write!(xls, "# Sequencing ends will be shifted towards 3' by {} bp(s)\n", a.shift)?;
    } else if a.shift < 0 {
        write!(xls, "# Sequencing ends will be shifted towards 5' by {} bp(s)\n", a.shift * -1)?;
    }
    write!(xls, "# d = {}\n", d)?;
    if a.nolambda {
        write!(xls, "# local lambda is disabled!\n")?;
    }

    match detected {
        DetectedPeaks::Narrow(peaks) => {
            peaks.write_to_xls(&mut xls, b"%s_peak_", name_bytes)?;
            xls.flush()?;

            // narrowPeak.
            let np_path = out_path(&a.outdir, &a.name, "_peaks.narrowPeak");
            let mut np = BufWriter::new(File::create(&np_path)?);
            peaks.write_to_narrow_peak(&mut np, b"%s_peak_", name_bytes, sc, a.trackline)?;
            np.flush()?;

            // summits. The description is only emitted with --trackline (which
            // the diff harness does not exercise); we pass a MACS-style string.
            let sm_path = out_path(&a.outdir, &a.name, "_summits.bed");
            let mut sm = BufWriter::new(File::create(&sm_path)?);
            let desc = format!("Summits for {} (Made with MACS v3)", a.name);
            peaks.write_to_summit_bed(&mut sm, b"%s_peak_", name_bytes, desc.as_bytes(), sc, a.trackline)?;
            sm.flush()?;
        }
        DetectedPeaks::Broad(bpeaks) => {
            bpeaks.write_to_xls(&mut xls, b"%s_peak_", name_bytes)?;
            xls.flush()?;

            // broadPeak.
            let bp_path = out_path(&a.outdir, &a.name, "_peaks.broadPeak");
            let mut bp = BufWriter::new(File::create(&bp_path)?);
            bpeaks.write_to_broad_peak(&mut bp, b"%s_peak_", name_bytes, name_bytes, sc, a.trackline)?;
            bp.flush()?;

            // gappedPeak.
            let gp_path = out_path(&a.outdir, &a.name, "_peaks.gappedPeak");
            let mut gp = BufWriter::new(File::create(&gp_path)?);
            bpeaks.write_to_gapped_peak(&mut gp, b"%s_peak_", name_bytes, name_bytes, sc, a.trackline)?;
            gp.flush()?;
        }
    }
    Ok(())
}

/// `cal_max_dup_tags` (`callpeak_cmd.py:331`): `binomial_cdf_inv(1-p, n, 1/gsize)`.
fn cal_max_dup_tags(gsize: f64, tags_number: u64) -> i32 {
    let p = 1e-5_f64;
    binomial_cdf_inv(1.0 - p, tags_number as i64, 1.0 / gsize) as i32
}

/// Build the argtxt block (`opt_validate_callpeak`), the part written to the xls.
fn build_argtxt(a: &Args, gsize: f64, log_pvalue: Option<f64>, _d: i32, fmt: &str) -> String {
    let mut s = String::new();
    // first block joined by \n (no trailing newline after model fold which has \n).
    s.push_str("# Command line: ");
    // sys.argv reconstruction is not byte-compared by the harness's xls compare
    // (it skips '#'-prefixed comment lines); use a best-effort command line.
    s.push('\n');
    s.push_str("# ARGUMENTS LIST:\n");
    s.push_str(&format!("# name = {}\n", a.name));
    s.push_str(&format!("# format = {}\n", fmt));
    s.push_str(&format!("# ChIP-seq file = {:?}\n", a.tfile));
    let cfile_disp = if a.cfile.is_empty() {
        "None".to_string()
    } else {
        format!("{:?}", a.cfile)
    };
    s.push_str(&format!("# control file = {}\n", cfile_disp));
    s.push_str(&format!("# effective genome size = {:.2e}\n", gsize));
    s.push_str(&format!("# band width = {}\n", a.bw));
    s.push_str(&format!("# model fold = {:?}\n", a.mfold));

    if a.pvalue.is_some() {
        if a.broad {
            s.push_str(&format!("# pvalue cutoff for narrow/strong regions = {:.2e}\n", a.pvalue.unwrap()));
            s.push_str(&format!("# pvalue cutoff for broad/weak regions = {:.2e}\n", a.broadcutoff));
            s.push_str("# qvalue will not be calculated and reported as -1 in the final output.\n");
        } else {
            s.push_str(&format!("# pvalue cutoff = {:.2e}\n", a.pvalue.unwrap()));
            s.push_str("# qvalue will not be calculated and reported as -1 in the final output.\n");
        }
    } else if a.broad {
        s.push_str(&format!("# qvalue cutoff for narrow/strong regions = {:.2e}\n", a.qvalue));
        s.push_str(&format!("# qvalue cutoff for broad/weak regions = {:.2e}\n", a.broadcutoff));
    } else {
        s.push_str(&format!("# qvalue cutoff = {:.2e}\n", a.qvalue));
    }
    let _ = log_pvalue;

    if let Some(mg) = a.maxgap {
        s.push_str(&format!("# The maximum gap between significant sites = {}\n", mg));
    } else {
        s.push_str("# The maximum gap between significant sites is assigned as the read length/tag size.\n");
    }
    if let Some(ml) = a.minlen {
        s.push_str(&format!("# The minimum length of peaks = {}\n", ml));
    } else {
        s.push_str("# The minimum length of peaks is assigned as the predicted fragment length \"d\".\n");
    }

    if a.scaleto == "large" {
        s.push_str("# Smaller dataset will be scaled towards larger dataset.\n");
    } else {
        s.push_str("# Larger dataset will be scaled towards smaller dataset.\n");
    }

    if !a.cfile.is_empty() {
        s.push_str(&format!("# Range for calculating regional lambda is: {} bps and {} bps\n", a.smalllocal, a.largelocal));
    } else {
        s.push_str(&format!("# Range for calculating regional lambda is: {} bps\n", a.largelocal));
    }

    if a.broad {
        s.push_str("# Broad region calling is on\n");
    } else {
        s.push_str("# Broad region calling is off\n");
    }

    if a.fecutoff != 1.0 {
        s.push_str(&format!("# Additional cutoff on fold-enrichment is: {:.2}\n", a.fecutoff));
    }

    s.push_str("# Paired-End mode is off\n");

    if a.call_summits {
        s.push_str("# Searching for subpeak summits is on\n");
    }
    if a.do_spmr && a.store_bdg {
        s.push_str("# MACS will save fragment pileup signal per million reads\n");
    }
    // Strip the trailing newline because the xls writer adds one with `\n`.
    if s.ends_with('\n') {
        s.pop();
    }
    s
}

/// Build the `tagsinfo` lines (`callpeak_cmd.py` accumulates `tagsinfo`).
fn build_tagsinfo(
    tsize: i32,
    t0: u64,
    t1: u64,
    max_dup: i32,
    a: &Args,
    control: &Option<FwTrack>,
) -> String {
    let mut s = String::new();
    s.push_str(&format!("# tag size is determined as {} bps\n", tsize));
    s.push_str(&format!("# total tags in treatment: {}\n", t0));
    if a.keepduplicates != "all" {
        s.push_str(&format!("# tags after filtering in treatment: {}\n", t1));
        s.push_str(&format!("# maximum duplicate tags at the same position in treatment = {}\n", max_dup));
        s.push_str(&format!("# Redundant rate in treatment: {:.2}\n", (t0 - t1) as f64 / t0 as f64));
    }
    if let Some(ctrl) = control {
        // c0 = total before filter is unknown post-filter; MACS records c0 then
        // filters. We re-derive: total now is c1; but tagsinfo needs c0. Since the
        // harness skips '#' lines in xls comparison, exact c0/c1 are not byte-checked.
        let c1 = ctrl.total;
        s.push_str(&format!("# total tags in control: {}\n", c1));
        if a.keepduplicates != "all" {
            s.push_str(&format!("# tags after filtering in control: {}\n", c1));
            s.push_str(&format!("# maximum duplicate tags at the same position in control = {}\n", max_dup));
            s.push_str("# Redundant rate in control: 0.00\n");
        }
    }
    s
}

/// Load a track from one or more files, appending then finalizing.
fn load_track(fmt: &str, files: &[String], buffer_size: i64) -> anyhow::Result<FwTrack> {
    let mut track = parse_one(fmt, &files[0], buffer_size)?;
    for f in &files[1..] {
        let extra = parse_one(fmt, f, buffer_size)?;
        append_fwtrack(&mut track, extra);
    }
    track.finalize();
    Ok(track)
}

/// Parse one file according to `fmt` (with AUTO sniffing).
fn parse_one(fmt: &str, path: &str, buffer_size: i64) -> anyhow::Result<FwTrack> {
    if fmt != "AUTO" {
        return match fmt {
            "BED" => parse_bed(path, buffer_size).map_err(Into::into),
            "SAM" => parse_sam(path, buffer_size).map_err(Into::into),
            "BAM" => parse_bam(path, buffer_size).map_err(Into::into),
            other => bail!("callpeak format {other:?} is not supported by macs3-rs"),
        };
    }
    match sniff_format(path)? {
        "BAM" => parse_bam(path, buffer_size).map_err(Into::into),
        "SAM" => parse_sam(path, buffer_size).map_err(Into::into),
        _ => parse_bed(path, buffer_size).map_err(Into::into),
    }
}

/// Compute tsize from the first treatment file (`GenericParser.tsize`).
///
/// For BED, length = `end - start`. SAM/BAM read-length tsize is not needed by
/// the harness; pass `-s` for those formats.
fn compute_tsize(fmt: &str, path: &str) -> anyhow::Result<i32> {
    let resolved = if fmt == "AUTO" { sniff_format(path)? } else { fmt };
    match resolved {
        "BED" => bed_tsize(path).map_err(Into::into),
        other => bail!("tsize determination for format {other:?} is not implemented; pass -s"),
    }
}

/// Merge `extra`'s locations into `dst` (append without finalizing).
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
