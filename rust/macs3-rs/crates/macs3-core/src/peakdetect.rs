//! Peak-calling orchestration, ported from
//! `MACS3/Signal/PeakDetect.py::PeakDetect`.
//!
//! `PeakDetect` is the top-level coordinator behind `callpeak`: it decides
//! treatment/control scaling, builds the multi-scale control lambda parameters,
//! runs [`crate::caller::CallerFromAlignments`] and returns the peaks. SE only.

use std::io::Write;

use crate::caller::CallerFromAlignments;
use crate::peak_io::{BroadPeakIO, PeakIO};
use crate::track_fw::FwTrack;

/// Either a narrow or broad peak result, depending on `--broad`.
#[derive(Debug)]
pub enum DetectedPeaks {
    /// Narrow-peak result (default mode).
    Narrow(PeakIO),
    /// Broad-peak result (`--broad`).
    Broad(BroadPeakIO),
}

/// Inputs to a peak-calling run, gathered from the validated `options`.
pub struct PeakDetectParams {
    /// Effective genome size in bp (`gsize`).
    pub gsize: f64,
    /// Fragment length `d` (extsize or model d).
    pub d: i32,
    /// `--max-gap` resolved (tsize when unset).
    pub maxgap: i32,
    /// `--min-length` resolved (`d` when unset).
    pub minlen: i32,
    /// `--slocal` region size (`smalllocal`).
    pub sregion: i32,
    /// `--llocal` region size (`largelocal`).
    pub lregion: i32,
    /// `--shift` (`end_shift`).
    pub end_shift: i32,
    /// `--nolambda`.
    pub nolambda: bool,
    /// `--scale-to large` => true; default (`small`) => false (`tolarge`).
    pub tolarge: bool,
    /// `-p` given => Some(-log10 p); else None (use q).
    pub log_pvalue: Option<f64>,
    /// `-q` => -log10 q (used when `log_pvalue` is None).
    pub log_qvalue: f64,
    /// `--broad`.
    pub broad: bool,
    /// `--broad-cutoff` as -log10 (`log_broadcutoff`).
    pub log_broadcutoff: f64,
    /// `--call-summits`.
    pub call_summits: bool,
    /// pseudocount (always 1.0 for callpeak).
    pub pseudocount: f32,
    /// `-B`.
    pub store_bdg: bool,
    /// `--SPMR`.
    pub do_spmr: bool,
    /// `--trackline`.
    pub trackline: bool,
    /// bedGraph filename prefix (`name`).
    pub name: String,
}

/// Build the [`CallerFromAlignments`] following `PeakDetect.__call_peaks_*`,
/// including the scaling / lambda decisions. Returns the configured caller plus
/// the resolved `(min_length, max_gap)`.
pub fn build_caller(
    treat: &mut FwTrack,
    control: Option<&mut FwTrack>,
    p: &PeakDetectParams,
) -> CallerFromAlignments {
    let d = p.d;
    let gsize = p.gsize;

    if let Some(control) = control {
        // ---- with control (__call_peaks_w_control) ----
        let treat_sum = treat.total as i64 * d as i64;
        let control_sum = control.total as i64 * d as i64;
        let ratio_treat2control = treat_sum as f64 / control_sum as f64;

        // tocontrol decision (from callpeak_cmd.py): scale-to small is default.
        // tolarge: if t1 > c1 => scale control up (tocontrol=False); else scale
        //          treatment up (tocontrol=True).
        // small (default): if t1 > c1 => scale treatment down (tocontrol=True);
        //                  else scale control down (tocontrol=False).
        let t1 = treat.total;
        let c1 = control.total;
        let tocontrol = if p.tolarge {
            !(t1 > c1)
        } else {
            t1 > c1
        };

        let (lambda_bg, treat_scale): (f32, f32) = if tocontrol {
            (
                (control_sum as f64 / gsize) as f32,
                (1.0 / ratio_treat2control) as f32,
            )
        } else {
            ((treat_sum as f64 / gsize) as f32, 1.0)
        };

        // ctrl_d_s / ctrl_scaling_factor_s.
        let mut ctrl_d_s: Vec<i32> = vec![d];
        let mut ctrl_scale_s: Vec<f32> = Vec::new();

        // d-window scale.
        let tmp_v = if !tocontrol { ratio_treat2control } else { 1.0 };
        ctrl_scale_s.push(tmp_v as f32);

        // slocal.
        if p.sregion != 0 {
            ctrl_d_s.push(p.sregion);
            let tmp_v = if !tocontrol {
                d as f64 / p.sregion as f64 * ratio_treat2control
            } else {
                d as f64 / p.sregion as f64
            };
            ctrl_scale_s.push(tmp_v as f32);
        }
        // llocal.
        if p.lregion != 0 && p.lregion > p.sregion {
            ctrl_d_s.push(p.lregion);
            let tmp_v = if !tocontrol {
                d as f64 / p.lregion as f64 * ratio_treat2control
            } else {
                d as f64 / p.lregion as f64
            };
            ctrl_scale_s.push(tmp_v as f32);
        }

        if p.nolambda {
            ctrl_d_s.clear();
            ctrl_scale_s.clear();
        }

        CallerFromAlignments::new(
            treat,
            Some(control),
            d,
            ctrl_d_s,
            treat_scale,
            ctrl_scale_s,
            p.end_shift,
            lambda_bg,
            p.pseudocount,
            p.store_bdg,
            p.name.clone(),
            p.do_spmr,
        )
    } else {
        // ---- without control (__call_peaks_wo_control) ----
        let treat_total = treat.total as i64;
        let lambda_bg = (d as f64 * treat_total as f64 / gsize) as f32;
        let treat_scale: f32 = 1.0;

        let (ctrl_d_s, ctrl_scale_s): (Vec<i32>, Vec<f32>) = if !p.nolambda {
            (vec![p.lregion], vec![(d as f64 / p.lregion as f64) as f32])
        } else {
            (Vec::new(), Vec::new())
        };

        CallerFromAlignments::new(
            treat,
            None,
            d,
            ctrl_d_s,
            treat_scale,
            ctrl_scale_s,
            p.end_shift,
            lambda_bg,
            p.pseudocount,
            p.store_bdg,
            p.name.clone(),
            p.do_spmr,
        )
    }
}

/// Run the caller producing narrow or broad peaks, writing the two `-B`
/// bedGraphs to `treat_out`/`ctrl_out` when requested
/// (`PeakDetect.call_peaks`).
#[allow(clippy::too_many_arguments)]
pub fn run_caller<W: Write>(
    caller: &mut CallerFromAlignments,
    p: &PeakDetectParams,
    treat_out: Option<&mut W>,
    ctrl_out: Option<&mut W>,
) -> DetectedPeaks {
    if p.trackline {
        caller.enable_trackline();
    }

    // scoring symbol + cutoff.
    let (symbol, cutoff): (u8, f32) = match p.log_pvalue {
        Some(lp) => (b'p', lp as f32),
        None => (b'q', p.log_qvalue as f32),
    };

    if p.broad {
        // lvl1 cutoff = main p/q cutoff; lvl2 cutoff = broad-cutoff.
        let bp = caller.call_broadpeaks(
            &[symbol],
            &[cutoff],
            &[p.log_broadcutoff as f32],
            p.minlen,
            p.maxgap,
            p.maxgap * 4,
            treat_out,
            ctrl_out,
        );
        DetectedPeaks::Broad(bp)
    } else {
        let np = caller.call_peaks(
            &[symbol],
            &[cutoff],
            p.minlen,
            p.maxgap,
            p.call_summits,
            treat_out,
            ctrl_out,
        );
        DetectedPeaks::Narrow(np)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a track with `n_plus` plus-strand tags on `chrom` (positions spread
    /// out so dedup never removes any), finalized.
    fn track_with(chrom: &[u8], positions: &[i32]) -> FwTrack {
        let mut t = FwTrack::new(50, String::new(), 100000);
        for (i, &p) in positions.iter().enumerate() {
            // alternate strands but keep unique positions.
            t.add_loc(chrom, p, (i % 2) as i32);
        }
        t.finalize();
        t
    }

    fn base_params(gsize: f64, d: i32) -> PeakDetectParams {
        PeakDetectParams {
            gsize,
            d,
            maxgap: 50,
            minlen: d,
            sregion: 200,
            lregion: 600,
            end_shift: 0,
            nolambda: false,
            tolarge: false,
            log_pvalue: Some(5.0),
            log_qvalue: 0.0,
            broad: false,
            log_broadcutoff: 0.0,
            call_summits: false,
            pseudocount: 1.0,
            store_bdg: false,
            do_spmr: false,
            trackline: false,
            name: "test".to_string(),
        }
    }

    /// Golden from `test_PeakDetect.py::test_call_peaks_with_control_scales_*`
    /// (with `tocontrol=False`, i.e. treat <= control or `--scale-to large` with
    /// treat <= control). Here we use treat(50) < control(100) and default
    /// scale-to small, so the larger (control) is scaled down → tocontrol=False.
    /// d=150, gsize=10000, slocal=200, llocal=600.
    #[test]
    fn build_caller_with_control_scales() {
        let d = 150;
        let gsize = 10000.0;
        // treat.total = 50, control.total = 100.
        let mut treat = track_with(b"chr1", &(0..50).map(|i| i * 100).collect::<Vec<_>>());
        let mut control = track_with(b"chr1", &(0..100).map(|i| i * 50).collect::<Vec<_>>());
        assert_eq!(treat.total, 50);
        assert_eq!(control.total, 100);

        let mut p = base_params(gsize, d);
        p.sregion = 200;
        p.lregion = 600;

        let caller = build_caller(&mut treat, Some(&mut control), &p);
        let (ctrl_d_s, ctrl_scale_s, lambda_bg, treat_scale) = caller.scaling_params();

        // tocontrol=False: lambda_bg = treat_sum/gsize = 50*150/10000 = 0.75.
        assert!((lambda_bg - 0.75).abs() < 1e-5, "lambda_bg={lambda_bg}");
        assert!((treat_scale - 1.0).abs() < 1e-5);
        assert_eq!(ctrl_d_s, &[150, 200, 600]);
        // ratio = treat_sum/control_sum = 50/100 = 0.5.
        // scales = [0.5, 150/200*0.5=0.375, 150/600*0.5=0.125].
        let expected = [0.5_f32, 0.375, 0.125];
        for (o, e) in ctrl_scale_s.iter().zip(expected) {
            assert!((o - e).abs() < 1e-5, "scale {o} != {e}");
        }
    }

    /// Without control, the only background track is the llocal window scaled by
    /// `d/lregion` (`__call_peaks_wo_control`).
    #[test]
    fn build_caller_without_control() {
        let d = 150;
        let gsize = 10000.0;
        let mut treat = track_with(b"chr1", &(0..50).map(|i| i * 100).collect::<Vec<_>>());

        let p = base_params(gsize, d);
        let caller = build_caller(&mut treat, None, &p);
        let (ctrl_d_s, ctrl_scale_s, lambda_bg, treat_scale) = caller.scaling_params();

        assert_eq!(ctrl_d_s, &[600]);
        assert!((ctrl_scale_s[0] - (150.0_f32 / 600.0)).abs() < 1e-6);
        assert!((treat_scale - 1.0).abs() < 1e-6);
        // lambda_bg = d * treat.total / gsize = 150*50/10000 = 0.75.
        assert!((lambda_bg - 0.75).abs() < 1e-5);
    }
}
