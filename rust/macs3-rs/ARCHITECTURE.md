# macs3-rs Architecture (Phase 0 foundation contract)

Bit-exact pure-Rust reimplementation of **MACS3 v3.0.4**. This document is the
**frozen contract** that the four parallel Phase-1 agents build on. Anything
marked FROZEN here must not change signature; fill in `todo!()` bodies in place.

Port source: `references/MACS/MACS3/` (Cython pure-python mode `.py` = spec).
Ground truth python: `macs3v304` conda env (`import MACS3`).

## Hard invariants
- All internal pileup / score values are `f32` (C `float`). Replicate MACS's
  `f64 -> f32` widening exactly. **No fast-math, no `ryu`.**
- Float output formatting goes through `fmt::*` (C `printf` `%.5f` / `%.6g`).
- Chromosome maps are `BTreeMap<Vec<u8>, _>` keyed by raw bytes so iteration is
  bytewise-sorted, matching MACS's `sorted(dict.keys())`.

## Per-file ownership (Phase 1)
| File | Owner | Status |
|---|---|---|
| `fmt.rs` | FROZEN (Phase 0) | done |
| `prob.rs` | FROZEN (Phase 0) | done |
| `constants.rs` | FROZEN (Phase 0) | done |
| `bedgraph_io.rs` | FROZEN (Phase 0) | done |
| `lib.rs`, `cli/main.rs`, `cli/cmd/mod.rs`, both `Cargo.toml` | FROZEN (Phase 0) | done |
| `pileup.rs`, `track_fw.rs`, `track_pe.rs`, `cmd/pileup.rs` | A | skeleton |
| `bedgraph.rs`, `peak_io.rs`, `cmd/{bdgopt,bdgpeakcall,bdgbroadcall}.rs` | B | skeleton |
| `scoretrack.rs`, `twocond.rs`, `cmd/{bdgcmp,bdgdiff}.rs` | C | skeleton |
| `peakmodel.rs`, `cmd/predictd.rs` | E | skeleton |
| `caller.rs`, `peakdetect.rs`, `cmd/callpeak.rs` | Phase 2 | skeleton |

Each owner fills `todo!()` bodies **in the same file**; they may add private
helpers and `pub` methods but must not change the frozen signatures below.

---

## Frozen public API

### `fmt` (`crates/macs3-core/src/fmt.rs`) — DONE
```rust
pub fn format_f5(v: f32) -> String;          // C printf("%.5f", (double)v)
pub fn write_f5(out: &mut String, v: f32);
pub fn format_g6(v: f32) -> String;          // C printf("%.6g", (double)v)
pub fn write_g6(out: &mut String, v: f32);
```

### `prob` (`crates/macs3-core/src/prob.rs`) — DONE (golden test passes)
```rust
pub const LSTEP: i32; pub const EXPTHRES: f64; pub const EXPSTEP: f64; pub const BIGX: f64;
pub fn logspace_add(logx: f64, logy: f64) -> f64;
pub fn factorial(n: u32) -> f64;
pub fn poisson_cdf(n: u32, lam: f64, lower: bool, log10: bool) -> f64;
pub fn poisson_cdf_inv(cdf: f64, lam: f64, maximum: i32) -> i32;
pub fn poisson_cdf_q_inv(cdf: f64, lam: f64, maximum: i32) -> i32;
pub fn chisq_pvalue_e(x: f64, df: u32) -> f64;
pub fn chisq_logp_e(x: f64, df: u32, log10: bool) -> f64;
pub fn binomial_pdf(x: i64, a: i64, b: f64) -> f64;
pub fn binomial_cdf(x: i64, a: i64, b: f64, lower: bool) -> f64;
pub fn binomial_sf(x: i64, a: i64, b: f64, lower: bool) -> f64;
pub fn binomial_cdf_inv(cdf: f64, a: i64, b: f64) -> i64;
```

### `constants` (`crates/macs3-core/src/constants.rs`) — DONE
`MACS_VERSION`, `MAX_PAIRNUM`, `MAX_LAMBDA`, `FESTEP`, `BUFFER_SIZE`,
`READ_BUFFER_SIZE`, `N_MP`, `EFFECTIVE_GS`, `effective_gsize(&str)->Option<f64>`,
and `DEFAULT_*` per-subcommand defaults. Extend (add) as needed downstream.

### `bedgraph_io` (`crates/macs3-core/src/bedgraph_io.rs`) — DONE (FROZEN)
```rust
pub fn read_bedgraph(path: &str) -> std::io::Result<BedGraphTrack>;
//   skips track/#/browse lines; atoi() for start/end; value parsed f64 then `as f32`;
//   transparent gzip via magic bytes; add_loc merges equal neighbours.
pub fn write_pv_exact<W: Write>(chrom: &[u8], pv: &[PosVal], out: &mut W, append: bool)
    -> std::io::Result<()>;
//   cPosValCalculation.c::write_pv_array_to_bedGraph: coalesce on EXACT `!=`,
//   pre_s=0, first interval [0, pv[0].pos), %.5f.
pub trait BreakPredicate { fn breaks(&self, pre: f32, cur: f32) -> bool; }
pub struct AlwaysBreak;   // BedGraphIO  (no coalesce)
pub struct Gt1e5;         // ScoreTrackII / CallPeakUnit:  |pre-cur| > 1e-5
pub struct Ge1e6;         // TwoConditionScores:           |pre-cur| >= 1e-6
pub fn write_bedgraph_predicate<W: Write, P: BreakPredicate>(
    chrom: &[u8], pos: &[i32], value: &[f32], pred: &P, out: &mut W) -> std::io::Result<()>;
pub fn bedgraph_trackline(name: &str, description: &str) -> String;
```

### `BedGraphTrack` (`bedgraph.rs`, owner B) — struct + add_loc FROZEN
```rust
pub struct ChromData { pub pos: Vec<i32>, pub val: Vec<f32> }
pub struct BedGraphTrack { /* private data */ pub maxvalue: f32, pub minvalue: f32, pub baseline_value: f32 }
impl BedGraphTrack {
  pub fn new(baseline_value: f32) -> Self;
  pub fn add_loc(&mut self, chromosome: &[u8], startpos: i32, endpos: i32, value: f32);
  pub fn add_loc_wo_merge(&mut self, chromosome: &[u8], startpos: i32, endpos: i32, value: f32);
  pub fn add_chrom_data(&mut self, chromosome: &[u8], pos: Vec<i32>, val: Vec<f32>);
  pub fn get_chr_names(&self) -> Vec<&[u8]>;                    // bytewise sorted
  pub fn get_data_by_chr(&self, chromosome: &[u8]) -> Option<&ChromData>;
  pub fn get_data_by_chr_mut(&mut self, chromosome: &[u8]) -> Option<&mut ChromData>;
  pub fn reset_baseline(&mut self, baseline_value: f32);
  pub fn total(&self) -> usize;
  pub fn is_empty(&self) -> bool;
  // STUB (Phase 1B): apply_func, p2q, call_peaks(-> PeakIO),
  //   call_broadpeaks(-> BroadPeakIO), cutoff_analysis
}
```

### `PosVal` + pileup fns (`pileup.rs`, owner A)
```rust
pub struct PosVal { pub pos: i32, pub val: f32 }   // FROZEN layout
impl PosVal { pub fn new(pos: i32, val: f32) -> Self; }
// STUB (Phase 1A): single_end_pileup, quick_pileup, max_over_two_pv_array
```

### Peak types (`peak_io.rs`, owner B) — struct + add FROZEN
```rust
pub struct NarrowPeak { chrom, start, end, length, summit, score, pileup, pscore, fc, qscore, name }
pub struct BroadPeak  { start, end, length, score, thick_start, thick_end, block_num,
                        block_sizes, block_starts, pileup, pscore, fc, qscore, name }
pub type GappedPeak = BroadPeak;   // gappedPeak is written from a BroadPeak (MACS has no distinct type)

pub struct PeakIO { /* private */ pub co_sorted: bool, pub total: i64, pub name: Vec<u8> }
impl PeakIO {
  pub fn new(name) -> Self; impl Default;     // default name b"MACS3"
  pub fn add(&mut self, chrom, start, end, summit, peak_score, pileup, pscore, fold_change, qscore, name);
  pub fn add_peak_content(&mut self, chrom, NarrowPeak);
  pub fn get_data_from_chrom(&mut self, chrom) -> &mut Vec<NarrowPeak>;
  pub fn peaks_by_chr(&self, chrom) -> Option<&Vec<NarrowPeak>>;
  pub fn get_chr_names(&self) -> Vec<&[u8]>;  // bytewise sorted
  pub fn sort(&mut self); pub fn iter(&self); pub fn is_empty(&self) -> bool;
  // STUB (Phase 1B): write_to_narrow_peak, write_to_xls, write_to_bed,
  //   write_to_summit_bed, filter_qscore, filter_fc
}
pub struct BroadPeakIO { /* private */ }
impl BroadPeakIO {
  pub fn new() -> Self; impl Default;
  pub fn add(&mut self, chrom, start, end, score, thick_start, thick_end, block_num,
             block_sizes, block_starts, pileup, pscore, fold_change, qscore, name);
  pub fn peaks_by_chr(&self, chrom); pub fn get_chr_names(&self); pub fn iter(&self);
  pub fn total(&self) -> i64; pub fn is_empty(&self) -> bool;
  // STUB (Phase 1B): write_to_gapped_peak, write_to_broad_peak, write_to_xls
}
```

### `ScoreTrack2` (`scoretrack.rs`, owner C)
```rust
pub struct ChromScore { pub pos: Vec<i32>, pub treat: Vec<f32>, pub ctrl: Vec<f32>, pub score: Vec<f32> }
pub struct ScoreTrack2 { /* private data */ pub treat_edm, ctrl_edm: f32,
   pub scoring_method, normalization_method: u8, pub pseudocount, cutoff: f32,
   pub trackline: bool, pub pvalue_stat: BTreeMap<u32,i64> }
impl ScoreTrack2 {
  pub fn new(treat_depth: f32, ctrl_depth: f32, pseudocount: f32) -> Self;
  pub fn set_pseudocount(&mut self, f32); pub fn enable_trackline(&mut self);
  pub fn add_chromosome(&mut self, chrom, chrom_max_len: usize);
  pub fn add(&mut self, chrom, endpos: i32, chip: f32, control: f32);
  pub fn finalize(&mut self);
  pub fn get_data_by_chr(&self, chrom) -> Option<&ChromScore>;
  pub fn get_data_by_chr_mut(&mut self, chrom) -> Option<&mut ChromScore>;
  pub fn get_chr_names(&self) -> Vec<&[u8]>;  pub fn total(&self) -> i64;
  // STUB (Phase 1C/2): change_score_method(u8), change_normalization_method(u8),
  //   make_pq_table()->BTreeMap<u32,f32>, call_peaks(->PeakIO),
  //   call_broadpeaks(->BroadPeakIO), cutoff_analysis(->String),
  //   write_bedgraph (uses bedgraph_io::Gt1e5)
}
```

### `TwoConditionScores` (`twocond.rs`, owner C)
```rust
pub struct ChromTwoCond { pub pos: Vec<i32>, pub llr_1, llr_2, llr_3: Vec<f32> }
pub struct TwoConditionScores { /* private data */ pub cond1_factor, cond2_factor, pseudocount, cutoff: f32,
   pub t1bdg, c1bdg, t2bdg, c2bdg: BedGraphTrack }
impl TwoConditionScores {
  pub fn new(t1bdg, c1bdg, t2bdg, c2bdg: BedGraphTrack, cond1_factor, cond2_factor, pseudocount: f32) -> Self;
  pub fn set_pseudocount(&mut self, f32);
  pub fn get_data_by_chr(&self, chrom) -> Option<&ChromTwoCond>;
  pub fn get_chr_names(&self) -> Vec<&[u8]>;
  pub fn push_interval(&mut self, chrom, endpos: i32, llr1, llr2, llr3: f32);
  // STUB (Phase 1C): build(), call_peaks(->(PeakIO,PeakIO,PeakIO)),
  //   write_bedgraph (uses bedgraph_io::Ge1e6)
}
```

### `FwTrack` (`track_fw.rs`, owner A) — simple paths implemented
```rust
pub const INT_MAX: i32 = i32::MAX;
pub struct ChromLoc { pub plus: Vec<i32>, pub minus: Vec<i32> }
pub struct FwTrack { /* private */ pub is_sorted: bool, pub total: u64, pub fw: i32,
   pub length: u64, pub buffer_size: i64, pub annotation: String }
impl FwTrack {
  pub fn new(fw: i32, anno: String, buffer_size: i64) -> Self;
  pub fn add_loc(&mut self, chrom, fiveendpos: i32, strand: i32);  // strand 0=+,1=-
  pub fn finalize(&mut self);  pub fn sort(&mut self);
  pub fn set_rlengths(&mut self, &BTreeMap<Vec<u8>,i32>);
  pub fn get_rlengths(&mut self) -> &BTreeMap<Vec<u8>,i32>;
  pub fn get_locations_by_chr(&self, chrom) -> Option<&ChromLoc>;
  pub fn get_chr_names(&self) -> Vec<&[u8]>;
  pub fn filter_dup(&mut self, maxnum: i32) -> u64;   // maxnum<0 short-circuit done; scan STUB (Phase 1A)
}
```

### `PeTrack` (`track_pe.rs`, owner A) — minimal skeleton (PE unused by macs-rna)
```rust
pub struct ChromPe { pub starts: Vec<i32>, pub ends: Vec<i32> }
pub struct PeTrack { /* private */ pub is_sorted: bool, pub total, length: u64,
   pub average_template_length: f32, pub buffer_size: i64, pub annotation: String }
impl PeTrack {
  pub fn new(anno: String, buffer_size: i64) -> Self;
  pub fn add_loc(&mut self, chrom, start: i32, end: i32);
  pub fn get_locations_by_chr(&self, chrom); pub fn get_chr_names(&self);
  // STUB (Phase 2): finalize, filter_dup
}
```

### `PeakModel` (`peakmodel.rs`, owner E)
```rust
pub struct NotEnoughPairs { pub value: String }   // impl Error/Display
pub struct PeakModelOptions { pub gsize: f64, pub umfold, lmfold, bw, d_min: i32 }
pub struct PeakModel { pub gz: f64, pub max_pairnum, umfold, lmfold, bw, d_min,
   tag_expansion_size, min_tags, max_tags, peaksize, d, scan_window: i32,
   pub alternative_d: Vec<i32>, pub plus_line, minus_line, shifted_line, xcorr, ycorr: Vec<f32>,
   pub summary: String }
impl PeakModel {
  pub fn new(opt: PeakModelOptions, max_pairnum: i32) -> Self;
  pub fn build(&mut self) -> Result<(), NotEnoughPairs>;   // STUB (Phase 1E)
}
```

### `CallerFromAlignments` (`caller.rs`) / `PeakDetect` (`peakdetect.rs`) — Phase 2
Opaque skeletons fixing the public names + `call_peaks` / `call_broadpeaks`
entry points (`todo!()`). Constructor shapes finalised in Phase 2 with the
concrete track types they consume.

---

## CLI surface (`crates/macs3-cli`)
`main.rs` (FROZEN) defines clap `Cli` (`name="macs3-rs"`) with 8 subcommands
dispatching to `cmd::<name>::run(&args)`. `cmd/mod.rs` (FROZEN) declares them.
Each `cmd/<name>.rs` has `#[derive(clap::Args)] pub struct Args` whose fields
use MACS `dest` names, and `pub fn run(&Args) -> anyhow::Result<()>` = `todo!()`.

Flags verified against `bin/macs3` argparse:
- `callpeak`: `-t -c -f -g -s --keep-dup -n --outdir -B --SPMR --trackline -q -p
  --scale-to --nolambda --slocal --llocal --max-gap --min-length --nomodel
  --shift --extsize --bw --d-min -m/--mfold --broad --broad-cutoff --call-summits
  --fe-cutoff --cutoff-analysis --buffer-size --verbose`
- `predictd`: `-i -f -g -s --bw --d-min -m --rfile --outdir --buffer-size --verbose`
- `pileup`: `-i -o --outdir -f -B/--both-direction --extsize --buffer-size --verbose`
  (note: 3.0.4 pileup has `-B` = both-direction and no `--shift`/`--scale-to`)
- `bdgcmp`: `-t -c -S -p -m(1+) --o-prefix -o(1+) --outdir --verbose`
- `bdgopt`: `-i -m -p(0+) -o --outdir --verbose`
- `bdgpeakcall`: `-i -c -l -g --call-summits --cutoff-analysis --no-trackline
  --o-prefix -o --outdir --verbose`  (`trackline()` accessor for store_false)
- `bdgbroadcall`: `-i -c -C -l -g -G --no-trackline --o-prefix -o --outdir --verbose`
- `bdgdiff`: `--t1 --c1 --t2 --c2 -C -l -g --d1 --d2 --o-prefix -o(=3) --outdir --verbose`

## Gate status (Phase 0 exit)
- `cargo build` (workspace): OK, 0 warnings.
- `cargo test -p macs3-core`: 31 passed (fmt golden, prob golden, bedgraph_io
  incl. `write_pv_exact` cross-checked against `cPosValCalculation.c::main`).

---

## Phase 1A addendum — parser API (appended; ports `IO/Parser.py`)

Single-end parsers build an `FwTrack` of 5' ends (strand `0=+`, `1=-`). None
call `finalize()`; the caller finalizes after appending all inputs (mirrors
`load_tag_files_options`). `parse_bam` additionally calls `FwTrack::set_rlengths`
with the BAM-header lengths, as `BAMParser.build_fwtrack` does.

```rust
// crates/macs3-core/src/parser.rs
pub fn parse_bed(path: &str, buffer_size: i64) -> std::io::Result<FwTrack>;
pub fn parse_sam(path: &str, buffer_size: i64) -> std::io::Result<FwTrack>;
pub fn parse_bam(path: &str, buffer_size: i64) -> std::io::Result<FwTrack>;
pub fn parse_bedpe(path: &str, buffer_size: i64) -> std::io::Result<PeTrack>; // STUB
pub fn parse_bampe(path: &str, buffer_size: i64) -> std::io::Result<PeTrack>; // STUB
```

5'-end / filter semantics (bit-exact with `Parser.py`):
- **BED** (`BEDParser`): skip leading `track`/`browser`/`#` lines. `+` -> col1
  (leftmost, strand 0); `-` -> col2 (rightmost, strand 1); missing strand col ->
  default `+`. Lines with `pos < 0` or empty chrom are dropped.
- **SAM** (`SAMParser`): skip `@` header lines. Filter
  `(flag & 2820) || (flag & 1 && (flag & 136 || !(flag & 2)))`. `+`: `pos-1`;
  `-` (flag&16): `pos-1 + Σ len(M/D/N/X/=)` from CIGAR. `.fa` suffix stripped
  from refname. (Note: macs3 v3.0.4's own SAM path raises `TypeError` — a stock
  bug — so SAM can't be byte-diffed against it; logic verified by unit tests.)
- **BAM** (`BAMParser`): BGZF-decompressed via `noodles-bgzf`, then raw-byte
  parsed at MACS's exact offsets (`bam_fw_binary_parse`): `n_cigar_op,flag` @
  bytes 12..16, `refID,pos` @ 0..8, CIGAR @ `32 + l_read_name`. Same flag mask;
  minus strand advances `pos` by reference-consuming CIGAR ops (codes 0/2/3/7/8).
  `set_rlengths(header_lengths)` is applied.

### `FwTrack` additions (Phase 1A, non-breaking — added methods/fields only)
- `pub fn rlengths_iter_order(&self) -> Vec<&[u8]>` — chromosome order that
  `get_rlengths()` / `pileup_and_write_se` iterate (`list(chrlengths.keys())`):
  first-seen insertion order when `set_rlengths` was never called (BED/SAM),
  else `sorted(valid) ++ sorted(missed)` (BAM). `pileup` uses this for the
  "emit chromosomes in rlengths insertion order" gotcha.
- `filter_dup(maxnum>=0)` now implements the de-dup scan, replicating the quirk
  that strands with `<= 1` element are left untouched and **not** counted into
  `total` (matches `FixWidthTrack.py::filter_dup`).

## Phase 1A gate results
- `cargo build`: 0 errors / 0 warnings in owned files
  (`parser.rs`, `pileup.rs`, `track_fw.rs`, `track_pe.rs`, `cmd/pileup.rs`).
- `cargo test -p macs3-core` (parser:: pileup:: track_fw::): all pass —
  golden values ported from `test_Pileup.py` (SE/quick/max) and
  `test_FixWidthTrack.py` (`filter_dup` 17->16->14->11).
- `diff_macs.py --only pileup`: **PASS** for `pileup_chip` and `pileup_ctrl`,
  byte-identical, worst delta `0.0`.
- Bonus e2e: `pileup -f BAM` on `yeast_500k_SRR1822137.bam` is byte-identical to
  macs3 (475 333 lines), validating BAM parse + `set_rlengths` ordering.
