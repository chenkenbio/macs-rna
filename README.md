# macs-rna

Strand-specific MACS3 peak calling for RNA-seq data (e.g., m6A-MeRIP-seq, eCLIP).

## Problem

Running `macs3 callpeak` independently on each strand after splitting a BAM introduces two issues:

1. **Unequal scaling**: Strands with different read depths get different statistical power, leading to inconsistent peak detection sensitivity.
2. **Incorrect pileup from PE data**: MACS3 with `-f BAM` treats every read (R1+R2) as an independent tag. R2 is offset from the actual signal site, blurring peak resolution. For spliced RNA-seq, `-f BAMPE` is even worse — it uses fragment coordinates spanning introns.

## Solution

`macs-rna` wraps MACS3 with RNA-seq-aware defaults:

- **R1-only by default** (`--read 1`): Only R1 reads enter peak calling, giving sharper peaks. R2 can be selected for eCLIP (`--read 2`).
- **Global scaling**: Uses the total reads across both strands as the effective library size, equalizing sensitivity.
- **No BAMPE**: Excluded from format choices to prevent incorrect spliced-read pileup.

The step-by-step pipeline:

1. Split BAM by strand (+ filter to R1/R2) via samtools flag filters
2. `macs3 callpeak -B` per strand — generates bedGraph pileup and lambda tracks
3. `macs3 bdgopt -m multiply -p {total/strand}` — scales bedGraphs to total library size
4. `macs3 bdgcmp -m ppois` — computes Poisson p-value scores
5. `macs3 bdgopt -m p2q` — BH-corrected q-values (when using `-q` cutoff)
6. `macs3 bdgpeakcall` / `bdgbroadcall` — final peak calling with equalized sensitivity

When `--SPMR` is passed, the scaling is done via `bdgcmp -S {total_reads_in_millions}` on SPMR-normalized bedGraphs instead.

## Installation

```bash
pip install -e .
```

Requires `macs3>=3.0.0`, `samtools`, and optionally `bedGraphToBigWig` (for bigwig conversion).

## Usage

```bash
# m6A-MeRIP-seq (R1-only, default)
macs-rna callpeak \
  -t ip.bam -c input.bam \
  --libtype FR \
  -f BAM -g hs --keep-dup all \
  --nomodel --extsize 150 \
  -q 0.05 \
  -n experiment --outdir results/ \
  --chrom-sizes hg38

# eCLIP (R2 at crosslink site)
macs-rna callpeak \
  -t eclip_ip.bam -c eclip_input.bam \
  --libtype FR --read 2 \
  -f BAM -g hs --keep-dup 1 \
  --nomodel --extsize 150 \
  -q 0.05 \
  -n eclip_exp --outdir results/
```

### Key arguments

| Argument | Description |
|----------|-------------|
| `--libtype {FR,RF,F,R}` | **(required)** Library strandedness type |
| `--read {1,2,both}` | Which read to use for PE data (default: `1`). Use `2` for eCLIP |
| `--chrom-sizes` | Chromsizes file or shortcut (`hg38`, `mm10`) for bdg-to-bigwig conversion |
| `--min-mapq` | Minimum mapping quality for strand split (default: 0) |
| `--primary` | Only keep primary alignments |
| `--mode {m6A-MeRIP,seCLIP}` | Preset recipe for common experiment types |
| `--dry-run` | Print all commands without executing |

Most `macs3 callpeak` arguments (`-t`, `-c`, `-f`, `-g`, `--keep-dup`, `--nomodel`, `--extsize`, `-q`/`-p`, `--broad`, etc.) are supported. Notable differences from macs3:
- `-f` only accepts `AUTO`, `BAM`, `SAM` (no BAMPE/BED/BEDPE — strand splitting requires BAM/SAM)
- `-t` and `-c` accept a single file each (no multi-file pooling)
- `--call-summits` is not available (incompatible with the step-by-step rescaling pipeline)
- `--primary` excludes both secondary (0x100) and supplementary (0x800) alignments

### Dry-run

Preview the full command sequence without running anything:

```bash
macs-rna callpeak \
  -t ip.bam --libtype FR -f BAM \
  --keep-dup all --nomodel --extsize 150 \
  -q 0.05 -n test --outdir out/ --dry-run
```

## Output

```
{outdir}/
├── {name}_peaks.narrowPeak              # Combined peaks (both strands, strand column set)
├── {name}_fwd_peaks.narrowPeak          # Forward strand peaks
├── {name}_rev_peaks.narrowPeak          # Reverse strand peaks
├── {name}_fwd_treat_pileup.bw           # Forward signal (if --chrom-sizes + -B)
├── {name}_rev_treat_pileup.bw           # Reverse signal
├── {name}_fwd_control_lambda.bw         # Forward control lambda
└── {name}_rev_control_lambda.bw         # Reverse control lambda
```

- When `--SPMR` is used, bedGraph/bigWig files are labeled `.SPMR.bw` / `.SPMR.bdg`
- When `--chrom-sizes` is provided, bdg files are converted to bigwig and originals are removed

## How scaling works

By default, both the treatment pileup and control lambda bedGraphs are multiplied by `total_treatment_reads / strand_treatment_reads` using `bdgopt -m multiply`. This inflates the per-strand counts to match the total library size, giving the Poisson test equal statistical power regardless of per-strand read distribution.

With `--SPMR`, SPMR normalizes pileup to "per million strand reads," and `bdgcmp -S {total_reads_M}` converts back to effective counts at the total library size.

## License

MIT
