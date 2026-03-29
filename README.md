# macs-rna

Strand-specific MACS3 peak calling for RNA-seq data (e.g., m6A-MeRIP-seq, eCLIP).

## Problem

Running `macs3 callpeak` independently on each strand after splitting a BAM introduces a scaling inconsistency: strands with different read depths get different statistical power, leading to unequal sensitivity in peak detection.

## Solution

`macs-rna` splits BAMs by strand but uses the **total reads across both strands** as the effective library size for the Poisson test. This is achieved via a step-by-step MACS3 pipeline:

1. `macs3 callpeak --SPMR -B` per strand — generates SPMR-normalized bedGraphs
2. `macs3 bdgcmp -S {total_reads_in_millions}` — re-scores with global scaling
3. `macs3 bdgopt -m p2q` — BH-corrected q-values (when using `-q` cutoff)
4. `macs3 bdgpeakcall` / `bdgbroadcall` — final peak calling with equalized sensitivity

## Installation

```bash
pip install -e .
```

Requires `macs3>=3.0.0`, `samtools`, and optionally `bedGraphToBigWig` (for bigwig conversion).

## Usage

```bash
macs-rna callpeak \
  -t ip.bam -c input.bam \
  --libtype FR \
  -f BAM -g hs --keep-dup all \
  --nomodel --extsize 150 \
  -q 0.05 \
  -n experiment --outdir results/ \
  --chrom-sizes hg38
```

### Key arguments

| Argument | Description |
|----------|-------------|
| `--libtype {FR,RF,F,R}` | **(required)** Library strandedness type |
| `--chrom-sizes` | Chromsizes file or shortcut (`hg38`, `mm10`) for bdg-to-bigwig conversion |
| `--min-mapq` | Minimum mapping quality for strand split (default: 0) |
| `--primary` | Only keep primary alignments |
| `--mode {m6A-MeRIP,seCLIP}` | Preset recipe for common experiment types |
| `--dry-run` | Print all commands without executing |

All standard `macs3 callpeak` arguments (`-t`, `-c`, `-f`, `-g`, `--keep-dup`, `--nomodel`, `--extsize`, `-q`/`-p`, `--broad`, etc.) are supported.

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
├── {name}_fwd_treat_pileup.SPMR.bw     # Forward signal (if --chrom-sizes + -B/--SPMR)
├── {name}_rev_treat_pileup.SPMR.bw     # Reverse signal
├── {name}_fwd_control_lambda.SPMR.bw   # Forward control lambda
└── {name}_rev_control_lambda.SPMR.bw   # Reverse control lambda
```

- BedGraph files are labeled `.SPMR.bdg` to explicitly indicate normalization type
- When `--chrom-sizes` is provided, bdg files are converted to bigwig and originals are removed

## How scaling works

SPMR normalizes pileup to "per million strand reads." When both strands are re-scaled via `bdgcmp -S {total_reads_M}`, the effective library size becomes identical for both strands, giving the Poisson test equal statistical power regardless of per-strand read distribution.

## License

MIT
