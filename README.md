# CRAMboozle

**CRAMboozle** is a tool and Snakemake workflow for de-identifying sequencing data stored in BAM or CRAM format, protecting the genetic privacy of donor individuals while preserving essential alignment data for analysis.

CRAMboozle is a derivative of [**BAMboozle**](https://github.com/sandberg-lab/dataprivacy) by Christoph Ziegenhain and Rickard Sandberg (Karolinska Institutet). The original tool is described in:

> Ziegenhain, C., Sandberg, R. BAMboozle removes genetic variation from human sequence data for open data sharing. *Nat Commun* **12**, 6216 (2021). [https://doi.org/10.1038/s41467-021-26152-8](https://doi.org/10.1038/s41467-021-26152-8)

## What's different from BAMboozle?

CRAMboozle extends the original BAMboozle (v0.5.0) with the following changes:

| Feature | BAMboozle | CRAMboozle |
|---|---|---|
| **Input formats** | BAM only | BAM or CRAM (auto-detected) |
| **Output format** | BAM | CRAM v3.1 by default (configurable) |
| **Compression** | Default BAM | CRAM v3.1 with level 0–9 and advanced codecs (LZMA, BZIP2, FQZ, TOK, ARITH) |
| **CPU allocation** | Manual `--p` flag | Auto-detects available CPUs (still overridable) |
| **Batch processing** | Single file at a time | Snakemake workflow for parallel multi-sample processing on SLURM clusters |
| **Reference validation** | None | Validates CRAM–reference compatibility before processing |

The core de-identification logic (sequence replacement, tag sanitization, splice handling) is preserved from the original BAMboozle.

## How it works

CRAMboozle replaces observed read sequences with reference genome sequence and sanitizes auxiliary tags, removing donor-specific genetic variation. For full details on the de-identification strategy (SNP replacement, insertion/deletion handling, clipping, splicing, tag sanitization, strict mode), see the [BAMboozle paper](https://doi.org/10.1038/s41467-021-26152-8) and [original README](https://github.com/sandberg-lab/dataprivacy#readme).

## Requirements

- Python 3.10+
- [pysam](https://pysam.readthedocs.io/) 0.21.0+ (for CRAM v3.1 support)
- [Snakemake](https://snakemake.readthedocs.io/) (for batch workflow; optional for single-file use)
- Reference genome FASTA file, indexed with `samtools faidx`

## Quick start (single file)

```bash
python CRAMboozle.py \
    --input sample.bam \
    --out sample_deidentified.cram \
    --fa /path/to/reference.fasta
```

Key options:

| Flag | Description |
|---|---|
| `--input` | Input BAM or CRAM file |
| `--out` | Output file path (CRAM by default) |
| `--fa` | Reference genome FASTA (indexed) |
| `--p N` | Number of processes (default: all available CPUs) |
| `--strict` | Also sanitize mapping quality & extra auxiliary tags |
| `--keepunmapped` | Keep unmapped reads in output |
| `--keepsecondary` | Keep secondary alignments in output |
| `--force-bam` | Force BAM output instead of CRAM |
| `--compression-level` | CRAM compression level 0–9 (default: 9) |

## Snakemake workflow (batch processing)

The included Snakemake workflow processes many samples in parallel on a SLURM cluster.

### 1. Configure your samples

Edit `config/samples.yaml`:

```yaml
samples:
  sample1: /path/to/sample1.bam
  sample2: /path/to/sample2.cram
  sample3: /path/to/sample3.bam
```

> **Tip:** See `config/generate_samples.txt` for a one-liner to auto-generate this file from a directory of CRAMs.

### 2. Configure settings

Edit `config/config.yaml` to set:
- `reference_genome` — path to your indexed FASTA reference
- `results_dir` — output directory (default: `results/`)
- `strict_mode`, `keep_unmapped`, `keep_secondary` — processing options

Edit `config/cluster_slurm.yaml` to match your cluster's partition names, memory, and CPU limits.

### 3. Run

```bash
# Dry run first (recommended)
snakemake -s CRAMboozle.snakefile \
    --cluster-config config/cluster_slurm.yaml \
    --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" \
    -j 40 --latency-wait 60 --keep-going -np

# Remove -np to execute
```

For local execution (small datasets, no cluster):
```bash
snakemake -s CRAMboozle.snakefile -j 4
```

## Output

For each sample, the workflow produces:
- `{sample}_deidentified.cram` — de-identified CRAM file
- `{sample}_deidentified.cram.crai` — CRAM index
- `logs/{sample}_cramboozle.log` — processing log

## Repository structure

```
CRAMboozle/
├── CRAMboozle.py          # Core de-identification script
├── CRAMboozle.snakefile   # Snakemake workflow
├── config/
│   ├── config.yaml        # Pipeline settings
│   ├── samples.yaml       # Sample manifest
│   └── cluster_slurm.yaml # SLURM cluster resources
└── results/               # Output directory
```

## Troubleshooting

- **Reference mismatch**: CRAMboozle validates the reference before processing CRAM files. If validation fails, check your CRAM headers (`samtools view -H file.cram | grep '^@SQ'`) and ensure the correct reference build is specified.
- **Missing index**: The reference FASTA must be indexed (`samtools faidx reference.fasta`). Input files should be coordinate-sorted and indexed, though CRAMboozle will attempt to sort and index them if not.
- **Memory issues**: Increase `mem` in `config/cluster_slurm.yaml`.

## License

This project is licensed under the **GNU General Public License v3.0** — the same license as the original BAMboozle. See [LICENSE](LICENSE) for details.

## Attribution

CRAMboozle is derived from [BAMboozle](https://github.com/sandberg-lab/dataprivacy) by Christoph Ziegenhain and Rickard Sandberg, originally released under GPL-3.0. If you use CRAMboozle, please cite the original BAMboozle paper:

> Ziegenhain, C., Sandberg, R. BAMboozle removes genetic variation from human sequence data for open data sharing. *Nat Commun* **12**, 6216 (2021). [https://doi.org/10.1038/s41467-021-26152-8](https://doi.org/10.1038/s41467-021-26152-8)

CRAMboozle modifications by Robert Patton ([Ha Lab](https://gavinhalab.org/), Fred Hutch Cancer Center).