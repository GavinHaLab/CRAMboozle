# CRAMboozle Snakemake Workflow

A Snakemake workflow for de-identifying BAM and CRAM files using CRAMboozle.py.

## Overview

This workflow processes multiple BAM/CRAM files in parallel, de-identifying sequencing data by removing identifying information while preserving the essential alignment data for analysis.

## Requirements

- Python 3.7+
- Snakemake
- pysam
- samtools
- Reference genome FASTA file (indexed with `samtools faidx`)

## Setup

### 1. Configure your samples

Edit `config/samples.yaml` to specify your samples and their input file paths:

```yaml
samples:
  patient001: /path/to/patient001.bam
  patient002: /path/to/patient002.cram
  control001: /path/to/control001.bam
```

Edit `config/config.yaml` to specify:
- Reference genome location
- Output directory
- Processing options

### 2. Configure cluster settings (Fred Hutch)

The `config/cluster_slurm.yaml` is pre-configured for Fred Hutch SLURM cluster.
Modify if using a different cluster system.

## Running the Workflow

### Local execution (small datasets)
```bash
snakemake -s CRAMboozle.snakefile -j 4
```

### Cluster execution (Fred Hutch)
```bash
# Load required modules
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml Python/3.7.4-foss-2019b-fh1

# Run workflow
snakemake -s CRAMboozle.snakefile \
    --latency-wait 60 \
    --keep-going \
    --cluster-config config/cluster_slurm.yaml \
    --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" \
    -j 40
```

### Dry run (recommended first)
Add `-np` flag to the end of any command to see what would be executed without running it.

## Output Files

For each sample, the workflow generates:
- `{sample}_deidentified.cram` - De-identified CRAM file
- `{sample}_deidentified.cram.crai` - CRAM index file
- `logs/{sample}_cramboozle.log` - Processing log
- `CRAMboozle_summary.txt` - Summary of all processed files

## Configuration Options

In `config/config.yaml`:

- `strict_mode`: Enable additional tag sanitization (default: false)
- `keep_unmapped`: Keep unmapped reads in output (default: false)  
- `keep_secondary`: Keep secondary alignments (default: false)
- `cramboozle.ncpus`: CPUs per job (default: 8, auto-detected if available)

## File Formats

- **Input**: BAM or CRAM files (auto-detected by extension)
- **Output**: CRAM files by default (more compressed than BAM)

## Monitoring

- Check cluster job status: `squeue -u $USER`
- View logs in: `results/logs/`
- Monitor progress: `snakemake -s CRAMboozle.snakefile --summary`

## Troubleshooting

1. **Missing reference**: Ensure reference genome FASTA is indexed
   ```bash
   samtools faidx /path/to/reference.fasta
   ```

2. **Permission errors**: Check file permissions and paths

3. **Memory issues**: Increase memory in `cluster_slurm.yaml`

4. **Failed jobs**: Check individual log files in `results/logs/`

## Example Directory Structure

```
CRAMboozle/
├── CRAMboozle.py
├── CRAMboozle.snakefile
├── config/
│   ├── config.yaml
│   └── cluster_slurm.yaml
└── results/
    ├── sample1_deidentified.cram
    ├── sample1_deidentified.cram.crai
    ├── logs/
    └── CRAMboozle_summary.txt
```