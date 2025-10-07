# CRAMboozle.snakefile
# Robert Patton, rpatton@fredhutch.org (Ha Lab)
# v1.1, 10/07/2025

"""
CRAMboozle Snakemake workflow for de-identifying BAM/CRAM files

# before running snakemake at Fred Hutch, do in tmux terminal:
ml snakemake/7.32.3-foss-2022b
ml Python/3.10.8-GCCcore-12.2.0
ml Pysam/0.21.0-GCC-12.2.0

# command to run snakemake (remove -np at end when done validating):
snakemake -s CRAMboozle.snakefile --latency-wait 60 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} -J {cluster.JobName}" -j 40 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"
configfile: "config/cluster_slurm.yaml"

rule all:
    input:
        expand("{results_dir}/{samples}_deidentified.cram", results_dir=config['results_dir'], samples=config['samples'].keys()),
        expand("{results_dir}/{samples}_deidentified.cram.crai", results_dir=config['results_dir'], samples=config['samples'].keys())

rule cramboozle:
    input:
        alignment_file = lambda wildcards: config["samples"][wildcards.samples]
    output:
        deidentified_file = "{results_dir}/{samples}_deidentified.cram",
        index_file = "{results_dir}/{samples}_deidentified.cram.crai"
    params:
        sample_name = "{samples}",
        reference_genome = config['reference_genome'],
        results_dir = config['results_dir'],
        strict_flag = "--strict" if config.get('strict_mode', False) else "",
        unmapped_flag = "--keepunmapped" if config.get('keep_unmapped', False) else "",
        secondary_flag = "--keepsecondary" if config.get('keep_secondary', False) else ""
    log:
        "{results_dir}/logs/{samples}_cramboozle.log"
    shell:
        """
        mkdir -p {params.results_dir}/logs
        python CRAMboozle.py \
            --input {input.alignment_file} \
            --out {output.deidentified_file} \
            --fa {params.reference_genome} \
            {params.strict_flag} \
            {params.unmapped_flag} \
            {params.secondary_flag} \
            2>&1 | tee {log}
        """

rule generate_summary:
    input:
        deidentified_files = expand("{results_dir}/{samples}_deidentified.cram", results_dir=config['results_dir'], samples=config['samples'].keys())
    output:
        summary = "{results_dir}/CRAMboozle_summary.txt"
    params:
        results_dir = config['results_dir']
    run:
        import os
        import datetime
        
        with open(output.summary, 'w') as f:
            f.write("CRAMboozle De-identification Summary\n")
            f.write("=" * 40 + "\n")
            f.write(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total samples processed: {len(input.deidentified_files)}\n\n")
            
            f.write("Processed files:\n")
            for file_path in input.deidentified_files:
                sample_name = os.path.basename(file_path).replace('_deidentified.cram', '')
                file_size = os.path.getsize(file_path) if os.path.exists(file_path) else 0
                f.write(f"  - {sample_name}: {file_path} ({file_size:,} bytes)\n")
            
            f.write(f"\nAll files saved to: {params.results_dir}\n")

