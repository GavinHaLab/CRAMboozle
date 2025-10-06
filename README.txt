Robert D Patton, 10/06/2025
rpatton@fredhutch.org
This is a modified version of BAMboozle for use by the Gavin Ha lab at Fred Hutch Cancer Center
Changes from the original (v0.5.0) include:

1. Adapted to output CRAM format by default, and take in either BAM or CRAM format (dynamic)
2. Removed the --p ('Number of processes to use') option; now uses all available CPUs by default
3. Integrated into a simple snakemake to speed up the pipeline process