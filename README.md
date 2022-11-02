# ddRAD tools
A collection of Bash scripts for processing ddRAD data.

## Dependencies
Requires a (Linux) computing cluster with SLURM and the following software:
- Stacks (version 2.59): https://catchenlab.life.illinois.edu/stacks/
- Bowtie2 (version 2.4.4): https://github.com/BenLangmead/bowtie2
- SAMtools (version 1.10): https://github.com/samtools/samtools
- VCFtools (version 0.1.16): https://vcftools.github.io/

The scripts should work if you use the listed versions. They probably work for later versions, but no guarantees.

## How to use
Put the scripts in a folder (e.g., `/ddrad`) along with the following:
- A `/raw` folder containing raw reads to be cleaned.
- A `/genomes` folder containing a Bowtie2-indexed reference genome.
- An `/info` folder containing ancillary files such as a population map or ddRAD library info.

The rest of the folders should be generated automatically. I provide the restriction enzymes in columns 4 and 5 in a `lib-info` 
file (without an extension) that contains some metadata.