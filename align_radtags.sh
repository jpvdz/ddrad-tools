#!/bin/bash

#SBATCH --time=02:00:00
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --job-name=bowtie2_align_radtags
#SBATCH --array=1-1
#SBATCH --output=%x-%A_%a.out

# load modules
module load Bowtie2/2.4.4-GCC-11.2.0
module load SAMtools/1.10-GCC-9.3.0

# read command line arguments
while getopts 'g:s:o:' OPTION; do
  case $OPTION in
    g)  genome_prefix=${OPTARG}
        ;;
    s)  sample_file=${OPTARG}
        ;;
    o)  output=${OPTARG}
        ;;
    ?)  printf "How to run %s: \n\
                -g : genome filename prefix\n\
                -s : name of sample map file\n\
                -o : output directory name" ${0##*/} >&2
    exit 2
    ;;
  esac
done

# read sample names
samples=()
while read -r sample
do
  samples+=(${sample})
done < ${dir}/info/${sample_file}

# create output directory
mkdir -p ${dir}/aligned/${output}

# get individual names from samples array
sample_name="${samples[${SLURM_ARRAY_TASK_ID}-1]}"

# run alignment
bowtie2 -q -x ${dir}/genomes/${genome_prefix} \
        -1 ${dir}/cleaned/${sample_name}.1.fq.gz \
        -2 ${dir}/cleaned/${sample_name}.2.fq.gz \
        -U ${dir}/cleaned/${sample_name}.rem.1.fq.gz,${dir}/cleaned/${sample_name}.rem.2.fq.gz \
        --very-sensitive \
        -p 12 \
        -S ${dir}/aligned/${output}/${sample_name}.sam \
        &> ${dir}/aligned/${output}/${sample_name}.log

# convert sam to bam
samtools    view -b -q 10 -h -o ${dir}/aligned/${output}/${sample_name}.bam ${dir}/aligned/${output}/${sample_name}.sam

# sort bam
samtools    sort -o ${dir}/aligned/${output}/${sample_name}.bam ${dir}/aligned/${output}/${sample_name}.bam

# calculate statistics for bam
samtools    stats ${dir}/aligned/${output}/${sample_name}.bam -c 1,1000,1 \
            >${dir}/aligned/${output}/${sample_name}.stats

# calculate depth
samtools    depth ${dir}/aligned/${output}/${sample_name}.bam -Q 10 \
            >${dir}/aligned/${output}/${sample_name}.depth

# calculate flagstats
samtools    flagstat ${dir}/aligned/${output}/${sample_name}.bam \
            >${dir}/aligned/${output}/${sample_name}.flagstat
