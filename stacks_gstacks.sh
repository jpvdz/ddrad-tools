#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --job-name=stacks_gstacks
#SBATCH --output=%x-%j.out

# load modules
module load Stacks/2.59-foss-2020a

# set stacks directory to script folder
stacks_dir=$(dirname (realpath $0))

# read command line arguments
while getopts 'P:a:o:M:Q:' OPTION; do
  case $OPTION in
    P)  popmap_file=${OPTARG}
        ;;
    a)  alignment=${OPTARG}
        ;;
    o)  output_name=${OPTARG}
        ;;
    M)  model=${OPTARG}
        ;;
    Q)  mapq=${OPTARG}
        ;;
    ?) printf "How to run %s: \n\
                -P : name of popmap file\n\
                -a : alignment directory name (under 'alignments')\n\
                -o : output directory name\n\
                -M : the model to use for variant calling and genotyping\n\
                -Q : minimum mapping quality" ${0##*/} >&2
        exit 2
        ;;
  esac
done

# set output path
out_dir=${stacks_dir}/output/${output_name}
mkdir -p ${out_dir}

# set popmap directory
popmap=${stacks_dir}/popmap/${popmap_file}

# run gstacks
gstacks -I ${alignment} \
        -O ${out_dir} \
        -M ${popmap} \
        --model ${model} \
        --min-mapq ${mapq}
