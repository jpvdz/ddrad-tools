#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --partition=regular
#SBATCH --job-name=process_radtags
#SBATCH --output=log/%x-%j.out

# load modules
module load Stacks/2.59-foss-2020a

# read command line arguments
while getopts 'L:o:' OPTION; do
  case $OPTION in
    L)  library="$OPTARG"
        ;;
    o)  output="$OPTARG"
        ;;
    ?)  printf "How to run %s: \n\
                -L : library name\n\
                -o : output directory name" ${0##*/} >&2
        exit 2
        ;;
  esac
done
shift $(($OPTIND - 1))

# set working directory to script folder
dir=$(dirname (realpath $0))

# create output dir
mkdir -p ${dir}/cleaned/${output}

# read restriction enzymes from lib-info file
re_one="$(grep ${library} ${dir}/info/lib-info | cut -f 4)"
re_two="$(grep ${library} ${dir}/info/lib-info | cut -f 5)"

# remove poor quality reads
process_radtags -i gzfastq \
                -P -p ${dir}/raw/ \
                -b ${dir}/info/${library}_barcodes --inline_index \
                --renz-1 ${re_one} --renz-2 ${re_two} \
                -o ${dir}/cleaned/${output} \
                -c -q
