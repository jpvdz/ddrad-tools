#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --job-name=stacks_populations
#SBATCH --output=%x-%j.out

# load modules
module load Stacks/2.59-foss-2020a
module load VCFtools/0.1.16-GCC-10.2.0

# set working directory to script folder
dir=$(dirname (realpath $0))

# read command line arguments
while getopts 'P:i:o:Q:m:r:p:R:k:t:' OPTION; do
  case $OPTION in
    P)  popmapfile=${OPTARG}
        ;;
    i)  input_name=${OPTARG}
        ;;
    o)  output_name=${OPTARG}
        ;;
    Q)  mapq=${OPTARG}
        ;;
    m)  maf=${OPTARG}
        ;;
    r)  min_perc_ind=${OPTARG}
        ;;
    p)  min_perc_pop=${OPTARG}
        ;;
    R)  min_perc_ind_pop=${OPTARG}
        ;;
    k)  max_het=${OPTARG}
        ;;
    t)  thin_window_size=${OPTARG}
        ;;
    ?)  printf "How to run %s: \n\
                -P : name of popmap list file\n\
                -i : alignment directory name (under 'alignments')\n\
                -o : output directory name\n\
                -Q : minimum mapping quality\n\
                -m : minor allele frequency\n\
                -r : minimum percentage of individuals within a population
                     required to process a locus\n\
                -p : minimum percentage of populations required to process a\n\
                     locus\n\
                -R : minimum percentage of individuals across populations\n\
                     required to process a locus\n\
                -k : maximum heterozygosity\n\
                -t : thinning window size" ${0##*/} >&2
        ;;
  esac
done

# set input and output path
input_dir=${dir}/output/${input_name}
out_dir=${input_dir}/${output_name}
mkdir -p ${out_dir}

# set popmap directory
popmap=${dir}/popmap/${popmap_file}

# run populations
populations -P ${input_dir} \
            -M ${popmap} \
            -O ${out_dir} \
            -r ${min_perc_ind} \
            -p ${min_perc_pop} \
            -R ${min_perc_ind_pop} \
            --min-maf ${maf} \
            --max-obs-het ${max_het} \
            --fasta_samples \
            --fstats \
            --hwe \
            --vcf \
            --genepop \
            --radpainter

# go to output directory
cd ${out_dir}

# stats and missing data vcf
vcf-stats populations.snps.vcf > populations.snps.stats
vcftools --vcf populations.snps.vcf --missing-indv --out populations.snps
vcftools --vcf populations.snps.vcf --site-mean-depth --out populations.snps
vcftools --vcf populations.snps.vcf --geno-depth --out populations.snps
vcftools --vcf populations.snps.vcf --depth --out populations.snps
vcftools --vcf populations.snps.vcf --het --out populations.snps
vcftools --vcf populations.snps.vcf --indv-freq-burden --out populations.snps
vcftools --vcf populations.snps.vcf --singletons --out populations.snps

# do thinning using VCFtools
vcftools --vcf populations.snps.vcf --thin ${thin_window_size} --recode --stdout \
  > populations.snps.t${thin_window_size}.vcf
populations -V populations.snps.t${thin_window_size}.vcf \
            -M ${popmap} \
            -O ${out_dir} \
            --structure \
            --treemix

# filter loci out of HWE (P < 0.05) using VCFtools
vcftools --vcf populations.snps.t${thin_window_size}.vcf --hwe 0.05 --recode --stdout \
  > populations.snps.t${thin_window_size}.hwe.vcf
populations -V populations.snps.t${thin_window_size}.hwe.vcf \
            -M ${popmap} \
            -O ${out_dir} \
            --structure \
            --treemix
