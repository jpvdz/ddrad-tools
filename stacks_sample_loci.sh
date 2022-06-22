#!/bin/bash

#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12gb
#SBATCH --job-nam=stacks_sample_loci
#SBATCH --output=%x-%j.out

# load modules
module load Stacks/2.59-foss-2020a
module load VCFtools/0.1.16-GCC-10.2.0

# set stacks directory to script folder
stacks_dir=$(dirname (realpath $0))

# read command line arguments
while getopts 'P:a:o:M:Q:m:r:p:k:' OPTION; do
  case $OPTION in
    P)  popmap_list="$OPTARG"
        ;;
    a)  alignment="$OPTARG"
        ;;
    o)  output_name="$OPTARG"
        ;;
    M)  model="$OPTARG"
        ;;
    Q)  mapq="$OPTARG"
        ;;
    m)  maf="$OPTARG"
        ;;
    r)  min_perc_ind="$OPTARG"
        ;;
    p)  num_pops="$OPTARG"
        ;;
    k)  max_het="$OPTARG"
        ;;
    ?)  printf "How to run %s: \n\
                -P : name of popmap list file\n\
                -a : alignment directory name (under 'alignments')\n\
                -o : output directory name\n\
                -M : the model to use for variant calling and genotyping\n\
                -Q : minimum mapping quality\n\
                -m : minor allele frequency\n\
                -r : minimum percentage of individuals within a population
                     required to process a locus\n\
                -p : minimum number of populations a locus must be present in\n\
                -k : maximum heterozygosity" ${0##*/} >&2
        exit 2
        ;;
  esac
done

# set output path
out_dir=${stacks_dir}/${output_name}
mkdir -p ${outdir}

# declare empty arrays for popmap file
popmaps=()
pops=()

# read filenames and pops from popmap list
while read -r filename pop; do
  popmaps+=(${filename})
  pops+=(${pop})
done < ${stacks_dir}/info/${popmap_list}

# run gstacks
gstacks -I ${stacks_dir}/alignments/${alignment} \
        -O ${out_dir} \
        -M ${stacks_dir}/info/${popmaps[0]} \
        --model ${model} \
        --min-mapq ${mapq}

# create sub-directories for output of STACKS populations
counter=1
for file in ${popmaps[@]}; do

  # create sub-directory for each popmap
  out_subdir=${out_dir}/${pops[${counter}-1]}
  mkdir -p ${out_subdir}

  if [ "${counter}" -eq 1 ]; then
    # run on all pops
    populations -P ${out_dir} \
                -O ${out_subdir} \
                -M ${stacks_dir}/info/${popmaps[${counter}-1]} \
                -r ${min_perc_ind} \
                --max_obs_het ${max_het} \
                --min_maf ${maf} \
                -p ${num_pops} \
                --fasta_loci \
                --fstats \
                --vcf \
                --structure \
                --write_single_snp \
                --ordered_export

    # change to output directory
    cd ${out_subdir}

    # make directories for thinned and filtered output
    mkdir -p filtered
    mkdir -p thinned

    # stats and missing data vcf
    vcf-stats populations.snps.vcf > populations.snps.stats
    vcftools --vcf populations.snps.vcf --missing-indv --out populations.snps
    vcftools --vcf populations.snps.vcf --site-mean-depth --out populations.snps
    vcftools --vcf populations.snps.vcf --geno-depth --out populations.snps

    # exclude sites below 30x coverage and above 300x coverage
    vcftools --vcf populations.snps.vcf --min-meanDP 30 --max-meanDP 300 --recode --stdout > populations.snps.filtered.vcf

    # stats and missing data filtered vcf
    vcf-stats populations.snps.filtered.vcf > populations.snps.filtered.stats
    vcftools --vcf populations.snps.filtered.vcf --missing-indv --out populations.snps.filtered
    vcftools --vcf populations.snps.filtered.vcf --site-mean-depth --out populations.snps.filtered
    vcftools --vcf populations.snps.filtered.vcf --geno-depth --out populations.snps.filtered

    # re-run populations on filtered vcf file
    populations -V populations.snps.filtered.vcf \
                -O filtered \
                -M ${stacks_dir}/info/${popmaps[0]} \
                --fstats

    # thin SNPs within 100kb window
    vcftools --vcf populations.snps.filtered.vcf --thin 100000 --recode --stdout > populations.snps.t100000.vcf

    # stats and missing data thinned vcf
    vcf-stats populations.snps.t100000.vcf > populations.snps.t100000.stats
    vcftools --vcf populations.snps.t100000.vcf --missing-indv --out populations.snps.t100000
    vcftools --vcf populations.snps.t100000.vcf --site-mean-depth -out populations.snps.t100000
    vcftools --vcf populations.snps.t100000.vcf --geno-depth --out populations.snps.t100000

    # re-run populations on thinned vcf file
    populations -V populations.snps.t100000.vcf \
                -O thinned \
                -M ${stacks_dir}/info/${popmaps[0]} \
                --fstats

    # create whitelist of thinned SNPs
    vcf-query populations.snps.t100000.vcf -f '%ID\n' | cut -d ':' -f1 | sort -n > ${out_dir}/whitelist.tsv
    let counter+=1
    continue
  else
    populations -P ${out_dir} \
                -O ${out_subdir} \
                -M ${stacks_dir}/info/${popmaps[${counter}-1]} \
                -W ${out_dir}/whitelist.tsv \
                -r ${min_perc_ind} \
                --max_obs_het ${max_het} \
                -p 1 \
                --fasta_loci \
                --fstats \
                --vcf \
                --structure \
                --write_single_snp \
                --ordered_export

    # randomly select loci
    vcf-query ${out_subdir}/populations.snps.vcf -f '%ID\n' | sort -R | head -n6000 | sort -n \
    > ${out_subdir}/selected_snps.tsv
  fi
  let counter+=1
done

# concatenate whitelists, sort and keep unique IDs;
# then randomly sample loci from unique IDs, sort and store in whitelist
cat ${out_dir}/${pops[1]}/selected_snps.tsv ${out_dir}/${pops[2]}/selected_snps.tsv \
  ${out_dir}/${pops[3]}/selected_snps.tsv | cut -d ':' -f1 | sort | uniq | shuf | head -n10000 | sort -n \
  > ${out_dir}/target_snp_cat_ids.tsv

# re-run populations using randomly sampled SNPs
mkdir ${out_dir}/final
populations -P ${out_dir} \
            -O ${out_dir}/final \
            -W ${out_dir}/target_snp_cat_ids.tsv \
            -M ${stacks_dir}/info/${popmaps[0]} \
            -r ${min_perc_ind} \
            --max_obs_het ${max_het} \
            --min_maf ${maf} \
            -p 1 \
            --fasta_loci \
            --fstats \
            --vcf \
            --structure \
            --write_single_snp \
            --ordered_export

# extract chromosome, position and ID for each randomly selected SNP from final VCF
vcf-query ${out_dir}/final/populations.snps.vcf -f '%CHROM\t%POS\t%ID\n' \
  > ${out_dir}/target_snp_chrom_pos.tsv
