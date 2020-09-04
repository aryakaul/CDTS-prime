#!/bin/bash
set -e

extract_highconf_lowcov_regions() {
    start=$(date +%s)
    mkdir -p $1/qc
    mkdir -p $1/coverage
    log "Downloading and splitting illumina's platinum genome"
    if [ ! -f $1/qc/confidentRegions_hg38.starch ]; then
        wget -q -N -P $1/qc "https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/ConfidentRegions.bed.gz"
        wget -q -N -P $1/qc "https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/small_variants/ConfidentRegions.bed.gz.tbi"
        zcat $1/qc/ConfidentRegions.bed.gz | sed 's/^chr//g' | sort-bed - > $1/qc/confidentRegions_hg38.bed
        while read -r chroms; do
            bedextract $chroms $1/qc/confidentRegions_hg38.bed > $1/qc/confidentRegions_hg38_chrom"$chroms".bed &
        done < $1/chromList.txt
        starch $1/qc/confidentRegions_hg38.bed > $1/qc/confidentRegions_hg38.starch &
        wait
        rm $1/qc/ConfidentRegions.bed.gz*
        rm $1/qc/confidentRegions_hg38.bed
    else
        log "  Already downloaded... skipping."
    fi
    
    log "Computing low coverage regions"
    if [ ! -s $1/coverage/low_coverage_regions.starch ]; then
        zcat $2/coverage/gnomad.genomes.r3.0.coverage.summary.tsv.bgz \
            | sed '1d' \
            | sed 's/^chr//g' \
            | sed 's/:/\t/g' \
            | awk -v OFS='\t' -v thresh=$3 '{if($6 < thresh) {print $1, $2-1, $2 }}' \
            | sort-bed - \
            | bedops -m - \
            > $1/coverage/low_coverage_regions.bed
        while read -r chroms; do
            bedextract $chroms $1/coverage/low_coverage_regions.bed \
                > $1/coverage/low_coverage_chrom"$chroms".bed &
        done < $1/chromList.txt
        starch $1/coverage/low_coverage_regions.bed > $1/coverage/low_coverage_regions.starch &
        wait
    else
        log "  Already computed... skipping."
    fi
    
    log "Intersecting low coverage and high quality regions"
    while read -r chroms; do
        if [ ! -s $1/qc/lowcov_highqual_chrom"$chroms".starch ]; then
            bedops --intersect $1/coverage/low_coverage_chrom"$chroms".bed $1/qc/confidentRegions_hg38_chrom"$chroms".bed \
                | starch - \
                > $1/qc/lowcov_highqual_chrom"$chroms".starch &
        else
            continue
        fi
    done < $1/chromList.txt
    wait
   
    log "Splitting into autosomal + pseudoautosomal (PAR)  and chrX non-PAR"
    echo -e "X\t10001\t2781479\nX\t155701383\t156030895" > $1/PAR_hg38.bed
    echo -e "X\t0\t10001\nX\t2781479\t155701383\nX\t156030895\t156040895" > $1/chrXnotPAR_hg38.bed
    [ ! -s $1/qc/confidentRegions_nolowcov_hg38_chrXnotPAR.starch ] &&
        bedops --intersect $1/qc/confidentRegions_hg38.starch $1/chrXnotPAR_hg38.bed \
            | bedops --difference - $1/qc/lowcov_highqual_chromX.starch \
            | starch - \
            > $1/qc/confidentRegions_nolowcov_hg38_chrXnotPAR.starch

    [ ! -s $1/qc/confidentRegions_nolowcov_hg38_autosomesAndPAR.starch ] &&
        bedops --difference $1/qc/confidentRegions_hg38.starch $1/chrXnotPAR_hg38.bed \
            | bedops --difference - <(sort-bed $1/qc/lowcov_highqual_chrom*.starch) \
            | starch - \
            > $1/qc/confidentRegions_nolowcov_hg38_autosomesAndPAR.starch
    
    end=$(date +%s)
    runtime=$((end-start))
    log "Time Taken" "$runtime seconds"
}


find_passing_snvs() {
    start=$(date +%s)
    mkdir -p $1/af
    while read -r chroms; do
        log "Filtering chrom $chroms"
        if [ ! -s $1/af/gnomad_PASS_SNV_chr${chroms}_af.starch ]; then
            # changing this to something more agnostic is prefereable
            tabix $2/vcf/gnomad.genomes.r3.0.sites.chr${chroms}.vcf.bgz chr${chroms} \
                | sed 's/^chr//g' \
                | awk -v OFS='\t' -F'\t' '{if ($7=="PASS" && length($4)==1 && \
                (length($5)==1 || (length($5)==3 && $5~",") || (length($5)==5 && $5~",.,"))) print $1, $2-1, $2, $4, $5, $8}' \
                | sed -e 's/\tAC=/\t/' -e 's/;AF=/\t/' -e 's/;AN=/\t/' -e 's/;.*//'  \
                | awk -v OFS='\t' -F'\t' '{split($5,alt,","); split($6,ac,","); split($8,af,","); \
                    for (i=1;i<=length(alt);i++) print $1, $2, $3, $4, alt[i], ac[i], af[i], $7}' \
                | python $3/reformat_ACAFAN_1linePerPos_addFlags.py \
                | starch - \
                > $1/af/gnomad_PASS_SNV_chr${chroms}_af.starch &
        else
            log "  Already filtered... skipping."
        fi
    done < $1/chromList.txt
    wait
    end=$(date +%s)
    runtime=$(((end-start)/60))
    log "Time Taken" "$runtime minutes"
}

find_passing_snvs_specific() {
    start=$(date +%s)
    CHOICE=$4
    mkdir -p $1/af/$CHOICE
    while read -r chroms; do
        if [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af.starch ]; then
            log "Filtering chrom $chroms for $CHOICE"
            # changing this to something more agnostic is prefereable
            tabix $2/vcf/gnomad.genomes.r3.0.sites.chr${chroms}.vcf.bgz chr${chroms} \
                | sed 's/^chr//g' \
                | awk -v OFS='\t' -F'\t' '{if ($7=="PASS" && length($4)==1 && \
                (length($5)==1 || (length($5)==3 && $5~",") || (length($5)==5 && $5~",.,"))) print $1, $2-1, $2, $4, $5, $8}' \
                | sed -e "s/;AC_$CHOICE=/\t/" -e "s/;AF_$CHOICE=/\t/" -e "s/;AN_$CHOICE=/\t/" \
                | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $7, $8, $9}' \
                | sed -e 's/;.*//' \
                | awk '{if ($6!=0){print $0}}' \
                | awk -v OFS='\t' -F'\t' '{split($5,alt,","); split($6,ac,","); split($8,af,","); \
                    for (i=1;i<=length(alt);i++) print $1, $2, $3, $4, alt[i], ac[i], af[i], $7}' \
                | awk -F'\t' -v OFS='\t' '{if(NF==8){print $0}}' \
                | python $3/reformat_ACAFAN_1linePerPos_addFlags.py \
                | starch - \
                > $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af.starch &
        else
            log "  Already filtered... skipping."
        fi
    done < $1/chromList.txt
    wait
    end=$(date +%s)
    runtime=$(((end-start)/60))
    log "Time Taken" "$runtime minutes"
}

find_high_missingness_specific() {
    start=$(date +%s)
    CHOICE=$2
    nball=1000
    while read -r chroms; do
        log "Filtering chrom $chroms"
        if [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_high_missingness.starch ]; then
            unstarch $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af.starch \
                | awk -v OFS='\t' -F'\t' -v nb=$nball '{if($13 < nb) {print $1,$2,$3}}' \
                | bedops -m - \
                | starch - \
                > $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_high_missingness.starch &
        else
            log "  Already filtered... skipping."
        fi
    done < $1/chromList.txt
    wait

    log "Splitting into pseudoautosomal regions and nonPAR"
    [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chrXPAR_high_missingness.starch ] &&
        bedops --intersect $1/af/$CHOICE/gnomad_PASS_SNV_chrX_high_missingness.starch $1/PAR_hg38.bed \
            | starch - \
            > $1/af/$CHOICE/gnomad_PASS_SNV_chrXPAR_high_missingness.starch
    
    [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chrXnotPAR_high_missingness.starch ] &&
        bedops --intersect $1/af/$CHOICE/gnomad_PASS_SNV_chrX_high_missingness.starch $1/chrXnotPAR_hg38.bed \
            | starch - \
            > $1/af/$CHOICE/gnomad_PASS_SNV_chrXnotPAR_high_missingness.starch

    [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_autosomesAndPAR_high_missingness.starch ] &&
        starchcat $1/af/$CHOICE/gnomad_PASS_SNV_chr[1-9]*_high_missingness.starch $1/af/$CHOICE/gnomad_PASS_SNV_chrXPAR_high_missingness.starch \
            > $1/af/$CHOICE/gnomad_PASS_SNV_autosomesAndPAR_high_missingness.starch 

    log "Removing positions with high missingness"
    mkdir -p $1/qc/$CHOICE
    [ ! -s $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_autosomesAndPAR.starch ] &&
        bedops --difference $1/qc/confidentRegions_nolowcov_hg38_autosomesAndPAR.starch $1/af/$CHOICE/gnomad_PASS_SNV_autosomesAndPAR_high_missingness.starch \
        | starch - \
        > $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_autosomesAndPAR.starch 
    
    [ ! -s $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_chrXnotPAR.starch ] &&
        bedops --difference $1/qc/confidentRegions_nolowcov_hg38_chrXnotPAR.starch $1/af/$CHOICE/gnomad_PASS_SNV_chrXnotPAR_high_missingness.starch \
        | starch - \
        > $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_chrXnotPAR.starch
    
    end=$(date +%s)
    runtime=$(((end-start)/60))
    log "Time Taken" "$runtime minutes"
}

extract_snvs() {
    start=$(date +%s)
    while read -r chroms; do
        log "Extracting chrom $chroms"
        if [ ! -s $1/af/gnomad_PASS_SNV_chr${chroms}_af_highconf.starch ]; then
            bedmap --echo --skip-unmapped $1/af/gnomad_PASS_SNV_chr${chroms}_af.starch $1/qc/confidentRegions_nolowcov_nohighmissingness_hg38_autosomesAndPAR.starch \
            | starch - \
            > $1/af/gnomad_PASS_SNV_chr${chroms}_af_highconf.starch &
        else
            log "  Already downloaded... skipping."
        fi
    done < $1/chromList.txt
    wait
   
    log "Recomputing for autosomal regions and non-PAR regions"
    [ ! -s $1/af/gnomad_PASS_SNV_chrXPAR_af_highconf.starch ] &&
        bedmap --echo --skip-unmapped $1/af/gnomad_PASS_SNV_chrX_af.starch $1/qc/confidentRegions_nolowcov_nohighmissingness_hg38_autosomesAndPAR.starch \
        | starch - \
        > $1/af/gnomad_PASS_SNV_chrXPAR_af_highconf.starch 
    
    [ ! -s $1/af/gnomad_PASS_SNV_chrXnotPAR_af_highconf.starch ] &&
        bedmap --echo --skip-unmapped $1/af/gnomad_PASS_SNV_chrX_af.starch $1/qc/confidentRegions_nolowcov_nohighmissingness_hg38_chrXnotPAR.starch \
        | starch - \
        > $1/af/gnomad_PASS_SNV_chrXnotPAR_af_highconf.starch 
    
    end=$(date +%s)
    runtime=$(((end-start)/60))
    log "Time Taken" "$runtime minutes"
}

extract_snvs_specific() {
    start=$(date +%s)
    CHOICE=$2
    while read -r chroms; do
        log "Extracting chrom $chroms"
        if [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af_highconf.starch ]; then
            bedmap --echo --skip-unmapped $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af.starch $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_autosomesAndPAR.starch \
            | starch - \
            > $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af_highconf.starch &
        else
            log "  Already downloaded... skipping."
        fi
    done < $1/chromList.txt
    wait
   
    log "Recomputing for autosomal regions and non-PAR regions"
    [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chrXPAR_af_highconf.starch ] &&
        bedmap --echo --skip-unmapped $1/af/$CHOICE/gnomad_PASS_SNV_chrX_af.starch $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_autosomesAndPAR.starch \
        | starch - \
        > $1/af/$CHOICE/gnomad_PASS_SNV_chrXPAR_af_highconf.starch 
    
    [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chrXnotPAR_af_highconf.starch ] &&
        bedmap --echo --skip-unmapped $1/af/$CHOICE/gnomad_PASS_SNV_chrX_af.starch $1/qc/$CHOICE/confidentRegions_nolowcov_nohighmissingness_hg38_chrXnotPAR.starch \
        | starch - \
        > $1/af/$CHOICE/gnomad_PASS_SNV_chrXnotPAR_af_highconf.starch 
    
    end=$(date +%s)
    runtime=$(((end-start)/60))
    log "Time Taken" "$runtime minutes"
}

update_allele_specific() {
    start=$(date +%s)
    CHOICE=$2
    seq 1 22 | cat - <(echo -e "XPAR\nXnotPAR") > $1/chromList_PARsplit.txt
    while read -r chroms; do
        log "Extracting chrom $chroms"
        if [ ! -s $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af_highconf_wrongref_idx_newref.starch ]; then
            unstarch $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af_highconf.starch \
                | awk -v OFS='\t' -F'\t' '{if($4 != $10){ \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $2-3, $2-2, $10, 1; \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $2-2, $2-1, $10, 2; \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $2-1, $2,   $10, 3; \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $2,   $3,   $10, 4; \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $3,   $3+1, $10, 5; \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $3+1, $3+2, $10, 6; \
                printf "%s\t%s\t%s\t%s\t%s\n", $1, $3+2, $3+3, $10, 7;}}' \
                | sort-bed - \
                | starch - \
                > $1/af/$CHOICE/gnomad_PASS_SNV_chr${chroms}_af_highconf_wrongref_idx_newref.starch &
        else
            log "  Already extracted... skipping."
        fi
    done < $1/chromList_PARsplit.txt
    wait
    end=$(date +%s)
    runtime=$((end-start))
    log "Time Taken" "$runtime seconds"
}

log() {
    lcol='\033[1;33m' lcol2='\033[1;36m' lclr='\033[m'
    printf '%b%s %b%s%b %s\n' \
        "$lcol" "${3:-->}" "${lclr}${2:+$lcol2}" "$1" "$lclr" "$2" >&2
}

die() {
    log "$1" "$2" "${3:-ERROR}"
    exit 2
}

main() {
    OUTDIR=""
    COVERAGE=0.033
    WD="."
    OPTIONS="male female"
    for i in "$@"; do
    case $i in
        -o|--outdir)
        OUTDIR="$2"
        shift
        shift
        ;;
        -g|--gnomaddir)
        GNOMADDIR="$2"
        shift
        shift
        ;;
        -c|--coverage)
        COVERAGE="$2"
        shift
        shift
        ;;
        -s|--scriptdir)
        WD="$2"
        shift
        shift
        ;;
        -ch|--choice)
        OPTIONS="$2"
        shift
        shift
        ;;
        -h|--help)
        usage
        exit
        ;;
        *)
        ;;
    esac
    done
    [ -z $GNOMADDIR ] && GNOMADDIR="$OUTDIR"
    log "Performing QC" "..." "STEP 2"
    
    log "Computing high-confidence, low coverage regions." "..." "STEP 2A"
    extract_highconf_lowcov_regions $OUTDIR $GNOMADDIR $COVERAGE && log "SUCCESSFUL" "High Confident, Low Coverage Starch files generated" "STEP 2A" || die "Step 2A failed"
  
    for CHOICE in $OPTIONS; do
        log "Computing SNVs that pass filters for $CHOICE" "..." "STEP 2B"
        find_passing_snvs_specific $OUTDIR $GNOMADDIR $WD $CHOICE && log "SUCCESSFUL" "SNVs filtered" "STEP 2B" || die "Step 2B failed"
    done

    log "Computing SNVs that are present in more than 1000 alleles" "..." "STEP 2C"
    for CHOICE in $OPTIONS; do 
        find_high_missingness_specific $OUTDIR $CHOICE && log "SUCCESSFUL" "SNVs filtered" "STEP 2C" || die "Step 2C failed"
    done

    log "Extract SNVs in high confidence, high coverage regions" "..." "STEP 2D"
    for CHOICE in $OPTIONS; do 
        extract_snvs_specific $OUTDIR $CHOICE && log "SUCCESSFUL" "SNVs filtered" "STEP 2D" || die "Step 2D failed"
    done
    
    log "Updating hg38 reference allele" "..." "STEP 2E"
    for CHOICE in $OPTIONS; do 
        update_allele_specific $OUTDIR $CHOICE && log "SUCCESSFUL" "Reference updated." "STEP 2E" || die "Step 2E failed"
    done
    
    log "QC Complete" "..." "STEP 2"
    
    
}
main "$@"
