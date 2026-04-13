#!/bin/bash
# added RepeatFilter count as the final count for all other indel filters, these include: NoAlt in Normal, SB tolerant, SB zero and RepeatFilter
# place this script in the top project folder


echo -e "SampleID\tPASS\tpon_snvs\tfings\tvaf\tPASS_indel\tpon_indel\tADD\tVAF\t"
# Get unique sample IDs from all .vcf files in snvs_Filtered and Filtered
sample_ids=$(find snvs_Filtered Filtered indel_Filtered -maxdepth 1 -type f -name "*.vcf" -printf "%f\n" \
  | sed -E 's/(_fings|\.vaf|_shared_snvs|_shared_indel|_pon_indel|_pon_snvs|\.RF|\.AF|\.SBtolerant|\.NoAlt|\.SBzero)\.vcf$//' \
  | sort -u)

# echo $sample_ids

# Count function
count_or_na() {
    local dir="$1"
    local sid="$2"
    local pattern="$3"
    local file="${dir}/${sid}${pattern}.vcf"

    if [[ -f "$file" ]]; then
        grep -vc '^#' "$file"
    else
        echo "NA"
    fi
}

# Loop through samples
for sid in $sample_ids; do
    fings=$(count_or_na "snvs_Filtered" "$sid" "_fings")
    vaf=$(count_or_na "snvs_Filtered" "$sid" ".vaf")
    shared_snvs=$(count_or_na "Filtered" "$sid" "_shared_snvs")
    shared_indel=$(count_or_na "Filtered" "$sid" "_shared_indel")
    pon_indel=$(count_or_na "Filtered" "$sid" "_pon_indel")
    pon_snvs=$(count_or_na "Filtered" "$sid" "_pon_snvs")
    add_indel=$(count_or_na "indel_Filtered" "$sid" ".RF")
    vaf_indel=$(count_or_na "indel_Filtered" "$sid" ".AF")
    echo -e "${sid}\t${shared_snvs}\t${pon_snvs}\t${fings}\t${vaf}\t${shared_indel}\t${pon_indel}\t${add_indel}\t${vaf_indel}"
done
