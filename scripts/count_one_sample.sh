#!/bin/bash
# Usage: bash count_one_sample.sh <sampleID>

sid="$1"
outfile="counts.tsv"

if [[ -z "$sid" ]]; then
    echo "Usage: bash $0 <sampleID>"
    exit 1
fi

# Output header if file does not exist yet
if [[ ! -f "$outfile" ]]; then
    echo -e "SampleID\tPASS\tpon_snvs\tfings\tvaf\tPASS_indel\tpon_indel\tADD\tVAF" > "$outfile"
fi

# Count helper
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

# Directories
snv_dir="snvs_Filtered"
flt_dir="Filtered"
indel_dir="indel_Filtered"

# Collect all counts
fings=$(count_or_na "$snv_dir" "$sid" "_fings")
vaf=$(count_or_na "$snv_dir" "$sid" ".vaf")
shared_snvs=$(count_or_na "$flt_dir" "$sid" "_shared_snvs")
shared_indel=$(count_or_na "$flt_dir" "$sid" "_shared_indel")
pon_indel=$(count_or_na "$flt_dir" "$sid" "_pon_indel")
pon_snvs=$(count_or_na "$flt_dir" "$sid" "_pon_snvs")
add_indel=$(count_or_na "$indel_dir" "$sid" ".RF")
vaf_indel=$(count_or_na "$indel_dir" "$sid" ".AF")

# Append results to file
echo -e "${sid}\t${shared_snvs}\t${pon_snvs}\t${fings}\t${vaf}\t${shared_indel}\t${pon_indel}\t${add_indel}\t${vaf_indel}" >> "$outfile"
