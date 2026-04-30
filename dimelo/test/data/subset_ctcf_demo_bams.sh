#!/usr/bin/env bash
set -euo pipefail

# ---- inputs ----
SRC_DIR="/Users/ngamarra/sherlock_oak/data/20250110_NG_one_pot"
OUT_DIR="/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/test/data"

PEAK_BED="${OUT_DIR}/ctcf_demo_peak.bed"
NOT_PEAK_BED="${OUT_DIR}/ctcf_demo_not_peak.bed"

BAMS=(
  "barcode17.merged.sorted.bam"
  "barcode18.merged.sorted.bam"
)

# ---- checks ----
command -v samtools >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found"; exit 1; }

mkdir -p "${OUT_DIR}"

for f in "${PEAK_BED}" "${NOT_PEAK_BED}"; do
  [[ -s "$f" ]] || { echo "ERROR: missing BED: $f"; exit 1; }
done

# ---- make combined target BED ----
TARGET_BED="${OUT_DIR}/ctcf_demo_peak_and_not_peak.sorted.merged.bed"

cat "${PEAK_BED}" "${NOT_PEAK_BED}" \
  | awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && NF >= 3 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "${TARGET_BED}"

echo "Wrote target regions:"
echo "  ${TARGET_BED}"
wc -l "${TARGET_BED}"

# ---- subset each BAM ----
for bam_name in "${BAMS[@]}"; do
  in_bam="${SRC_DIR}/${bam_name}"

  [[ -s "$in_bam" ]] || { echo "ERROR: missing BAM: $in_bam"; exit 1; }

  # Ensure source BAM is indexed. This may write to the SSHFS-mounted OAK dir.
  if [[ ! -s "${in_bam}.bai" && ! -s "${in_bam%.bam}.bai" ]]; then
    echo "Indexing source BAM: ${in_bam}"
    samtools index "${in_bam}"
  fi

  stem="${bam_name%.bam}"
  tmp_bam="${OUT_DIR}/${stem}.ctcf_demo.tmp.bam"
  out_bam="${OUT_DIR}/${stem}.ctcf_demo.sorted.bam"

  echo
  echo "Subsetting:"
  echo "  input:  ${in_bam}"
  echo "  output: ${out_bam}"

  # -L keeps all alignments overlapping the BED intervals.
  # -b outputs BAM.
  samtools view -b -L "${TARGET_BED}" "${in_bam}" \
    | samtools sort -o "${tmp_bam}" -

  mv "${tmp_bam}" "${out_bam}"
  samtools index "${out_bam}"

  echo "Done:"
  samtools idxstats "${out_bam}" | awk '$3 > 0 {print}'
done

echo
echo "Finished. Outputs in:"
echo "  ${OUT_DIR}"
ls -lh "${OUT_DIR}"/*.ctcf_demo.sorted.bam*
