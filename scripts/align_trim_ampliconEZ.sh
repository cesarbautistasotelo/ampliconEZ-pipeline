#!/bin/bash
#BSUB -J align_trim_ampliconEZ
#BSUB -q short
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -W 2:00
#BSUB -o align_trim_ampliconEZ.out
#BSUB -e align_trim_ampliconEZ.err

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 REF FQ_DIR OUT_DIR [PRIMER_TSV] [TRIM_MODE]"
    exit 1
fi

REF="$1"
FQ_DIR="$2"
OUT_DIR="$3"
PRIMER_TSV="${4:-}"
TRIM_MODE="${5:-none}"

module load bwa/0.7.17
module load samtools/1.3
module load cutadapt/4.1

mkdir -p "$OUT_DIR"
[[ "$TRIM_MODE" != "none" ]] && mkdir -p "$OUT_DIR/trimmed"

case "$TRIM_MODE" in
    none|trimmed|both) ;;
    *) echo "[ERROR] Invalid TRIM_MODE"; exit 1 ;;
esac

declare -a GENES
declare -A FWD_PRIMERS
declare -A REV_PRIMERS

if [[ "$TRIM_MODE" != "none" ]]; then
    while read -r gene fwd rev _rest; do
        # Skip empty lines
        [[ -z "$gene" ]] && continue
        # Skip commented lines or header row starting with "gene"
        [[ "$gene" == \#* ]] && continue
        [[ "$gene" == "gene" ]] && continue

        GENES+=("$gene")
        FWD_PRIMERS["$gene"]="$fwd"
        REV_PRIMERS["$gene"]="$rev"
    done < "$PRIMER_TSV"

    echo "[INFO] Loaded ${#GENES[@]} primer entries from $PRIMER_TSV"
fi

get_gene() {
    local sample="$1"
    for g in "${GENES[@]}"; do
        [[ "$sample" == *"$g"* ]] && echo "$g" && return
    done
    echo ""
}

for R1 in "$FQ_DIR"/*_R1_001.fastq.gz; do
    [[ ! -f "$R1" ]] && continue

    SAMPLE=$(basename "$R1" | sed 's/_R1_001.fastq.gz//')
    R2="$FQ_DIR/${SAMPLE}_R2_001.fastq.gz"

    UNTRIM_BAM="$OUT_DIR/${SAMPLE}_untrimmed.bam"
    TRIM_R1="$OUT_DIR/trimmed/${SAMPLE}_R1_trimmed.fastq.gz"
    TRIM_R2="$OUT_DIR/trimmed/${SAMPLE}_R2_trimmed.fastq.gz"
    TRIM_BAM="$OUT_DIR/${SAMPLE}_trimmed.bam"

    echo "=== Processing $SAMPLE ==="

    if [[ "$TRIM_MODE" == "none" || "$TRIM_MODE" == "both" ]]; then
        if [[ ! -f "$UNTRIM_BAM" || ! -f "${UNTRIM_BAM}.bai" ]]; then
            bwa mem -t 4 "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "$UNTRIM_BAM"
            samtools index "$UNTRIM_BAM"
        else
            echo "[SKIP] $SAMPLE untrimmed exists"
        fi
    fi

    if [[ "$TRIM_MODE" != "none" ]]; then

        gene=$(get_gene "$SAMPLE")
        if [[ -z "$gene" ]]; then
            echo "[WARN] No gene match for $SAMPLE"
            continue
        fi

        FWD="${FWD_PRIMERS[$gene]}"
        REV="${REV_PRIMERS[$gene]}"

        if [[ ! -f "$TRIM_BAM" || ! -f "${TRIM_BAM}.bai" ]]; then
            cutadapt -g "$FWD" -G "$REV" \
                -o "$TRIM_R1" -p "$TRIM_R2" \
                "$R1" "$R2" > "$OUT_DIR/trimmed/${SAMPLE}_cutadapt.log"

            bwa mem -t 4 "$REF" "$TRIM_R1" "$TRIM_R2" | samtools sort -@ 4 -o "$TRIM_BAM"
            samtools index "$TRIM_BAM"
        else
            echo "[SKIP] $SAMPLE trimmed exists"
        fi
    fi

done

echo "[DONE] Output in $OUT_DIR"
