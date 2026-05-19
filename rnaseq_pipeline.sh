#!/usr/bin/bash
# =============================================================================
# RNA-Seq Pipeline — TFM MUBBC (UAM)
# Author : Ismael de la Iglesia San Sebastián
# Description: End-to-end RNA-Seq pipeline using STAR for alignment and
#              splice-junction-aware quantification (GeneCounts + chimeric
#              read detection). Reference genome: GRCh38 / Ensembl release 108.
#
# Pipeline stages (uncomment to activate):
#   1. Reference genome & annotation download + STAR index generation
#   2. FASTQ lane merging
#   3. Pre-trimming QC  (FastQC + MultiQC)
#   4. Adapter / quality trimming  (Trim Galore)
#   5. Post-trimming QC (FastQC + MultiQC)
#   6. First-pass STAR alignment  (splice-junction discovery)
#   7. Second-pass STAR alignment (final BAMs + GeneCounts + chimeric reads)
#
# Requirements: STAR ≥2.7, FastQC, MultiQC, Trim Galore, rename (util-linux)
# Usage       : bash rnaseq_pipeline.sh
# =============================================================================

set -euo pipefail          # Exit on error, unset variable, or pipe failure
IFS=$'\n\t'                # Safer word splitting

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
HOME_DIR=$(pwd)
NTHREADS=16

STAR_INDEX_DIR="./star_index/hg38"
GENOME_DIR="${STAR_INDEX_DIR}/genome"
ANNOTATION_DIR="${STAR_INDEX_DIR}/annotation"
SJ_INDEX_DIR="${STAR_INDEX_DIR}/SJ_index"

RESULTS_DIR="./RESULTS"
FASTQC_PRE_DIR="${RESULTS_DIR}/fastqc-pre"
FASTQC_POST_DIR="${RESULTS_DIR}/fastqc-post"
TRIMMED_DIR="${RESULTS_DIR}/trimmedseq"
STAR_FIRST_DIR="${RESULTS_DIR}/STAR"
STAR_FINAL_DIR="${RESULTS_DIR}/STAR_final"

ENSEMBL_FASTA="https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
ENSEMBL_GTF="http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die()  { log "ERROR: $*" >&2; exit 1; }

require_tool() {
    command -v "$1" &>/dev/null || die "Required tool not found: $1"
}

# ---------------------------------------------------------------------------
# Stage 1 — Reference download & STAR genome index
# ---------------------------------------------------------------------------
stage1_build_index() {
    log "Stage 1: Building STAR genome index"

    mkdir -p "${GENOME_DIR}" "${ANNOTATION_DIR}"

    log "  Downloading genome FASTA..."
    wget -q "${ENSEMBL_FASTA}" -P "${GENOME_DIR}"
    gunzip "${GENOME_DIR}"/*.gz

    log "  Downloading GTF annotation..."
    wget -q "${ENSEMBL_GTF}" -P "${ANNOTATION_DIR}"
    gunzip "${ANNOTATION_DIR}"/*.gz

    log "  Running STAR genomeGenerate..."
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX_DIR}" \
        --genomeFastaFiles "${GENOME_DIR}"/* \
        --sjdbGTFfile "${ANNOTATION_DIR}"/* \
        --runThreadN "${NTHREADS}"

    log "Stage 1 complete."
}

# ---------------------------------------------------------------------------
# Stage 2 — Merge paired-end lanes
# ---------------------------------------------------------------------------
stage2_merge_lanes() {
    log "Stage 2: Merging FASTQ lanes"

    ls -1 ./*R1*.gz \
        | awk -F 'L' '{print $1}' \
        | sort -u > ID

    while IFS= read -r sample; do
        log "  Merging R1 for ${sample}"
        cat "${sample}L001_R1_001.fastq.gz" "${sample}L002_R1_001.fastq.gz" \
            > "${sample}L001_R1_001_merged.fastq.gz"
        rm  "${sample}L001_R1_001.fastq.gz" "${sample}L002_R1_001.fastq.gz"

        log "  Merging R2 for ${sample}"
        cat "${sample}L001_R2_001.fastq.gz" "${sample}L002_R2_001.fastq.gz" \
            > "${sample}L001_R2_001_merged.fastq.gz"
        rm  "${sample}L001_R2_001.fastq.gz" "${sample}L002_R2_001.fastq.gz"
    done < ID

    log "Stage 2 complete."
}

# ---------------------------------------------------------------------------
# Stage 3 — Pre-trimming QC
# ---------------------------------------------------------------------------
stage3_fastqc_pre() {
    log "Stage 3: Pre-trimming QC"
    mkdir -p "${FASTQC_PRE_DIR}"

    fastqc ./*_merged.fastq.gz -t "${NTHREADS}" -q -o "${FASTQC_PRE_DIR}"
    multiqc "${FASTQC_PRE_DIR}" -o "${FASTQC_PRE_DIR}" \
            --title 'RNA-Seq QC (pre-trimming)'

    log "Stage 3 complete."
}

# ---------------------------------------------------------------------------
# Stage 4 — Adapter & quality trimming (Trim Galore)
# ---------------------------------------------------------------------------
stage4_trim() {
    log "Stage 4: Trimming reads"
    mkdir -p "${TRIMMED_DIR}"

    find ./*_R1_001_merged.fastq.gz > l1
    find ./*_R2_001_merged.fastq.gz > l2

    paste l1 l2 | while IFS=$'\t' read -r R1 R2; do
        log "  Trimming: ${R1}  ${R2}"
        trim_galore "${R1}" "${R2}" \
            --paired \
            -j "${NTHREADS}" \
            --length 15 \
            -o "${TRIMMED_DIR}"
    done

    cd "${TRIMMED_DIR}"
    rename 's/_val_1\.fq\.gz/.fastq.gz/' ./*.fq.gz
    rename 's/_val_2\.fq\.gz/.fastq.gz/' ./*.fq.gz
    rm -f ./*_report.txt
    cd "${HOME_DIR}"

    log "Stage 4 complete."
}

# ---------------------------------------------------------------------------
# Stage 5 — Post-trimming QC
# ---------------------------------------------------------------------------
stage5_fastqc_post() {
    log "Stage 5: Post-trimming QC"
    mkdir -p "${FASTQC_POST_DIR}"

    fastqc "${TRIMMED_DIR}"/*_merged.fastq.gz \
        -t "${NTHREADS}" -q -o "${FASTQC_POST_DIR}"
    multiqc "${FASTQC_POST_DIR}" -o "${FASTQC_POST_DIR}" \
            --title 'RNA-Seq QC (post-trimming)'

    log "Stage 5 complete."
}

# ---------------------------------------------------------------------------
# Stage 6 — First-pass STAR (splice-junction discovery)
# ---------------------------------------------------------------------------
stage6_star_first_pass() {
    log "Stage 6: First-pass STAR alignment"
    mkdir -p "${STAR_FIRST_DIR}"

    find "${TRIMMED_DIR}"/*_R1_001_merged.fastq.gz > l1
    find "${TRIMMED_DIR}"/*_R2_001_merged.fastq.gz > l2

    paste l1 l2 | while IFS=$'\t' read -r R1 R2; do
        outdir=$(echo "${R2}" | cut -f4 -d"/" | cut -f1 -d"_")
        log "  Aligning (pass 1): ${outdir}"

        STAR \
            --runThreadN "${NTHREADS}" \
            --genomeDir "${STAR_INDEX_DIR}" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --outSAMtype None \
            --sjdbGTFfile "${ANNOTATION_DIR}"/* \
            --limitOutSJcollapsed 5000000 \
            --outFileNamePrefix "${STAR_FIRST_DIR}/${outdir}_"
    done

    # Merge and filter high-confidence novel splice junctions
    cat "${STAR_FIRST_DIR}"/*.tab \
        | awk '($5 > 0 && $7 > 2 && $6 == 0)' \
        | cut -f1-6 \
        | sort -u \
        > "${STAR_FIRST_DIR}/SJ_out_total.tab"

    log "  Novel junctions retained: $(wc -l < "${STAR_FIRST_DIR}/SJ_out_total.tab")"

    log "  Building splice-junction-aware index..."
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${SJ_INDEX_DIR}" \
        --genomeFastaFiles "${GENOME_DIR}"/* \
        --sjdbGTFfile "${ANNOTATION_DIR}"/* \
        --sjdbFileChrStartEnd "${STAR_FIRST_DIR}/SJ_out_total.tab"

    log "Stage 6 complete."
}

# ---------------------------------------------------------------------------
# Stage 7 — Second-pass STAR (final BAMs + GeneCounts + chimeric reads)
# ---------------------------------------------------------------------------
stage7_star_final() {
    log "Stage 7: Final STAR alignment"
    mkdir -p "${STAR_FINAL_DIR}"

    paste l1 l2 | while IFS=$'\t' read -r R1 R2; do
        outdir=$(echo "${R2}" | cut -f4 -d"/" | cut -f1 -d"_")
        log "  Aligning (pass 2): ${outdir}"

        STAR \
            --runThreadN "${NTHREADS}" \
            --genomeDir "${SJ_INDEX_DIR}" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbGTFfile "${ANNOTATION_DIR}"/* \
            --quantMode GeneCounts \
            \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 8 \
            --chimOutJunctionFormat 1 \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG -4 \
            --chimMultimapNmax 20 \
            --chimNonchimScoreDropMin 10 \
            --chimMainSegmentMultNmax 1 \
            \
            --alignMatesGapMax 1000000 \
            --alignIntronMax 1000000 \
            --alignSJDBoverhangMin 10 \
            --alignEndsProtrude 20 ConcordantPair \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 100 \
            --outFilterScoreMinOverLread 0.25 \
            --outFilterMatchNminOverLread 0.25 \
            --limitOutSJcollapsed 50000000 \
            \
            --outFileNamePrefix "${STAR_FINAL_DIR}/${outdir}_"

        log "  Done: ${outdir}"
    done

    log "Stage 7 complete."
}

# ---------------------------------------------------------------------------
# Main — activate stages as needed
# ---------------------------------------------------------------------------
main() {
    require_tool STAR
    require_tool fastqc
    require_tool multiqc
    require_tool trim_galore

    # Uncomment the stages you need to run:
    # stage1_build_index
    # stage2_merge_lanes
    # stage3_fastqc_pre
    # stage4_trim
    # stage5_fastqc_post
    # stage6_star_first_pass
    stage7_star_final

    log "Pipeline finished successfully."
}

main "$@"
