#!/bin/bash

set -euo pipefail

readonly VERSION="v1.0.0"
readonly HELP_TEXT="Usage: $0 -tumor <tumor_bam_cram> -normal <normal_bam_cram> -ref_fasta <ref_fasta> -variants_vcf <variants_to_prep> -stm_repo <sequence_tube_map_repo>"

if [ $# -eq 0 ]; then
  echo "$0 $VERSION"
  echo "$HELP_TEXT"
  exit 1
fi

while [[ $# -gt 1 ]]; do
  key=${1}
  case ${key} in
    -tumor)
      TUMOR=${2};;
    -normal)
      NORMAL=${2};;
    -ref_fasta)
      REF_FASTA=${2};;
    -variants_vcf)
      VARIANTS_VCF=${2};;
    -stm_repo)
      STM_REPO=${2};;
    *)
      echo "$HELP_TEXT"
      exit 65
    ;;
  esac
  shift
  shift
done

function log() {
  echo "[$(date +'%Y%m%dT%H%M%S')] :: ${*}" >&2
}

function ensure_exists() {
  hash "$1" &> /dev/null
  if [ $? -eq 1 ]; then
    echo "$1 not found, but is necessary to continue"
    exit 1
  fi
}

ensure_exists "Lancet2"
ensure_exists "samtools"
ensure_exists "vg"
ensure_exists "bcftools"
ensure_exists "jq"


TUMOR="$(realpath "${TUMOR}")"
NORMAL="$(realpath "${NORMAL}")"
REF_FASTA=$(realpath "${REF_FASTA}")
STM_REPO=$(realpath "${STM_REPO}")
OUT_DIR="${STM_REPO}/exampleData"

bcftools query -f '%CHROM\t%REF\t%ALT\t%POS\t%INFO/TYPE\n' "${VARIANTS_VCF}" | while IFS= read -r line
do
  CHROM="$(echo "${line}" | cut -f1)"
  REF="$(echo "${line}" | cut -f2)"
  ALT="$(echo "${line}" | cut -f3)"
  VAR_POS="$(echo "${line}" | cut -f4)"
  VAR_TYPE="$(echo "${line}" | cut -f5)"
  WIN_START=$((VAR_POS-500))
  WIN_END=$((VAR_POS+500))

  REGION="${CHROM}:${WIN_START}-${WIN_END}"
  CURR_NAME="${CHROM}_${WIN_START}_${WIN_END}_${REF}_${ALT}_${VAR_TYPE}"
  WORK_DIR="${OUT_DIR}/workdir/${CURR_NAME}"
  CHUNK_DIR="${OUT_DIR}/chunks/${CURR_NAME}"
  rm -rf "${WORK_DIR}" && mkdir -p "${WORK_DIR}"
  rm -rf "${CHUNK_DIR}" && mkdir -p "${CHUNK_DIR}"

  log "Running Lancet2 to serialize window graph for ${REGION}"
  Lancet2 pipeline --window-size 1000 --padding 0 --tumor "${TUMOR}" --normal "${NORMAL}" --reference "${REF_FASTA}" \
    --region "${CHROM}:${VAR_POS}" --out-vcfgz "${WORK_DIR}/${CURR_NAME}.lancet_calls.vcf.gz" --graphs-dir "${WORK_DIR}/graphs"

  TMR_FQ="${WORK_DIR}/tumor_reads.fastq"
  log "Extracting tumor reads from ${REGION}"
  samtools view --uncompressed "${TUMOR}" "${REGION}" -o tmp.bam && \
  samtools collate -Ou tmp.bam | samtools fastq >| "${TMR_FQ}" && rm -f tmp.bam

  NML_FQ="${WORK_DIR}/normal_reads.fastq"
  log "Extracting normal reads from ${REGION}"
  samtools view --uncompressed "${NORMAL}" "${REGION}" -o tmp.bam && \
  samtools collate -Ou tmp.bam | samtools fastq >| "${NML_FQ}" && rm -f tmp.bam

  UNCHOPPED_GFA="${WORK_DIR}/${CURR_NAME}.unchopped.gfa"
  log "Unchopping Lancet GFA graph and creating giraffe indexes for ${REGION}"
  vg mod --unchop "$(ls "${WORK_DIR}/graphs/poa_graph/*.gfa")" | sed 's/ref0/'"${CHROM}"'/g' >| "${UNCHOPPED_GFA}" && \
  vg autoindex --workflow giraffe --gfa "${UNCHOPPED_GFA}" --prefix "${WORK_DIR}/${CHROM}_${WIN_START}_${WIN_END}"

  GBZ="${WORK_DIR}/${CHROM}_${WIN_START}_${WIN_END}.giraffe.gbz"
  MIN="${WORK_DIR}/${CHROM}_${WIN_START}_${WIN_END}.min"
  DIST="${WORK_DIR}/${CHROM}_${WIN_START}_${WIN_END}.dist"

  NML_NAME="$(bcftools view --header-only "${WORK_DIR}/${CURR_NAME}.lancet_calls.vcf.gz" | tail -n1 | cut -f10)"
  TMR_NAME="$(bcftools view --header-only "${WORK_DIR}/${CURR_NAME}.lancet_calls.vcf.gz" | tail -n1 | cut -f11)"

  log "Running giraffe to align tumor reads to the ${REGION} graph"
  vg giraffe --gbz-name "${GBZ}" --minimizer-name "${MIN}" --dist-name "${DIST}" --fastq-in "${TMR_FQ}" --sample "${TMR_NAME}" >| tmp.gam && \
  vg gamsort --index "${WORK_DIR}/tumor.gam.gai" tmp.gam >| "${WORK_DIR}/tumor.gam" && rm -f tmp.gam

  log "Running giraffe to align normal reads to the ${REGION} graph"
  vg giraffe --gbz-name "${GBZ}" --minimizer-name "${MIN}" --dist-name "${DIST}" --fastq-in "${NML_FQ}" --sample "${NML_NAME}" >| tmp.gam && \
  vg gamsort --index "${WORK_DIR}/normal.gam.gai" tmp.gam >| "${WORK_DIR}/normal.gam" && rm -f tmp.gam

  log "Creating local VG chunk for ${REGION} to visualize in SequenceTubeMap"
  "${STM_REPO}/prep_local_chunk.sh" -x "${GBZ}" -r "${REGION}" -o "${CHUNK_DIR}" \
    -g "${WORK_DIR}/normal.gam" -p '{"mainPalette": "blues", "auxPalette": "blues"}' \
    -g "${WORK_DIR}/tumor.gam" -p '{"mainPalette": "reds", "auxPalette": "reds"}' \
    -d "${REF} -> ${ALT} variant (${VAR_TYPE}) at ${CHROM}:${VAR_POS}" >> "${OUT_DIR}/index.bed"

  vg paths --drop-paths --paths-by "path_cover_" --xg "${CHUNK_DIR}/chunk.vg" >"${CHUNK_DIR}/chunk.vg.new" && \
  mv "${CHUNK_DIR}/chunk.vg.new" "${CHUNK_DIR}/chunk.vg"
done
