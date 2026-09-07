#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  03_run_dcc_typing.sh --sample NAME --r1 READ1.fastq.gz --r2 READ2.fastq.gz \
    --subspecies ABS|MAS --varscan-jar VarScan.v2.3.9.jar [--threads 8] [--outdir DIR]
EOF
}

SAMPLE=""
R1=""
R2=""
SUBSPECIES=""
VARSCAN_JAR=""
THREADS=8
OUTDIR="."

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2 ;;
    --r1) R1="$2"; shift 2 ;;
    --r2) R2="$2"; shift 2 ;;
    --subspecies) SUBSPECIES="$2"; shift 2 ;;
    --varscan-jar) VARSCAN_JAR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -z "$SAMPLE" || -z "$R1" || -z "$R2" || -z "$SUBSPECIES" || -z "$VARSCAN_JAR" ]]; then
  usage >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

case "$SUBSPECIES" in
  ABS)
    REF="${PACKAGE_DIR}/references/ABS_GZ002_NZ_CP034181.1.fna"
    MARKERS="${PACKAGE_DIR}/markers/ABS_DCC_markers.tsv"
    ;;
  MAS)
    REF="${PACKAGE_DIR}/references/MAS_CCUG48898_NZ_AP014547.1.fna"
    MARKERS="${PACKAGE_DIR}/markers/MAS_DCC_markers.tsv"
    ;;
  *)
    echo "--subspecies must be ABS or MAS" >&2
    exit 1
    ;;
esac

for command in bwa samtools java python3; do
  command -v "$command" >/dev/null || { echo "Missing dependency: $command" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
PREFIX="${OUTDIR}/${SAMPLE}"

[[ -f "${REF}.bwt" ]] || bwa index "$REF"
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

bwa mem -t "$THREADS" -c 100 -M \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -o "${PREFIX}.sorted.bam" -

samtools index "${PREFIX}.sorted.bam"

samtools mpileup -q 30 -Q 20 -B -O -f "$REF" "${PREFIX}.sorted.bam" \
  > "${PREFIX}.pileup"

java -jar "$VARSCAN_JAR" mpileup2cns "${PREFIX}.pileup" \
  --min-coverage 3 \
  --min-avg-qual 20 \
  --min-var-freq 0.75 \
  --min-reads2 2 \
  --strand-filter 0 \
  > "${PREFIX}.cns"

python3 "${SCRIPT_DIR}/01_cns_to_typesnp.py" \
  "${PREFIX}.cns" -o "${PREFIX}.typesnp.tsv"

python3 "${SCRIPT_DIR}/02_assign_dcc.py" \
  "$MARKERS" "${PREFIX}.typesnp.tsv" \
  --subspecies "$SUBSPECIES" --sample "$SAMPLE" \
  > "${PREFIX}.dcc_assignment.tsv"

printf 'DCC assignment written to %s\n' "${PREFIX}.dcc_assignment.tsv"
