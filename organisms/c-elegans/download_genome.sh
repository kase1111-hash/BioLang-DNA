#!/usr/bin/env bash
# download_genome.sh — Download C. elegans WBcel235 reference genome from NCBI
#
# Cryptanalytic analogy: This is the acquisition of the ciphertext.
# Before any frequency analysis or pattern matching can begin, we need
# the raw encoded material in all its forms: the full sequence (genome),
# the known annotations (partial decryption), the protein outputs
# (known plaintext), and the coding regions isolated for frequency analysis.
#
# Usage: ./download_genome.sh [--no-decompress]
# Run from the repository root or organisms/c-elegans/ directory.

set -euo pipefail

# --- Configuration ---
ACCESSION="GCF_000002985.6"
ASSEMBLY="WBcel235"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/${ACCESSION}_${ASSEMBLY}"
DATASETS_API="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${ACCESSION}/download"

# Resolve data directory relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"

DECOMPRESS=true
if [[ "${1:-}" == "--no-decompress" ]]; then
    DECOMPRESS=false
fi

# --- Functions ---
log() { echo "[$(date '+%H:%M:%S')] $*"; }

download_with_retry() {
    local url="$1"
    local output="$2"
    local max_retries=4
    local delay=2

    for attempt in $(seq 1 $max_retries); do
        log "Downloading $(basename "$output") (attempt $attempt/$max_retries)..."
        if wget -q --show-progress -O "$output" "$url" 2>&1; then
            log "  ✓ Downloaded $(basename "$output") ($(du -h "$output" | cut -f1))"
            return 0
        fi
        if [[ $attempt -lt $max_retries ]]; then
            log "  Retrying in ${delay}s..."
            sleep $delay
            delay=$((delay * 2))
        fi
    done

    log "  ✗ Failed to download $(basename "$output") after $max_retries attempts"
    return 1
}

# --- Main ---
log "=== C. elegans Genome Acquisition ==="
log "Assembly: ${ASSEMBLY} (${ACCESSION})"
log "Target:   ${DATA_DIR}"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

# Method 1: Try NCBI Datasets API (single zip with everything)
log ""
log "--- Method 1: NCBI Datasets API ---"
if command -v datasets &>/dev/null; then
    log "Using NCBI datasets CLI..."
    datasets download genome accession "$ACCESSION" \
        --include genome,gff3,protein,cds \
        --filename ncbi_dataset.zip && \
    unzip -o ncbi_dataset.zip -d ncbi_dataset && \
    NCBI_DIR="ncbi_dataset/ncbi_dataset/data/${ACCESSION}"
    if [[ -d "$NCBI_DIR" ]]; then
        cp "${NCBI_DIR}/"*_genomic.fna genome.fna 2>/dev/null || true
        cp "${NCBI_DIR}/"*.gff annotations.gff3 2>/dev/null || true
        cp "${NCBI_DIR}/"*protein.faa proteome.faa 2>/dev/null || true
        cp "${NCBI_DIR}/"*cds_from_genomic.fna cds.fna 2>/dev/null || true
        rm -rf ncbi_dataset ncbi_dataset.zip
        log "✓ All files acquired via datasets CLI"
        exit 0
    fi
fi

# Method 2: Direct FTP download (individual files)
log ""
log "--- Method 2: Direct FTP download ---"

FAILED=0

# 1. Full genome FASTA
download_with_retry \
    "${BASE_URL}/${ACCESSION}_${ASSEMBLY}_genomic.fna.gz" \
    "genome.fna.gz" || FAILED=$((FAILED + 1))

# 2. Gene annotations (GFF3)
download_with_retry \
    "${BASE_URL}/${ACCESSION}_${ASSEMBLY}_genomic.gff.gz" \
    "annotations.gff3.gz" || FAILED=$((FAILED + 1))

# 3. Predicted proteome
download_with_retry \
    "${BASE_URL}/${ACCESSION}_${ASSEMBLY}_protein.faa.gz" \
    "proteome.faa.gz" || FAILED=$((FAILED + 1))

# 4. Coding sequences (CDS) — isolated for frequency analysis
download_with_retry \
    "${BASE_URL}/${ACCESSION}_${ASSEMBLY}_cds_from_genomic.fna.gz" \
    "cds.fna.gz" || FAILED=$((FAILED + 1))

if [[ $FAILED -gt 0 ]]; then
    log ""
    log "⚠ $FAILED download(s) failed. You may need to download manually from:"
    log "  ${BASE_URL}/"
    log ""
    log "Or use the NCBI Datasets API:"
    log "  ${DATASETS_API}?include_annotation_type=GENOME_FASTA,GENOME_GFF,PROT_FASTA,CDS_FASTA"
    exit 1
fi

# Decompress
if $DECOMPRESS; then
    log ""
    log "--- Decompressing ---"
    for gz in genome.fna.gz annotations.gff3.gz proteome.faa.gz cds.fna.gz; do
        if [[ -f "$gz" ]]; then
            gunzip -f "$gz"
            log "  ✓ Decompressed $gz"
        fi
    done
fi

# Verify
log ""
log "--- Verification ---"
for f in genome.fna annotations.gff3 proteome.faa cds.fna; do
    if [[ -f "$f" ]]; then
        SIZE=$(du -h "$f" | cut -f1)
        LINES=$(wc -l < "$f")
        log "  ✓ $f — ${SIZE}, ${LINES} lines"
    elif [[ -f "${f}.gz" ]]; then
        SIZE=$(du -h "${f}.gz" | cut -f1)
        log "  ✓ ${f}.gz — ${SIZE} (compressed)"
    else
        log "  ✗ $f — MISSING"
    fi
done

# Quick stats
log ""
log "--- Quick Stats ---"
if [[ -f genome.fna ]]; then
    CHROMS=$(grep -c "^>" genome.fna)
    BASES=$(grep -v "^>" genome.fna | tr -d '\n' | wc -c)
    log "  Genome: $CHROMS sequences, ~$((BASES / 1000000))M bases"
fi
if [[ -f cds.fna ]]; then
    CDS_COUNT=$(grep -c "^>" cds.fna)
    log "  CDS: $CDS_COUNT coding sequences"
fi
if [[ -f proteome.faa ]]; then
    PROT_COUNT=$(grep -c "^>" proteome.faa)
    log "  Proteome: $PROT_COUNT proteins"
fi
if [[ -f annotations.gff3 ]]; then
    GENE_COUNT=$(grep -c $'\tgene\t' annotations.gff3 || echo "0")
    log "  Annotations: $GENE_COUNT genes"
fi

log ""
log "=== Acquisition complete ==="
