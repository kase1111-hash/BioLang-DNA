#!/usr/bin/env python3
"""
generate_test_data.py — Generate synthetic C. elegans-like test data.

Creates a small but realistic genome, GFF3, and CDS FASTA for testing
the frequency_analyzer.py pipeline when real genome data is unavailable.

The synthetic data mimics key properties of C. elegans:
- 6 chromosomes (I-V + X), scaled down to ~100kb each
- Realistic gene density with CDS, exon, and gene features
- C. elegans-like codon bias (AT-rich organism)
- Intergenic, regulatory, and coding regions

Usage:
    python generate_test_data.py --outdir organisms/c-elegans/data
"""

import argparse
import os
import random


# C. elegans has a slight AT bias; weight nucleotide probabilities accordingly
BASE_WEIGHTS = {"A": 0.32, "T": 0.32, "C": 0.18, "G": 0.18}

CHROMOSOMES = ["I", "II", "III", "IV", "V", "X"]
CHR_SIZE = 100_000  # 100kb per chromosome (scaled down from ~15-21 Mb)

# Simplified codon table for generating biased CDS
PREFERRED_CODONS = {
    "F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC"],
    "I": ["ATT", "ATC"], "M": ["ATG"], "V": ["GTT", "GTC", "GTA"],
    "S": ["TCT", "TCC", "TCA", "AGT", "AGC"],
    "P": ["CCT", "CCA"], "T": ["ACT", "ACC", "ACA"],
    "A": ["GCT", "GCC", "GCA"], "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"], "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"], "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"], "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"], "W": ["TGG"],
    "R": ["AGA", "AGG", "CGT", "CGC"],
    "G": ["GGA", "GGT", "GGC"],
}
AMINO_ACIDS = list(PREFERRED_CODONS.keys())
STOP_CODONS = ["TAA", "TAG", "TGA"]


def weighted_base():
    r = random.random()
    cumulative = 0
    for base, weight in BASE_WEIGHTS.items():
        cumulative += weight
        if r <= cumulative:
            return base
    return "T"


def random_sequence(length):
    return "".join(weighted_base() for _ in range(length))


def random_cds(n_codons):
    """Generate a random CDS with realistic codon bias."""
    seq = "ATG"  # start codon
    for _ in range(n_codons - 2):
        aa = random.choice(AMINO_ACIDS)
        codon = random.choice(PREFERRED_CODONS[aa])
        seq += codon
    seq += random.choice(STOP_CODONS)
    return seq


def generate_test_data(outdir):
    os.makedirs(outdir, exist_ok=True)
    random.seed(42)

    genome = {}
    genes = []
    cds_records = []
    gene_counter = 0

    for chrom in CHROMOSOMES:
        seqid = f"NC_00{CHROMOSOMES.index(chrom)+1}.1"
        seq = list(random_sequence(CHR_SIZE))

        # Place genes along the chromosome
        pos = 3000  # leave room for first regulatory region
        while pos < CHR_SIZE - 5000:
            gene_counter += 1
            gene_name = f"gene_{gene_counter:04d}"
            strand = random.choice(["+", "-"])

            # Gene with 2-5 exons
            n_exons = random.randint(2, 5)
            exon_starts = []
            exon_ends = []
            current = pos

            for _ in range(n_exons):
                exon_len = random.randint(100, 500)
                exon_starts.append(current)
                exon_ends.append(current + exon_len)
                current += exon_len + random.randint(50, 300)  # intron

            gene_start = exon_starts[0]
            gene_end = exon_ends[-1]

            # Generate CDS and embed in genome
            total_cds_len = sum(e - s for s, e in zip(exon_starts, exon_ends))
            n_codons = total_cds_len // 3
            if n_codons < 10:
                pos = gene_end + random.randint(1000, 3000)
                continue

            cds_seq = random_cds(n_codons)
            cds_offset = 0

            for es, ee in zip(exon_starts, exon_ends):
                exon_len = ee - es
                chunk = cds_seq[cds_offset:cds_offset + exon_len]
                for j, base in enumerate(chunk):
                    if es + j < CHR_SIZE:
                        seq[es + j] = base
                cds_offset += exon_len

            # Record gene features
            genes.append({
                "seqid": seqid,
                "gene_name": gene_name,
                "gene_start": gene_start + 1,  # 1-based
                "gene_end": gene_end,
                "strand": strand,
                "exon_starts": [s + 1 for s in exon_starts],
                "exon_ends": exon_ends,
            })

            # CDS record
            full_cds = cds_seq[:cds_offset]
            cds_records.append((gene_name, seqid, full_cds))

            # Add some tandem repeats in intergenic regions
            if random.random() < 0.3:
                repeat_pos = gene_end + random.randint(100, 500)
                unit = random_sequence(random.randint(2, 8))
                copies = random.randint(5, 20)
                repeat_seq = unit * copies
                for j, base in enumerate(repeat_seq):
                    if repeat_pos + j < CHR_SIZE:
                        seq[repeat_pos + j] = base

            pos = gene_end + random.randint(1000, 3000)

        genome[seqid] = "".join(seq)

    # Write genome FASTA
    genome_path = os.path.join(outdir, "genome.fna")
    with open(genome_path, "w") as fh:
        for seqid, seq in genome.items():
            fh.write(f">{seqid} synthetic C. elegans chromosome\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    print(f"Wrote genome: {genome_path} ({len(genome)} chromosomes, {sum(len(s) for s in genome.values()):,} bp)")

    # Write GFF3
    gff_path = os.path.join(outdir, "annotations.gff3")
    with open(gff_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for g in genes:
            fh.write(f"{g['seqid']}\tsynthetic\tgene\t{g['gene_start']}\t{g['gene_end']}\t.\t{g['strand']}\t.\tID={g['gene_name']};Name={g['gene_name']}\n")
            for i, (es, ee) in enumerate(zip(g["exon_starts"], g["exon_ends"])):
                fh.write(f"{g['seqid']}\tsynthetic\texon\t{es}\t{ee}\t.\t{g['strand']}\t.\tParent={g['gene_name']}\n")
                fh.write(f"{g['seqid']}\tsynthetic\tCDS\t{es}\t{ee}\t.\t{g['strand']}\t0\tParent={g['gene_name']}\n")
    print(f"Wrote annotations: {gff_path} ({len(genes)} genes)")

    # Write CDS FASTA
    cds_path = os.path.join(outdir, "cds.fna")
    with open(cds_path, "w") as fh:
        for name, seqid, seq in cds_records:
            fh.write(f">{name} [locus_tag={name}] [chromosome={seqid}]\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    print(f"Wrote CDS: {cds_path} ({len(cds_records)} sequences)")

    return genome_path, gff_path, cds_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic test data for frequency_analyzer.py")
    parser.add_argument("--outdir", default="organisms/c-elegans/data", help="Output directory")
    args = parser.parse_args()
    generate_test_data(args.outdir)
