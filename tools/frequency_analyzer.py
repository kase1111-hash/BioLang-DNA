#!/usr/bin/env python3
"""
frequency_analyzer.py — Statistical Fingerprinting Pipeline for GeneLang

Cryptanalytic Analogy:
    This tool performs the foundational frequency analysis on a genome, analogous
    to the first step in breaking any cipher: counting character frequencies.
    Just as letter-frequency analysis reveals the structure of substitution
    ciphers, codon and k-mer frequency analysis reveals the "dialect" and
    structural preferences encoded in a genome's regulatory logic.

    The genome is treated as ciphertext. By measuring statistical properties
    across different functional regions (coding, regulatory, intergenic), we
    build the baseline fingerprint needed for all downstream pattern recognition.

Sub-analyses:
    3a. Codon Frequency Table — Letter frequency analysis for the genetic code
    3b. K-mer Distribution    — N-gram analysis across functional regions
    3c. Information Density   — Entropy mapping to find information-dense zones
    3d. Repeat Element Catalog — Key reuse detection (tandem repeats)

Usage:
    python frequency_analyzer.py codon  --genome FASTA --gff GFF3 --cds CDS_FASTA --outdir DIR
    python frequency_analyzer.py kmer   --genome FASTA --gff GFF3 --outdir DIR
    python frequency_analyzer.py density --genome FASTA --gff GFF3 --outdir DIR
    python frequency_analyzer.py repeats --genome FASTA --gff GFF3 --outdir DIR
    python frequency_analyzer.py all    --genome FASTA --gff GFF3 --cds CDS_FASTA --outdir DIR
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from itertools import product
from math import log2, sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASES = "ACGT"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]  # all 64 codons

# Standard genetic code: codon -> amino acid
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Group codons by amino acid for RSCU calculation
AA_CODON_FAMILIES = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    AA_CODON_FAMILIES[aa].append(codon)

UPSTREAM_SIZE = 2000  # bp upstream of gene start for regulatory regions


# ---------------------------------------------------------------------------
# FASTA / GFF3 Parsing
# ---------------------------------------------------------------------------

def parse_fasta(filepath):
    """Parse a FASTA file into a dict of {sequence_id: sequence_string}.

    Handles multi-line sequences and strips whitespace. Sequence IDs are
    taken from the first whitespace-delimited token after '>'.
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            elif current_id is not None:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq).upper()

    return sequences


def parse_gff3(filepath):
    """Parse a GFF3 file into a list of feature dicts.

    Returns features with keys: seqid, source, type, start, end, score,
    strand, phase, attributes (as dict).
    """
    features = []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, val = item.split("=", 1)
                    attrs[key] = val

            features.append({
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),  # 1-based inclusive
                "end": int(parts[4]),    # 1-based inclusive
                "score": parts[5],
                "strand": parts[6],
                "phase": parts[7],
                "attributes": attrs,
            })

    return features


# ---------------------------------------------------------------------------
# Region Classification
# ---------------------------------------------------------------------------

def classify_regions(features, genome_lengths):
    """Classify genomic regions into coding, regulatory, and intergenic.

    Coding: CDS features from GFF3.
    Regulatory: 2kb upstream of each gene start (strand-aware).
    Intergenic: everything not in coding or regulatory.

    Returns:
        coding_intervals: dict of seqid -> sorted list of (start, end) 0-based
        regulatory_intervals: dict of seqid -> sorted list of (start, end) 0-based
        intergenic_intervals: dict of seqid -> sorted list of (start, end) 0-based
    """
    coding = defaultdict(list)
    regulatory = defaultdict(list)

    # Collect CDS intervals (convert to 0-based)
    for f in features:
        if f["type"] == "CDS":
            coding[f["seqid"]].append((f["start"] - 1, f["end"]))

    # Collect regulatory intervals (2kb upstream of gene, strand-aware)
    for f in features:
        if f["type"] == "gene":
            seqid = f["seqid"]
            seq_len = genome_lengths.get(seqid, float("inf"))
            if f["strand"] == "+":
                reg_start = max(0, f["start"] - 1 - UPSTREAM_SIZE)
                reg_end = f["start"] - 1  # up to gene start
            else:
                reg_start = f["end"]  # past gene end
                reg_end = min(seq_len, f["end"] + UPSTREAM_SIZE)
            if reg_start < reg_end:
                regulatory[seqid].append((reg_start, reg_end))

    # Merge overlapping intervals for coding and regulatory
    coding_merged = {s: _merge_intervals(ivs) for s, ivs in coding.items()}
    regulatory_merged = {s: _merge_intervals(ivs) for s, ivs in regulatory.items()}

    # Compute intergenic as complement
    intergenic = {}
    for seqid, seq_len in genome_lengths.items():
        occupied = _merge_intervals(
            coding_merged.get(seqid, []) + regulatory_merged.get(seqid, [])
        )
        intergenic[seqid] = _complement_intervals(occupied, seq_len)

    return coding_merged, regulatory_merged, intergenic


def _merge_intervals(intervals):
    """Merge overlapping intervals. Input: list of (start, end). Returns sorted merged list."""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [sorted_ivs[0]]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _complement_intervals(intervals, total_length):
    """Compute the complement of merged intervals within [0, total_length)."""
    complement = []
    prev_end = 0
    for start, end in intervals:
        if start > prev_end:
            complement.append((prev_end, start))
        prev_end = max(prev_end, end)
    if prev_end < total_length:
        complement.append((prev_end, total_length))
    return complement


def extract_sequences_for_regions(genome, intervals):
    """Extract and concatenate sequences for a set of intervals per chromosome.

    Args:
        genome: dict of seqid -> sequence string
        intervals: dict of seqid -> list of (start, end) 0-based

    Returns:
        list of sequence strings (one per interval)
    """
    seqs = []
    for seqid, ivs in intervals.items():
        seq = genome.get(seqid, "")
        for start, end in ivs:
            chunk = seq[start:end]
            if chunk:
                seqs.append(chunk)
    return seqs


# ---------------------------------------------------------------------------
# 3a. Codon Frequency Analysis
# ---------------------------------------------------------------------------

def analyze_codons(cds_fasta_path, outdir):
    """Task 3a: Codon Frequency Table.

    Cryptanalytic Analogy:
        Codon frequency analysis is analogous to letter frequency analysis in
        substitution cipher cracking. Just as English uses 'E' ~12.7% and 'Z'
        ~0.07%, organisms prefer certain codons over synonymous alternatives.
        These preferences reveal the "dialect" of the organism's genetic code
        and hint at translational selection pressures — optimization choices
        the genome has made over evolutionary time.

    Computes:
        - Raw codon counts across all CDS regions
        - Relative Synonymous Codon Usage (RSCU) per codon
        - Chi-squared test for codon usage bias per amino acid family
        - CSV output + summary statistics + visualization
    """
    print("\n=== Task 3a: Codon Frequency Analysis ===")
    print("Analogy: Letter frequency analysis for the genetic substitution cipher\n")

    os.makedirs(outdir, exist_ok=True)

    # Parse CDS sequences
    cds_seqs = parse_fasta(cds_fasta_path)
    print(f"  Loaded {len(cds_seqs)} CDS sequences")

    # Count codons across all CDS
    codon_counts = Counter()
    total_codons = 0
    skipped_seqs = 0

    for seq_id, seq in cds_seqs.items():
        # Clean sequence: only standard bases
        seq = re.sub(r"[^ACGT]", "", seq.upper())
        if len(seq) < 3:
            skipped_seqs += 1
            continue
        # Count codons (non-overlapping, reading frame from start)
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if len(codon) == 3 and all(b in BASES for b in codon):
                codon_counts[codon] += 1
                total_codons += 1

    if skipped_seqs:
        print(f"  Skipped {skipped_seqs} sequences (too short or invalid)")
    print(f"  Total codons counted: {total_codons:,}")

    # Compute RSCU (Relative Synonymous Codon Usage)
    # RSCU_i = (observed_i * n_synonyms) / sum(observed for all synonyms)
    rscu = {}
    for aa, family_codons in AA_CODON_FAMILIES.items():
        family_total = sum(codon_counts.get(c, 0) for c in family_codons)
        n_syn = len(family_codons)
        for codon in family_codons:
            if family_total > 0:
                rscu[codon] = (codon_counts.get(codon, 0) * n_syn) / family_total
            else:
                rscu[codon] = 0.0

    # Chi-squared test for each amino acid family
    chi_squared_results = {}
    for aa, family_codons in AA_CODON_FAMILIES.items():
        if len(family_codons) <= 1:
            continue  # Met, Trp, stop — no synonyms to compare
        observed = [codon_counts.get(c, 0) for c in family_codons]
        total = sum(observed)
        if total == 0:
            continue
        expected = [total / len(family_codons)] * len(family_codons)
        chi2, p_value = stats.chisquare(observed, f_exp=expected)
        chi_squared_results[aa] = {
            "chi2": chi2,
            "p_value": p_value,
            "df": len(family_codons) - 1,
            "codons": family_codons,
            "observed": observed,
            "expected": expected,
        }

    # --- Output CSV ---
    csv_path = os.path.join(outdir, "codon_frequency.csv")
    rows = []
    for codon in CODONS:
        aa = CODON_TABLE[codon]
        count = codon_counts.get(codon, 0)
        freq = count / total_codons if total_codons > 0 else 0
        rows.append({
            "codon": codon,
            "amino_acid": aa,
            "count": count,
            "frequency": round(freq, 6),
            "rscu": round(rscu.get(codon, 0), 4),
        })

    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False)
    print(f"  Wrote codon frequency table: {csv_path}")

    # --- Chi-squared summary CSV ---
    chi2_path = os.path.join(outdir, "codon_chi_squared.csv")
    chi2_rows = []
    for aa, result in sorted(chi_squared_results.items()):
        chi2_rows.append({
            "amino_acid": aa,
            "codons": "/".join(result["codons"]),
            "chi_squared": round(result["chi2"], 4),
            "p_value": f"{result['p_value']:.2e}",
            "df": result["df"],
            "significant": "Yes" if result["p_value"] < 0.05 else "No",
        })
    pd.DataFrame(chi2_rows).to_csv(chi2_path, index=False)
    print(f"  Wrote chi-squared results: {chi2_path}")

    # --- Summary stats ---
    summary_path = os.path.join(outdir, "codon_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("=== Codon Frequency Analysis Summary ===\n\n")
        fh.write(f"Total CDS sequences analyzed: {len(cds_seqs)}\n")
        fh.write(f"Total codons counted: {total_codons:,}\n\n")

        fh.write("Top 10 most frequent codons:\n")
        for codon, count in codon_counts.most_common(10):
            aa = CODON_TABLE[codon]
            freq = count / total_codons * 100
            fh.write(f"  {codon} ({aa}): {count:>10,}  ({freq:.2f}%)  RSCU={rscu[codon]:.3f}\n")

        fh.write("\nTop 10 least frequent codons:\n")
        for codon, count in codon_counts.most_common()[-10:]:
            aa = CODON_TABLE[codon]
            freq = count / total_codons * 100
            fh.write(f"  {codon} ({aa}): {count:>10,}  ({freq:.2f}%)  RSCU={rscu[codon]:.3f}\n")

        sig_count = sum(1 for r in chi_squared_results.values() if r["p_value"] < 0.05)
        fh.write(f"\nChi-squared test: {sig_count}/{len(chi_squared_results)} amino acid families "
                 f"show significant codon usage bias (p < 0.05)\n")
    print(f"  Wrote summary: {summary_path}")

    # --- Visualization ---
    _plot_codon_frequencies(df, rscu, outdir)

    print("  Task 3a complete.\n")
    return df


def _plot_codon_frequencies(df, rscu, outdir):
    """Generate codon frequency visualizations."""
    # 1. Bar chart of all 64 codon frequencies, colored by amino acid
    fig, ax = plt.subplots(figsize=(18, 6))
    colors = plt.cm.tab20(np.linspace(0, 1, 20))
    aa_list = sorted(set(CODON_TABLE.values()))
    aa_color_map = {aa: colors[i % 20] for i, aa in enumerate(aa_list)}

    bar_colors = [aa_color_map[CODON_TABLE[c]] for c in df["codon"]]
    ax.bar(range(len(df)), df["frequency"], color=bar_colors, edgecolor="none")
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["codon"], rotation=90, fontsize=6)
    ax.set_xlabel("Codon")
    ax.set_ylabel("Frequency")
    ax.set_title("C. elegans Codon Usage Frequency (All 64 Codons)")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "codon_frequency_bar.png"), dpi=150)
    plt.close()

    # 2. RSCU heatmap grouped by amino acid
    fig, ax = plt.subplots(figsize=(14, 8))
    rscu_data = []
    labels = []
    for aa in sorted(AA_CODON_FAMILIES.keys()):
        for codon in sorted(AA_CODON_FAMILIES[aa]):
            rscu_data.append(rscu.get(codon, 0))
            labels.append(f"{codon} ({aa})")

    # Reshape for heatmap: 8 rows x 8 cols
    n = len(rscu_data)
    ncols = 8
    nrows = (n + ncols - 1) // ncols
    padded = rscu_data + [0] * (nrows * ncols - n)
    padded_labels = labels + [""] * (nrows * ncols - n)

    matrix = np.array(padded).reshape(nrows, ncols)
    label_matrix = np.array(padded_labels).reshape(nrows, ncols)

    sns.heatmap(matrix, annot=label_matrix, fmt="", cmap="YlOrRd",
                xticklabels=False, yticklabels=False, ax=ax,
                cbar_kws={"label": "RSCU"})
    ax.set_title("Relative Synonymous Codon Usage (RSCU) — C. elegans")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "codon_rscu_heatmap.png"), dpi=150)
    plt.close()

    print("  Wrote visualizations: codon_frequency_bar.png, codon_rscu_heatmap.png")


# ---------------------------------------------------------------------------
# 3b. K-mer Distribution Analysis
# ---------------------------------------------------------------------------

def analyze_kmers(genome, features, genome_lengths, outdir):
    """Task 3b: K-mer Distribution Analysis.

    Cryptanalytic Analogy:
        K-mer analysis is the n-gram frequency attack applied to genomic
        ciphertext. In classical cryptanalysis, digram and trigram frequencies
        reveal word boundaries, common affixes, and syntactic structure.
        Similarly, k-mer frequencies in different genomic regions reveal
        the distinct "vocabularies" used by coding sequences (constrained
        by the genetic code), regulatory regions (shaped by transcription
        factor binding preferences), and intergenic DNA (relatively
        unconstrained — the "noise floor" baseline).

    Computes for k=2 through k=8, separately for coding/regulatory/intergenic:
        - Frequency distribution of all k-mers
        - Shannon entropy (bits per base pair)
        - Kullback-Leibler divergence from uniform distribution
        - Top 20 overrepresented and underrepresented k-mers
    """
    print("\n=== Task 3b: K-mer Distribution Analysis ===")
    print("Analogy: N-gram frequency attack across functional genomic regions\n")

    kmer_outdir = os.path.join(outdir, "kmer")
    os.makedirs(kmer_outdir, exist_ok=True)

    # Classify regions
    print("  Classifying genomic regions...")
    coding_ivs, regulatory_ivs, intergenic_ivs = classify_regions(features, genome_lengths)

    # Extract sequences for each region type
    region_seqs = {
        "coding": extract_sequences_for_regions(genome, coding_ivs),
        "regulatory": extract_sequences_for_regions(genome, regulatory_ivs),
        "intergenic": extract_sequences_for_regions(genome, intergenic_ivs),
    }

    for rname, seqs in region_seqs.items():
        total_bp = sum(len(s) for s in seqs)
        print(f"  {rname}: {len(seqs)} segments, {total_bp:,} bp total")

    # Analyze k-mers for each region type and k value
    entropy_summary = []  # for the comparative heatmap

    for k in range(2, 9):
        n_possible = 4 ** k
        print(f"\n  --- k={k} ({n_possible} possible {k}-mers) ---")

        for region_name, seqs in region_seqs.items():
            # Count k-mers
            kmer_counts = Counter()
            total_kmers = 0
            for seq in seqs:
                clean = re.sub(r"[^ACGT]", "", seq)
                for i in range(len(clean) - k + 1):
                    kmer = clean[i:i + k]
                    kmer_counts[kmer] += 1
                    total_kmers += 1

            if total_kmers == 0:
                print(f"    {region_name}: no valid k-mers")
                continue

            # Frequency distribution
            freqs = {kmer: count / total_kmers for kmer, count in kmer_counts.items()}

            # Shannon entropy (bits per base pair)
            entropy = -sum(f * log2(f) for f in freqs.values() if f > 0)

            # KL divergence from uniform
            uniform_p = 1.0 / n_possible
            kl_div = sum(
                f * log2(f / uniform_p)
                for f in freqs.values() if f > 0
            )

            # Top 20 over/underrepresented
            expected_freq = 1.0 / n_possible
            enrichment = {
                kmer: freqs.get(kmer, 0) / expected_freq
                for kmer in kmer_counts
            }
            sorted_enrichment = sorted(enrichment.items(), key=lambda x: x[1], reverse=True)
            top_over = sorted_enrichment[:20]

            # For underrepresented, include absent k-mers
            all_kmers_set = set("".join(c) for c in product(BASES, repeat=k))
            for kmer in all_kmers_set:
                if kmer not in enrichment:
                    enrichment[kmer] = 0.0
            sorted_under = sorted(enrichment.items(), key=lambda x: x[1])
            top_under = sorted_under[:20]

            entropy_summary.append({
                "k": k,
                "region": region_name,
                "entropy": entropy,
                "kl_divergence": kl_div,
                "total_kmers": total_kmers,
                "unique_kmers": len(kmer_counts),
            })

            print(f"    {region_name}: entropy={entropy:.4f} bits, "
                  f"KL_div={kl_div:.4f}, unique={len(kmer_counts)}/{n_possible}")

            # --- Write CSV for this region + k ---
            csv_path = os.path.join(kmer_outdir, f"kmer_k{k}_{region_name}.csv")
            rows = []
            for kmer in sorted(kmer_counts.keys()):
                rows.append({
                    "kmer": kmer,
                    "count": kmer_counts[kmer],
                    "frequency": round(freqs[kmer], 8),
                    "enrichment_vs_uniform": round(enrichment[kmer], 4),
                })
            pd.DataFrame(rows).to_csv(csv_path, index=False)

            # --- Write top over/underrepresented ---
            top_path = os.path.join(kmer_outdir, f"kmer_k{k}_{region_name}_extremes.csv")
            top_rows = []
            for kmer, enrich in top_over:
                top_rows.append({"kmer": kmer, "enrichment": round(enrich, 4), "direction": "over"})
            for kmer, enrich in top_under:
                top_rows.append({"kmer": kmer, "enrichment": round(enrich, 4), "direction": "under"})
            pd.DataFrame(top_rows).to_csv(top_path, index=False)

    # --- Entropy summary CSV ---
    summary_df = pd.DataFrame(entropy_summary)
    summary_df.to_csv(os.path.join(kmer_outdir, "entropy_summary.csv"), index=False)
    print(f"\n  Wrote entropy summary: {os.path.join(kmer_outdir, 'entropy_summary.csv')}")

    # --- Comparative heatmap ---
    _plot_kmer_heatmap(summary_df, kmer_outdir)

    print("  Task 3b complete.\n")
    return summary_df


def _plot_kmer_heatmap(summary_df, outdir):
    """Generate comparative heatmap of entropy across region types and k values."""
    # Pivot for heatmap: rows=k, cols=region
    pivot = summary_df.pivot(index="k", columns="region", values="entropy")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Entropy heatmap
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="viridis", ax=axes[0],
                cbar_kws={"label": "Shannon Entropy (bits)"})
    axes[0].set_title("Shannon Entropy by Region Type and k")
    axes[0].set_ylabel("k-mer length (k)")

    # KL divergence heatmap
    pivot_kl = summary_df.pivot(index="k", columns="region", values="kl_divergence")
    sns.heatmap(pivot_kl, annot=True, fmt=".3f", cmap="magma", ax=axes[1],
                cbar_kws={"label": "KL Divergence from Uniform"})
    axes[1].set_title("KL Divergence by Region Type and k")
    axes[1].set_ylabel("k-mer length (k)")

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "kmer_entropy_heatmap.png"), dpi=150)
    plt.close()
    print("  Wrote visualization: kmer_entropy_heatmap.png")


# ---------------------------------------------------------------------------
# 3c. Information Density Map
# ---------------------------------------------------------------------------

def analyze_density(genome, outdir, window_size=1000, step_size=100):
    """Task 3c: Information Density Map.

    Cryptanalytic Analogy:
        Entropy mapping is the information-theoretic scan for structure in
        ciphertext. High-entropy regions contain maximum information density —
        these are the complex regulatory programs, the "encrypted payloads"
        of the genome. Low-entropy regions are repetitive and predictable —
        structural elements, tandem repeats, the "known plaintext" that
        provides leverage for breaking the cipher. By mapping entropy across
        the genome, we identify the hotspots of regulatory complexity where
        the most sophisticated "code" resides.

    Computes:
        - Sliding window Shannon entropy (dinucleotide-based) across each chromosome
        - Flags regions >2 SD above/below genome-wide mean
        - BED file of flagged regions
        - Whole-genome entropy profile plots per chromosome
    """
    print("\n=== Task 3c: Information Density Map ===")
    print("Analogy: Entropy scan for information-dense zones in the ciphertext\n")

    density_outdir = os.path.join(outdir, "density")
    os.makedirs(density_outdir, exist_ok=True)

    all_entropies = []
    window_records = []  # (chr, start, end, entropy)

    for seqid, seq in sorted(genome.items()):
        seq_len = len(seq)
        if seq_len < window_size:
            print(f"  Skipping {seqid} (length {seq_len} < window size {window_size})")
            continue

        print(f"  Processing {seqid} ({seq_len:,} bp)...")
        chr_entropies = []

        for start in range(0, seq_len - window_size + 1, step_size):
            window = seq[start:start + window_size]
            entropy = _dinucleotide_entropy(window)
            chr_entropies.append(entropy)
            all_entropies.append(entropy)
            window_records.append((seqid, start, start + window_size, entropy))

    if not all_entropies:
        print("  No windows to analyze.")
        return

    all_entropies = np.array(all_entropies)
    mean_entropy = np.mean(all_entropies)
    std_entropy = np.std(all_entropies)
    threshold_high = mean_entropy + 2 * std_entropy
    threshold_low = mean_entropy - 2 * std_entropy

    print(f"\n  Genome-wide entropy: mean={mean_entropy:.4f}, std={std_entropy:.4f}")
    print(f"  High threshold (>2 SD): {threshold_high:.4f}")
    print(f"  Low threshold (<2 SD):  {threshold_low:.4f}")

    # Flag regions and write BED file
    bed_path = os.path.join(density_outdir, "flagged_regions.bed")
    high_count = 0
    low_count = 0
    with open(bed_path, "w") as fh:
        fh.write("#chrom\tstart\tend\tentropy\tflag\n")
        for chrom, start, end, entropy in window_records:
            if entropy > threshold_high:
                fh.write(f"{chrom}\t{start}\t{end}\t{entropy:.4f}\thigh_entropy\n")
                high_count += 1
            elif entropy < threshold_low:
                fh.write(f"{chrom}\t{start}\t{end}\t{entropy:.4f}\tlow_entropy\n")
                low_count += 1

    print(f"  Flagged regions: {high_count} high-entropy, {low_count} low-entropy")
    print(f"  Wrote BED file: {bed_path}")

    # Write full entropy profile CSV
    profile_path = os.path.join(density_outdir, "entropy_profile.csv")
    with open(profile_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["chrom", "start", "end", "entropy"])
        for rec in window_records:
            writer.writerow(rec)
    print(f"  Wrote entropy profile: {profile_path}")

    # Summary stats
    summary_path = os.path.join(density_outdir, "density_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("=== Information Density Map Summary ===\n\n")
        fh.write(f"Window size: {window_size} bp, Step size: {step_size} bp\n")
        fh.write(f"Total windows analyzed: {len(all_entropies):,}\n")
        fh.write(f"Genome-wide mean entropy: {mean_entropy:.4f}\n")
        fh.write(f"Genome-wide std entropy:  {std_entropy:.4f}\n")
        fh.write(f"High-entropy threshold:   {threshold_high:.4f}\n")
        fh.write(f"Low-entropy threshold:    {threshold_low:.4f}\n")
        fh.write(f"High-entropy windows:     {high_count}\n")
        fh.write(f"Low-entropy windows:      {low_count}\n")

    # Plot entropy profiles per chromosome
    _plot_entropy_profiles(window_records, mean_entropy, threshold_high, threshold_low,
                           density_outdir)

    print("  Task 3c complete.\n")


def _dinucleotide_entropy(sequence):
    """Compute Shannon entropy of a sequence using dinucleotide frequencies.

    Returns entropy in bits per base pair.
    """
    dinucs = Counter()
    clean = re.sub(r"[^ACGT]", "", sequence)
    n = len(clean)
    if n < 2:
        return 0.0

    for i in range(n - 1):
        dinucs[clean[i:i + 2]] += 1

    total = sum(dinucs.values())
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in dinucs.values():
        p = count / total
        if p > 0:
            entropy -= p * log2(p)

    return entropy


def _plot_entropy_profiles(window_records, mean_ent, high_thresh, low_thresh, outdir):
    """Plot whole-genome entropy profiles per chromosome."""
    # Group by chromosome
    chr_data = defaultdict(list)
    for chrom, start, end, entropy in window_records:
        chr_data[chrom].append((start, entropy))

    # Plot each chromosome
    for chrom, data in sorted(chr_data.items()):
        positions = [d[0] for d in data]
        entropies = [d[1] for d in data]

        fig, ax = plt.subplots(figsize=(16, 4))
        ax.plot(positions, entropies, linewidth=0.3, color="steelblue", alpha=0.7)
        ax.axhline(y=mean_ent, color="black", linestyle="--", linewidth=0.8, label="Mean")
        ax.axhline(y=high_thresh, color="red", linestyle=":", linewidth=0.8, label="+2 SD")
        ax.axhline(y=low_thresh, color="blue", linestyle=":", linewidth=0.8, label="-2 SD")
        ax.fill_between(positions, entropies, high_thresh,
                         where=[e > high_thresh for e in entropies],
                         color="red", alpha=0.3)
        ax.fill_between(positions, entropies, low_thresh,
                         where=[e < low_thresh for e in entropies],
                         color="blue", alpha=0.3)
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Shannon Entropy (bits)")
        ax.set_title(f"Information Density Profile — {chrom}")
        ax.legend(loc="upper right", fontsize=8)
        plt.tight_layout()

        safe_name = chrom.replace("|", "_").replace("/", "_")
        plt.savefig(os.path.join(outdir, f"entropy_profile_{safe_name}.png"), dpi=150)
        plt.close()

    print(f"  Wrote {len(chr_data)} chromosome entropy profile plots")


# ---------------------------------------------------------------------------
# 3d. Repeat Element Catalog
# ---------------------------------------------------------------------------

def analyze_repeats(genome, features, genome_lengths, outdir,
                    min_unit=1, max_unit=50, min_copies=3):
    """Task 3d: Repeat Element Catalog.

    Cryptanalytic Analogy:
        Repeat detection is the search for "key reuse" in the genomic
        ciphertext. In classical cryptanalysis, repeated sequences (like the
        Kasiski examination for Vigenère ciphers) reveal the key length and
        structure. In DNA, tandem repeats serve structural, regulatory, and
        even computational roles — microsatellites in promoters act as
        tunable expression dials, while larger repeats mark transposon
        insertion sites (horizontal "code imports"). Mapping these repeats
        catalogs the genome's reused keys and imported code blocks.

    Identifies:
        - Tandem repeats (exact and approximate)
        - Microsatellites (1-6 bp units), minisatellites (7-50 bp units)
        - Cross-references with gene annotations
    """
    print("\n=== Task 3d: Repeat Element Catalog ===")
    print("Analogy: Kasiski examination — finding key reuse in the ciphertext\n")

    repeat_outdir = os.path.join(outdir, "repeats")
    os.makedirs(repeat_outdir, exist_ok=True)

    # Build gene index for cross-referencing
    gene_index = _build_gene_index(features, genome_lengths)

    # Classify regions for labeling
    coding_ivs, regulatory_ivs, intergenic_ivs = classify_regions(features, genome_lengths)

    all_repeats = []

    for seqid, seq in sorted(genome.items()):
        print(f"  Scanning {seqid} ({len(seq):,} bp) for tandem repeats...")
        repeats = _find_tandem_repeats(seq, min_unit, max_unit, min_copies)
        print(f"    Found {len(repeats)} tandem repeat regions")

        for start, end, unit, copies in repeats:
            # Find nearest gene
            nearest_gene, distance = _find_nearest_gene(
                seqid, start, end, gene_index
            )
            # Classify region
            region_type = _classify_position(
                seqid, start, end, coding_ivs, regulatory_ivs, intergenic_ivs
            )
            # Categorize repeat
            unit_len = len(unit)
            if unit_len <= 6:
                repeat_class = "microsatellite"
            elif unit_len <= 50:
                repeat_class = "minisatellite"
            else:
                repeat_class = "satellite"

            all_repeats.append({
                "chr": seqid,
                "start": start,
                "end": end,
                "repeat_unit": unit,
                "unit_length": unit_len,
                "copy_count": round(copies, 1),
                "total_length": end - start,
                "repeat_class": repeat_class,
                "nearest_gene": nearest_gene or "none",
                "distance_to_gene": distance,
                "region_type": region_type,
            })

    df = pd.DataFrame(all_repeats)

    # Write catalog CSV
    csv_path = os.path.join(repeat_outdir, "repeat_catalog.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n  Total repeats cataloged: {len(df)}")
    print(f"  Wrote repeat catalog: {csv_path}")

    # Summary statistics
    if not df.empty:
        summary_path = os.path.join(repeat_outdir, "repeat_summary.txt")
        with open(summary_path, "w") as fh:
            fh.write("=== Repeat Element Catalog Summary ===\n\n")
            fh.write(f"Total repeats found: {len(df)}\n\n")

            fh.write("By class:\n")
            for cls, count in df["repeat_class"].value_counts().items():
                fh.write(f"  {cls}: {count}\n")

            fh.write("\nBy region type:\n")
            for region, count in df["region_type"].value_counts().items():
                fh.write(f"  {region}: {count}\n")

            fh.write(f"\nMean repeat length: {df['total_length'].mean():.1f} bp\n")
            fh.write(f"Mean copy count: {df['copy_count'].mean():.1f}\n")

            fh.write("\nTop 20 most common repeat units:\n")
            for unit, count in df["repeat_unit"].value_counts().head(20).items():
                fh.write(f"  {unit}: {count}\n")

        print(f"  Wrote summary: {summary_path}")

        # Visualization
        _plot_repeat_distribution(df, repeat_outdir)

    print("  Task 3d complete.\n")
    return df


def _find_tandem_repeats(sequence, min_unit=1, max_unit=50, min_copies=3):
    """Find tandem repeats in a sequence using a sliding window approach.

    Returns list of (start, end, repeat_unit, copy_count).
    """
    repeats = []
    seq = re.sub(r"[^ACGT]", "N", sequence.upper())
    seq_len = len(seq)

    for unit_len in range(min_unit, min(max_unit + 1, seq_len // min_copies + 1)):
        i = 0
        while i <= seq_len - unit_len * min_copies:
            unit = seq[i:i + unit_len]
            if "N" in unit:
                i += 1
                continue

            # Count consecutive copies (allowing ~10% mismatch for approximate)
            copies = 1
            j = i + unit_len
            while j + unit_len <= seq_len:
                candidate = seq[j:j + unit_len]
                mismatches = sum(1 for a, b in zip(unit, candidate) if a != b)
                if mismatches <= max(1, unit_len // 10):  # allow ~10% mismatch
                    copies += 1
                    j += unit_len
                else:
                    break

            if copies >= min_copies:
                end = i + copies * unit_len
                # Check this region wasn't already captured by a smaller unit
                is_subrepeat = False
                for existing_start, existing_end, existing_unit, _ in repeats:
                    if i >= existing_start and end <= existing_end:
                        is_subrepeat = True
                        break
                if not is_subrepeat:
                    repeats.append((i, end, unit, copies))
                i = end  # skip past this repeat
            else:
                i += 1

    return repeats


def _build_gene_index(features, genome_lengths):
    """Build an index of genes for nearest-gene lookups.

    Returns dict of seqid -> sorted list of (midpoint, gene_name, start, end).
    """
    gene_index = defaultdict(list)
    for f in features:
        if f["type"] == "gene":
            name = f["attributes"].get("Name", f["attributes"].get("ID", "unknown"))
            mid = (f["start"] + f["end"]) // 2
            gene_index[f["seqid"]].append((mid, name, f["start"] - 1, f["end"]))
    # Sort by midpoint for binary search
    for seqid in gene_index:
        gene_index[seqid].sort()
    return gene_index


def _find_nearest_gene(seqid, start, end, gene_index):
    """Find the nearest gene to a given interval.

    Returns (gene_name, distance). Distance is 0 if overlapping.
    """
    genes = gene_index.get(seqid, [])
    if not genes:
        return None, -1

    mid = (start + end) // 2

    # Binary search for closest gene midpoint
    import bisect
    idx = bisect.bisect_left([g[0] for g in genes], mid)

    best_gene = None
    best_dist = float("inf")

    for i in [idx - 1, idx, idx + 1]:
        if 0 <= i < len(genes):
            _, name, g_start, g_end = genes[i]
            # Distance between intervals
            if start >= g_end:
                dist = start - g_end
            elif end <= g_start:
                dist = g_start - end
            else:
                dist = 0  # overlapping
            if dist < best_dist:
                best_dist = dist
                best_gene = name

    return best_gene, best_dist


def _classify_position(seqid, start, end, coding_ivs, regulatory_ivs, intergenic_ivs):
    """Classify a genomic position as coding, regulatory, or intergenic."""
    mid = (start + end) // 2
    for ivs, label in [
        (coding_ivs, "coding"),
        (regulatory_ivs, "regulatory"),
        (intergenic_ivs, "intergenic"),
    ]:
        for iv_start, iv_end in ivs.get(seqid, []):
            if iv_start <= mid < iv_end:
                return label
    return "unclassified"


def _plot_repeat_distribution(df, outdir):
    """Plot repeat distribution visualizations."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1. Distribution by class
    class_counts = df["repeat_class"].value_counts()
    axes[0].bar(class_counts.index, class_counts.values, color=["#2196F3", "#FF9800", "#4CAF50"])
    axes[0].set_title("Repeats by Class")
    axes[0].set_ylabel("Count")

    # 2. Distribution by region type
    region_counts = df["region_type"].value_counts()
    axes[1].bar(region_counts.index, region_counts.values, color=["#9C27B0", "#F44336", "#00BCD4", "#607D8B"])
    axes[1].set_title("Repeats by Genomic Region")
    axes[1].set_ylabel("Count")
    axes[1].tick_params(axis="x", rotation=30)

    # 3. Unit length distribution
    axes[2].hist(df["unit_length"], bins=50, color="steelblue", edgecolor="white")
    axes[2].set_title("Repeat Unit Length Distribution")
    axes[2].set_xlabel("Unit Length (bp)")
    axes[2].set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "repeat_distribution.png"), dpi=150)
    plt.close()
    print("  Wrote visualization: repeat_distribution.png")


# ---------------------------------------------------------------------------
# CLI Interface
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="frequency_analyzer",
        description=(
            "GeneLang Statistical Fingerprinting Pipeline — "
            "Cryptanalytic frequency analysis of genomic data.\n\n"
            "Treats the genome as ciphertext and applies classical frequency "
            "analysis techniques to characterize codon usage, k-mer distributions, "
            "information density, and repeat elements across functional regions."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", help="Analysis to run")

    # Common arguments
    def add_common_args(sub):
        sub.add_argument("--genome", required=True, help="Path to genome FASTA file")
        sub.add_argument("--gff", required=True, help="Path to GFF3 annotation file")
        sub.add_argument("--outdir", required=True, help="Output directory")

    # 3a: Codon frequency
    codon_parser = subparsers.add_parser(
        "codon",
        help="3a: Codon frequency analysis (RSCU + chi-squared)",
        description=(
            "Count all 64 codons across CDS regions, compute RSCU and "
            "chi-squared bias tests. Analogous to letter frequency analysis "
            "in substitution cipher cracking."
        ),
    )
    codon_parser.add_argument("--cds", required=True, help="Path to CDS FASTA file")
    codon_parser.add_argument("--outdir", required=True, help="Output directory")

    # 3b: K-mer distribution
    kmer_parser = subparsers.add_parser(
        "kmer",
        help="3b: K-mer distribution analysis across genomic regions",
        description=(
            "Compute k-mer frequencies (k=2..8) separately for coding, "
            "regulatory, and intergenic regions. Computes Shannon entropy, "
            "KL divergence, and identifies over/underrepresented k-mers."
        ),
    )
    add_common_args(kmer_parser)

    # 3c: Information density
    density_parser = subparsers.add_parser(
        "density",
        help="3c: Information density map (sliding window entropy)",
        description=(
            "Sliding window Shannon entropy analysis across each chromosome. "
            "Flags high-entropy (regulatory complexity) and low-entropy "
            "(repetitive/structural) regions."
        ),
    )
    density_parser.add_argument("--genome", required=True, help="Path to genome FASTA file")
    density_parser.add_argument("--outdir", required=True, help="Output directory")
    density_parser.add_argument("--window", type=int, default=1000, help="Window size in bp (default: 1000)")
    density_parser.add_argument("--step", type=int, default=100, help="Step size in bp (default: 100)")

    # 3d: Repeat catalog
    repeat_parser = subparsers.add_parser(
        "repeats",
        help="3d: Repeat element catalog (tandem repeat detection)",
        description=(
            "Identify tandem repeats, catalog microsatellites and "
            "minisatellites, and cross-reference with gene annotations."
        ),
    )
    add_common_args(repeat_parser)
    repeat_parser.add_argument("--min-unit", type=int, default=1, help="Minimum repeat unit length (default: 1)")
    repeat_parser.add_argument("--max-unit", type=int, default=50, help="Maximum repeat unit length (default: 50)")
    repeat_parser.add_argument("--min-copies", type=int, default=3, help="Minimum copy count (default: 3)")

    # All analyses
    all_parser = subparsers.add_parser(
        "all",
        help="Run all analyses (3a + 3b + 3c + 3d)",
    )
    add_common_args(all_parser)
    all_parser.add_argument("--cds", required=True, help="Path to CDS FASTA file")
    all_parser.add_argument("--window", type=int, default=1000, help="Window size for density map (default: 1000)")
    all_parser.add_argument("--step", type=int, default=100, help="Step size for density map (default: 100)")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    print("=" * 60)
    print("GeneLang — Statistical Fingerprinting Pipeline")
    print("Cryptanalytic frequency analysis of genomic data")
    print("=" * 60)

    if args.command == "codon":
        analyze_codons(args.cds, args.outdir)

    elif args.command == "kmer":
        genome = parse_fasta(args.genome)
        features = parse_gff3(args.gff)
        genome_lengths = {sid: len(seq) for sid, seq in genome.items()}
        analyze_kmers(genome, features, genome_lengths, args.outdir)

    elif args.command == "density":
        genome = parse_fasta(args.genome)
        analyze_density(genome, args.outdir, args.window, args.step)

    elif args.command == "repeats":
        genome = parse_fasta(args.genome)
        features = parse_gff3(args.gff)
        genome_lengths = {sid: len(seq) for sid, seq in genome.items()}
        analyze_repeats(genome, features, genome_lengths, args.outdir)

    elif args.command == "all":
        # Load data once
        print("\nLoading genome...")
        genome = parse_fasta(args.genome)
        genome_lengths = {sid: len(seq) for sid, seq in genome.items()}
        total_bp = sum(genome_lengths.values())
        print(f"  Loaded {len(genome)} sequences, {total_bp:,} total bp")

        print("Loading annotations...")
        features = parse_gff3(args.gff)
        print(f"  Loaded {len(features)} features")

        # Run all analyses
        analyze_codons(args.cds, args.outdir)
        analyze_kmers(genome, features, genome_lengths, args.outdir)
        analyze_density(genome, args.outdir, args.window, args.step)
        analyze_repeats(genome, features, genome_lengths, args.outdir)

        print("=" * 60)
        print("All analyses complete.")
        print(f"Results in: {args.outdir}")
        print("=" * 60)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
