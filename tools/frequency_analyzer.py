#!/usr/bin/env python3
"""
frequency_analyzer.py — Statistical Fingerprinting Pipeline for Genomic Analysis

Applies cryptanalytic frequency analysis to genomic data, producing codon usage
tables, k-mer distributions, information density maps, and repeat element
catalogs. Each analysis treats the genome as "ciphertext" and extracts
statistical signatures that distinguish coding, regulatory, and inert regions.

Subcommands:
    codon    — Codon frequency table with RSCU and chi-squared analysis
    kmer     — K-mer distribution across coding/regulatory/intergenic regions
    entropy  — Sliding-window information density map
    repeats  — Tandem repeat element catalog

Usage:
    python tools/frequency_analyzer.py codon   --cds organisms/c-elegans/data/cds.fna -o organisms/c-elegans/fingerprint/
    python tools/frequency_analyzer.py kmer    --genome organisms/c-elegans/data/genome.fna --gff organisms/c-elegans/data/annotations.gff3 -o organisms/c-elegans/fingerprint/
    python tools/frequency_analyzer.py entropy --genome organisms/c-elegans/data/genome.fna -o organisms/c-elegans/fingerprint/
    python tools/frequency_analyzer.py repeats --genome organisms/c-elegans/data/genome.fna --gff organisms/c-elegans/data/annotations.gff3 -o organisms/c-elegans/fingerprint/
"""

import argparse
import os
import sys
from collections import Counter, defaultdict
from math import log, log2

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from Bio import SeqIO

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GENETIC_CODE = {
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
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODON_FAMILIES = defaultdict(list)
for _codon, _aa in sorted(GENETIC_CODE.items()):
    CODON_FAMILIES[_aa].append(_codon)

AMINO_ACID_NAMES = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    "*": "Stop",
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def log_msg(msg: str) -> None:
    """Print a timestamped log message to stderr."""
    print(f"[freq] {msg}", file=sys.stderr)


def load_genome(path: str) -> dict[str, str]:
    """Load a genome FASTA as {record_id: uppercase_sequence}."""
    log_msg(f"Loading genome from {path}")
    genome = {}
    for rec in SeqIO.parse(path, "fasta"):
        genome[rec.id] = str(rec.seq).upper()
    total_bp = sum(len(s) for s in genome.values())
    log_msg(f"  Loaded {len(genome)} sequences, {total_bp:,} bp total")
    return genome


def shannon_entropy(counts: Counter, base: float = 2.0) -> float:
    """Compute Shannon entropy in bits from a counter of observations."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    H = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            H -= p * log(p) / log(base)
    return H


def kl_divergence_from_uniform(counts: Counter, alphabet_size: int) -> float:
    """KL divergence D_KL(P || U) where U is uniform over alphabet_size."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    q = 1.0 / alphabet_size
    D = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            D += p * log2(p / q)
    return D


def parse_gff3_features(path: str, feature_types: set[str] | None = None):
    """
    Parse a GFF3 file, yielding feature dicts.

    Each dict has keys: seqid, type, start (1-based), end (1-based inclusive),
    strand, attributes (raw string).
    """
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            ftype = parts[2]
            if feature_types and ftype not in feature_types:
                continue
            yield {
                "seqid": parts[0],
                "type": ftype,
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "attributes": parts[8],
            }


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent intervals. Input: list of (start, end)."""
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _subtract_intervals(target: list[tuple[int, int]],
                         exclude: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Subtract exclude intervals from target intervals. All 0-based half-open."""
    result = []
    for t_start, t_end in target:
        remaining = [(t_start, t_end)]
        for e_start, e_end in exclude:
            new_remaining = []
            for r_start, r_end in remaining:
                if e_end <= r_start or e_start >= r_end:
                    new_remaining.append((r_start, r_end))
                else:
                    if r_start < e_start:
                        new_remaining.append((r_start, e_start))
                    if e_end < r_end:
                        new_remaining.append((e_end, r_end))
            remaining = new_remaining
        result.extend(remaining)
    return result


def classify_genomic_regions(gff_path: str, genome: dict[str, str]):
    """
    Classify genome into coding, regulatory, and intergenic intervals.

    Returns dict mapping region type to {seqid: [(start, end), ...]} with
    0-based half-open intervals.

    - Coding: CDS features from GFF3
    - Regulatory: 2kb upstream of each gene (minus any coding overlap)
    - Intergenic: everything else
    """
    log_msg("Classifying genomic regions from GFF3...")

    coding_raw = defaultdict(list)
    gene_upstreams_raw = defaultdict(list)

    for feat in parse_gff3_features(gff_path, {"CDS", "gene"}):
        sid = feat["seqid"]
        if sid not in genome:
            continue

        if feat["type"] == "CDS":
            coding_raw[sid].append((feat["start"] - 1, feat["end"]))

        elif feat["type"] == "gene":
            seq_len = len(genome[sid])
            if feat["strand"] == "+":
                up_end = feat["start"] - 1
                up_start = max(0, up_end - 2000)
            else:
                up_start = feat["end"]
                up_end = min(seq_len, up_start + 2000)
            if up_start < up_end:
                gene_upstreams_raw[sid].append((up_start, up_end))

    # Merge coding intervals
    coding = {}
    for sid in coding_raw:
        coding[sid] = _merge_intervals(coding_raw[sid])

    # Regulatory = upstream minus coding
    regulatory = {}
    for sid in gene_upstreams_raw:
        merged_up = _merge_intervals(gene_upstreams_raw[sid])
        regulatory[sid] = _subtract_intervals(merged_up, coding.get(sid, []))

    # Intergenic = genome minus coding minus regulatory
    intergenic = {}
    for sid, seq in genome.items():
        used = _merge_intervals(coding.get(sid, []) + regulatory.get(sid, []))
        gaps = []
        prev_end = 0
        for start, end in used:
            if prev_end < start:
                gaps.append((prev_end, start))
            prev_end = max(prev_end, end)
        if prev_end < len(seq):
            gaps.append((prev_end, len(seq)))
        intergenic[sid] = gaps

    def total_bp(ivs):
        return sum(e - s for regions in ivs.values() for s, e in regions)

    log_msg(f"  Coding:     {total_bp(coding):>12,} bp")
    log_msg(f"  Regulatory: {total_bp(regulatory):>12,} bp")
    log_msg(f"  Intergenic: {total_bp(intergenic):>12,} bp")

    return {"coding": coding, "regulatory": regulatory, "intergenic": intergenic}


def extract_sequences(genome: dict[str, str],
                      intervals: dict[str, list[tuple[int, int]]]) -> str:
    """Extract and concatenate sequences for a set of intervals."""
    parts = []
    for sid, ivs in intervals.items():
        seq = genome.get(sid, "")
        for start, end in ivs:
            parts.append(seq[start:end])
    return "".join(parts)


def count_kmers(sequence: str, k: int) -> Counter:
    """Count all k-mers in a sequence, skipping windows containing N."""
    counts = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if "N" not in kmer:
            counts[kmer] += 1
    return counts


def _parse_gene_name(attributes: str) -> str:
    """Extract gene name from GFF3 attributes string."""
    for attr in attributes.split(";"):
        if attr.startswith("Name="):
            return attr.split("=", 1)[1]
    for attr in attributes.split(";"):
        if attr.startswith("gene="):
            return attr.split("=", 1)[1]
    for attr in attributes.split(";"):
        if attr.startswith("ID="):
            return attr.split("=", 1)[1]
    return "unknown"


# ---------------------------------------------------------------------------
# Task 3a: Codon Frequency Table
# ---------------------------------------------------------------------------

def cmd_codon(args):
    """
    Codon frequency analysis with RSCU and chi-squared tests.

    Cryptanalytic analogy: Codon frequency analysis is analogous to letter
    frequency analysis in substitution cipher cracking. Just as English uses
    'E' far more than 'Z', organisms prefer certain synonymous codons over
    others. This codon usage bias reveals the organism's "dialect" and
    reflects selection pressures on translation efficiency and accuracy.
    """
    log_msg("=== Task 3a: Codon Frequency Analysis ===")
    os.makedirs(args.outdir, exist_ok=True)

    codon_counts = Counter()
    total_seqs = 0
    skipped = 0

    for rec in SeqIO.parse(args.cds, "fasta"):
        seq = str(rec.seq).upper().replace("U", "T")
        remainder = len(seq) % 3
        if remainder != 0:
            seq = seq[:len(seq) - remainder]
        if len(seq) < 3:
            skipped += 1
            continue
        total_seqs += 1
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if "N" not in codon and codon in GENETIC_CODE:
                codon_counts[codon] += 1

    total_codons = sum(codon_counts.values())
    log_msg(f"  Processed {total_seqs:,} CDS sequences ({skipped} skipped)")
    log_msg(f"  Total codons counted: {total_codons:,}")

    # Compute RSCU and chi-squared per amino acid family
    rows = []
    chi2_results = []

    for aa in sorted(CODON_FAMILIES.keys()):
        codons = CODON_FAMILIES[aa]
        n_syn = len(codons)
        family_total = sum(codon_counts.get(c, 0) for c in codons)

        if n_syn >= 2 and family_total > 0:
            observed = [codon_counts.get(c, 0) for c in codons]
            expected_val = family_total / n_syn
            chi2, pval = scipy_stats.chisquare(observed, f_exp=[expected_val] * n_syn)
            chi2_results.append({
                "amino_acid": aa,
                "name": AMINO_ACID_NAMES.get(aa, aa),
                "n_synonymous": n_syn,
                "total_count": family_total,
                "chi2": chi2,
                "p_value": pval,
                "significant": pval < 0.001,
            })

        for codon in codons:
            count = codon_counts.get(codon, 0)
            freq = count / total_codons if total_codons > 0 else 0
            if family_total > 0 and n_syn > 0:
                rscu = (count * n_syn) / family_total
            else:
                rscu = 0.0
            rows.append({
                "codon": codon,
                "amino_acid": aa,
                "aa_name": AMINO_ACID_NAMES.get(aa, aa),
                "count": count,
                "frequency": freq,
                "rscu": rscu,
                "n_synonymous": n_syn,
                "family_total": family_total,
            })

    # Write codon frequency CSV
    df = pd.DataFrame(rows)
    codon_csv = os.path.join(args.outdir, "codon_frequency.csv")
    df.to_csv(codon_csv, index=False, float_format="%.6f")
    log_msg(f"  Wrote {codon_csv}")

    # Write chi-squared CSV
    chi2_df = pd.DataFrame(chi2_results)
    chi2_csv = os.path.join(args.outdir, "codon_chi_squared.csv")
    chi2_df.to_csv(chi2_csv, index=False, float_format="%.4f")
    log_msg(f"  Wrote {chi2_csv}")

    # Write summary
    summary_path = os.path.join(args.outdir, "codon_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Codon Frequency Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"CDS sequences analyzed: {total_seqs:,}\n")
        f.write(f"Total codons counted:   {total_codons:,}\n")
        f.write(f"Sequences skipped:      {skipped}\n\n")

        f.write("Most preferred codons (highest RSCU):\n")
        top = df[df["n_synonymous"] > 1].nlargest(10, "rscu")
        for _, row in top.iterrows():
            f.write(f"  {row['codon']} ({row['aa_name']}): "
                    f"RSCU = {row['rscu']:.3f}, count = {row['count']:,}\n")

        f.write("\nMost avoided codons (lowest RSCU, multi-codon families):\n")
        bottom = df[df["n_synonymous"] > 1].nsmallest(10, "rscu")
        for _, row in bottom.iterrows():
            f.write(f"  {row['codon']} ({row['aa_name']}): "
                    f"RSCU = {row['rscu']:.3f}, count = {row['count']:,}\n")

        f.write("\nChi-squared tests (codon usage uniformity):\n")
        for _, row in chi2_df.iterrows():
            sig = "***" if row["significant"] else ""
            f.write(f"  {row['name']:4s} ({row['amino_acid']}): "
                    f"chi2 = {row['chi2']:10.1f}, "
                    f"p = {row['p_value']:.2e} {sig}\n")

    log_msg(f"  Wrote {summary_path}")

    # Visualization
    _plot_codon_rscu(df, args.outdir)
    log_msg("=== Codon analysis complete ===")


def _plot_codon_rscu(df: pd.DataFrame, outdir: str):
    """Generate RSCU bar chart and frequency bar chart."""
    plot_df = df[df["n_synonymous"] > 1].copy()
    plot_df = plot_df.sort_values(["amino_acid", "codon"])

    aa_list = plot_df["amino_acid"].unique()
    colors = sns.color_palette("husl", len(aa_list))
    aa_colors = dict(zip(aa_list, colors))

    # RSCU bar chart
    fig, ax = plt.subplots(figsize=(18, 6))
    x = range(len(plot_df))
    bar_colors = [aa_colors[aa] for aa in plot_df["amino_acid"]]
    ax.bar(x, plot_df["rscu"].values, color=bar_colors, edgecolor="none", width=0.8)
    ax.axhline(y=1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(plot_df["codon"].values, rotation=90, fontsize=7,
                       fontfamily="monospace")
    ax.set_ylabel("RSCU")
    ax.set_title("Relative Synonymous Codon Usage — C. elegans CDS")
    ax.set_xlim(-0.5, len(plot_df) - 0.5)

    # Amino acid group labels at top
    prev_aa = None
    group_start = 0
    for i, aa in enumerate(plot_df["amino_acid"].values):
        if aa != prev_aa:
            if prev_aa is not None:
                mid = (group_start + i - 1) / 2
                ax.text(mid, ax.get_ylim()[1] * 0.95, prev_aa, ha="center",
                        va="top", fontsize=8, fontweight="bold",
                        color=aa_colors[prev_aa])
            group_start = i
            prev_aa = aa
    if prev_aa is not None:
        mid = (group_start + len(plot_df) - 1) / 2
        ax.text(mid, ax.get_ylim()[1] * 0.95, prev_aa, ha="center",
                va="top", fontsize=8, fontweight="bold",
                color=aa_colors[prev_aa])

    plt.tight_layout()
    path = os.path.join(outdir, "codon_rscu.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    log_msg(f"  Wrote {path}")

    # Raw frequency bar chart (all 64 codons)
    all_df = df.sort_values(["amino_acid", "codon"])
    all_aa = all_df["amino_acid"].unique()
    all_colors_map = dict(zip(all_aa, sns.color_palette("husl", len(all_aa))))

    fig2, ax2 = plt.subplots(figsize=(18, 6))
    bar_colors2 = [all_colors_map[aa] for aa in all_df["amino_acid"]]
    ax2.bar(range(len(all_df)), all_df["frequency"].values, color=bar_colors2,
            edgecolor="none", width=0.8)
    ax2.set_xticks(range(len(all_df)))
    ax2.set_xticklabels(all_df["codon"].values, rotation=90, fontsize=7,
                        fontfamily="monospace")
    ax2.set_ylabel("Frequency (fraction of all codons)")
    ax2.set_title("Codon Frequency Distribution — C. elegans CDS")
    plt.tight_layout()
    path2 = os.path.join(outdir, "codon_frequency.png")
    fig2.savefig(path2, dpi=150)
    plt.close(fig2)
    log_msg(f"  Wrote {path2}")


# ---------------------------------------------------------------------------
# Task 3b: K-mer Distribution Analysis
# ---------------------------------------------------------------------------

def cmd_kmer(args):
    """
    K-mer distribution analysis across coding, regulatory, and intergenic regions.

    Cryptanalytic analogy: K-mer analysis is analogous to n-gram frequency
    analysis in cryptanalysis. Just as bigrams like 'TH' and 'HE' reveal
    English even under substitution, the frequency of short DNA motifs differs
    systematically between coding, regulatory, and inert regions. These
    differences are the statistical fingerprint of distinct "sublanguages"
    within the genome.
    """
    log_msg("=== Task 3b: K-mer Distribution Analysis ===")
    os.makedirs(args.outdir, exist_ok=True)

    genome = load_genome(args.genome)
    regions = classify_genomic_regions(args.gff, genome)

    region_seqs = {}
    for rtype in ("coding", "regulatory", "intergenic"):
        region_seqs[rtype] = extract_sequences(genome, regions[rtype])
        log_msg(f"  {rtype}: {len(region_seqs[rtype]):,} bp extracted")

    entropy_rows = []

    for k in range(args.kmin, args.kmax + 1):
        alphabet_size = 4 ** k
        log_msg(f"  Analyzing k={k} ({alphabet_size} possible {k}-mers)...")

        for rtype in ("coding", "regulatory", "intergenic"):
            seq = region_seqs[rtype]
            counts = count_kmers(seq, k)
            total = sum(counts.values())
            if total == 0:
                continue

            H = shannon_entropy(counts)
            max_H = log2(alphabet_size)
            D_KL = kl_divergence_from_uniform(counts, alphabet_size)

            entropy_rows.append({
                "k": k, "region": rtype, "total_kmers": total,
                "unique_kmers": len(counts), "possible_kmers": alphabet_size,
                "shannon_entropy": H, "max_entropy": max_H,
                "normalized_entropy": H / max_H if max_H > 0 else 0,
                "kl_divergence": D_KL,
            })

            # Full k-mer CSV
            kmer_df = pd.DataFrame([
                {"kmer": kmer, "count": cnt, "frequency": cnt / total}
                for kmer, cnt in counts.most_common()
            ])
            csv_path = os.path.join(args.outdir, f"kmer_k{k}_{rtype}.csv")
            kmer_df.to_csv(csv_path, index=False, float_format="%.8f")

            # Top/bottom 20 summary
            top20 = counts.most_common(20)
            all_sorted = counts.most_common()
            bottom20 = all_sorted[-20:] if len(all_sorted) > 20 else all_sorted
            bottom20.reverse()

            summary_path = os.path.join(args.outdir,
                                        f"kmer_k{k}_{rtype}_summary.txt")
            with open(summary_path, "w") as f:
                f.write(f"K-mer Analysis: k={k}, region={rtype}\n")
                f.write("=" * 50 + "\n")
                f.write(f"Total {k}-mers:      {total:,}\n")
                f.write(f"Unique {k}-mers:     {len(counts):,} / {alphabet_size}\n")
                f.write(f"Shannon entropy:    {H:.4f} bits (max {max_H:.4f})\n")
                if max_H > 0:
                    f.write(f"Normalized entropy: {H / max_H:.4f}\n")
                f.write(f"KL divergence:      {D_KL:.4f}\n\n")

                f.write("Top 20 overrepresented:\n")
                for kmer, cnt in top20:
                    f.write(f"  {kmer}  {cnt:>10,}  ({cnt / total:.6f})\n")
                f.write("\nBottom 20 underrepresented:\n")
                for kmer, cnt in bottom20:
                    f.write(f"  {kmer}  {cnt:>10,}  ({cnt / total:.6f})\n")

    # Entropy summary CSV
    entropy_df = pd.DataFrame(entropy_rows)
    entropy_csv = os.path.join(args.outdir, "kmer_entropy_summary.csv")
    entropy_df.to_csv(entropy_csv, index=False, float_format="%.6f")
    log_msg(f"  Wrote {entropy_csv}")

    _plot_kmer_heatmap(entropy_df, args.outdir)
    log_msg("=== K-mer analysis complete ===")


def _plot_kmer_heatmap(entropy_df: pd.DataFrame, outdir: str):
    """Comparative heatmap: entropy and KL divergence across regions and k."""
    pivot_H = entropy_df.pivot(index="k", columns="region",
                               values="normalized_entropy")
    pivot_KL = entropy_df.pivot(index="k", columns="region",
                                values="kl_divergence")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    sns.heatmap(pivot_H, annot=True, fmt=".3f", cmap="YlOrRd",
                ax=axes[0], vmin=0, vmax=1)
    axes[0].set_title("Normalized Shannon Entropy")
    axes[0].set_ylabel("k (mer length)")

    sns.heatmap(pivot_KL, annot=True, fmt=".3f", cmap="YlOrRd", ax=axes[1])
    axes[1].set_title("KL Divergence from Uniform")
    axes[1].set_ylabel("k (mer length)")

    plt.suptitle("K-mer Statistical Fingerprint — C. elegans Genome Regions",
                 fontsize=13)
    plt.tight_layout()
    path = os.path.join(outdir, "kmer_entropy_heatmap.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    log_msg(f"  Wrote {path}")


# ---------------------------------------------------------------------------
# Task 3c: Information Density Map
# ---------------------------------------------------------------------------

def cmd_entropy(args):
    """
    Sliding-window entropy analysis across each chromosome.

    Cryptanalytic analogy: Information density mapping is analogous to index
    of coincidence analysis — scanning ciphertext for regions of unusually
    high or low randomness. High-entropy regions suggest dense regulatory
    logic (complex, non-repetitive "instructions"), while low-entropy regions
    indicate structural or repetitive elements (the "padding" of the genome).
    """
    log_msg("=== Task 3c: Information Density Map ===")
    os.makedirs(args.outdir, exist_ok=True)

    genome = load_genome(args.genome)
    window = args.window
    step = args.step

    all_entropies = []
    chr_profiles = {}

    for sid, seq in sorted(genome.items()):
        if len(seq) < window:
            continue
        log_msg(f"  Processing {sid} ({len(seq):,} bp)...")
        positions = []
        entropies = []
        for i in range(0, len(seq) - window + 1, step):
            w = seq[i:i + window]
            di_counts = Counter()
            for j in range(len(w) - 1):
                dinuc = w[j:j + 2]
                if "N" not in dinuc:
                    di_counts[dinuc] += 1
            H = shannon_entropy(di_counts)
            positions.append(i)
            entropies.append(H)

        chr_profiles[sid] = (positions, entropies)
        all_entropies.extend(entropies)

    if not all_entropies:
        log_msg("  No sequences long enough for windowed analysis")
        return

    mean_H = float(np.mean(all_entropies))
    std_H = float(np.std(all_entropies))
    log_msg(f"  Genome-wide entropy: mean={mean_H:.4f}, std={std_H:.4f}")

    # Flag outlier regions and write BED
    bed_rows = []
    for sid, (positions, entropies) in chr_profiles.items():
        for pos, H in zip(positions, entropies):
            z = (H - mean_H) / std_H if std_H > 0 else 0
            if abs(z) > 2.0:
                label = "high_entropy" if z > 0 else "low_entropy"
                bed_rows.append(f"{sid}\t{pos}\t{pos + window}\t"
                                f"{label}\t{H:.4f}\t{z:.2f}")

    bed_path = os.path.join(args.outdir, "entropy_flagged_regions.bed")
    with open(bed_path, "w") as f:
        f.write("#chrom\tstart\tend\tname\tentropy\tz_score\n")
        for row in bed_rows:
            f.write(row + "\n")
    log_msg(f"  Wrote {bed_path} ({len(bed_rows):,} flagged windows)")

    n_high = sum(1 for r in bed_rows if "high_entropy" in r)
    n_low = sum(1 for r in bed_rows if "low_entropy" in r)

    summary_path = os.path.join(args.outdir, "entropy_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Information Density Map Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Window size:          {window} bp\n")
        f.write(f"Step size:            {step} bp\n")
        f.write(f"Total windows:        {len(all_entropies):,}\n")
        f.write(f"Mean entropy:         {mean_H:.4f} bits\n")
        f.write(f"Std entropy:          {std_H:.4f} bits\n")
        f.write(f"High-entropy windows: {n_high:,} (>2 SD above mean)\n")
        f.write(f"Low-entropy windows:  {n_low:,} (>2 SD below mean)\n")
    log_msg(f"  Wrote {summary_path}")

    _plot_entropy_profiles(chr_profiles, mean_H, std_H, args.outdir)
    log_msg("=== Entropy analysis complete ===")


def _plot_entropy_profiles(chr_profiles, mean_H, std_H, outdir):
    """Per-chromosome entropy profile plots."""
    n_chrs = len(chr_profiles)
    if n_chrs == 0:
        return

    fig, axes = plt.subplots(n_chrs, 1, figsize=(16, 2.5 * n_chrs),
                             sharex=False)
    if n_chrs == 1:
        axes = [axes]

    for ax, (sid, (positions, entropies)) in zip(
            axes, sorted(chr_profiles.items())):
        pos_mb = [p / 1e6 for p in positions]
        ax.plot(pos_mb, entropies, linewidth=0.3, color="steelblue", alpha=0.7)
        ax.axhline(y=mean_H, color="black", linestyle="-", linewidth=0.5,
                   alpha=0.5)
        ax.axhline(y=mean_H + 2 * std_H, color="red", linestyle="--",
                   linewidth=0.5, alpha=0.5)
        ax.axhline(y=mean_H - 2 * std_H, color="blue", linestyle="--",
                   linewidth=0.5, alpha=0.5)
        ax.set_ylabel("Entropy\n(bits)", fontsize=8)
        ax.set_title(sid, fontsize=9, loc="left")
        ax.tick_params(labelsize=7)

    axes[-1].set_xlabel("Position (Mb)")
    plt.suptitle("Information Density Profile — C. elegans Chromosomes",
                 fontsize=12, y=1.01)
    plt.tight_layout()
    path = os.path.join(outdir, "entropy_profile.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log_msg(f"  Wrote {path}")


# ---------------------------------------------------------------------------
# Task 3d: Repeat Element Catalog
# ---------------------------------------------------------------------------

def cmd_repeats(args):
    """
    Tandem repeat element identification and annotation.

    Cryptanalytic analogy: Repeat detection is analogous to Kasiski examination
    — identifying "key reuse" in polyalphabetic ciphers. Repeated sequences
    reveal copy-paste mechanisms, structural constraints, and evolutionary
    history. Repeats near genes may serve regulatory functions; repeats in
    intergenic regions may be structural scaffolding or evolutionary debris.
    """
    log_msg("=== Task 3d: Repeat Element Catalog ===")
    os.makedirs(args.outdir, exist_ok=True)

    genome = load_genome(args.genome)

    # Parse gene positions for cross-referencing
    log_msg("  Loading gene annotations...")
    gene_data = defaultdict(list)
    for feat in parse_gff3_features(args.gff, {"gene"}):
        sid = feat["seqid"]
        if sid not in genome:
            continue
        name = _parse_gene_name(feat["attributes"])
        gene_data[sid].append(
            (feat["start"] - 1, feat["end"], name, feat["strand"])
        )

    for sid in gene_data:
        gene_data[sid].sort()

    # Find tandem repeats
    log_msg("  Scanning for tandem repeats...")
    repeat_rows = []
    min_unit = args.min_unit
    max_unit = args.max_unit
    min_copies = args.min_copies

    for sid, seq in sorted(genome.items()):
        log_msg(f"    {sid} ({len(seq):,} bp)...")
        found = _find_tandem_repeats(seq, min_unit, max_unit, min_copies)
        genes = gene_data.get(sid, [])

        for start, end, unit, copies in found:
            unit_len = len(unit)
            if unit_len <= 6:
                rtype = "microsatellite"
            elif unit_len <= 100:
                rtype = "minisatellite"
            else:
                rtype = "satellite"

            nearest_gene, gene_dist, gene_dir = _find_nearest_gene(
                start, end, genes)

            repeat_rows.append({
                "chrom": sid, "start": start, "end": end,
                "repeat_unit": unit, "unit_length": unit_len,
                "copy_count": copies, "total_length": end - start,
                "repeat_type": rtype, "nearest_gene": nearest_gene,
                "distance_to_gene": gene_dist, "gene_direction": gene_dir,
            })

    if not repeat_rows:
        log_msg("  No repeats found")
        return

    repeat_df = pd.DataFrame(repeat_rows)
    csv_path = os.path.join(args.outdir, "repeat_catalog.csv")
    repeat_df.to_csv(csv_path, index=False)
    log_msg(f"  Wrote {csv_path} ({len(repeat_rows):,} repeats)")

    # Summary
    summary_path = os.path.join(args.outdir, "repeat_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Repeat Element Catalog Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total repeats found:  {len(repeat_rows):,}\n")
        f.write(f"Unit size range:      {min_unit}-{max_unit} bp\n")
        f.write(f"Minimum copies:       {min_copies}\n\n")
        for rtype in ("microsatellite", "minisatellite", "satellite"):
            subset = repeat_df[repeat_df["repeat_type"] == rtype]
            f.write(f"{rtype}s: {len(subset):,}\n")
            if len(subset) > 0:
                f.write(f"  Total span:    {subset['total_length'].sum():,} bp\n")
                f.write(f"  Mean copies:   {subset['copy_count'].mean():.1f}\n")
                top_units = subset["repeat_unit"].value_counts().head(5)
                f.write(f"  Top units:     {', '.join(top_units.index)}\n")
            f.write("\n")
        near_gene = repeat_df[repeat_df["distance_to_gene"] <= 1000]
        f.write(f"Repeats within 1kb of a gene: {len(near_gene):,} "
                f"({100 * len(near_gene) / len(repeat_df):.1f}%)\n")
    log_msg(f"  Wrote {summary_path}")

    _plot_repeat_summary(repeat_df, args.outdir)
    log_msg("=== Repeat analysis complete ===")


def _find_tandem_repeats(seq: str, min_unit: int, max_unit: int,
                         min_copies: int) -> list[tuple[int, int, str, int]]:
    """
    Find tandem repeats using numpy-accelerated period matching.

    For each candidate period (unit_len), compares each base to the base
    unit_len positions ahead. Consecutive runs of matches indicate tandem
    repeats. Only reports repeats whose unit is primitive (not itself a
    repetition of a smaller unit).

    Returns list of (start, end, unit_string, copy_count).
    """
    n = len(seq)
    if n == 0:
        return []

    seq_bytes = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    N_byte = ord("N")
    not_n = seq_bytes != N_byte
    repeats = []

    max_possible = min(max_unit, n // min_copies)
    for unit_len in range(min_unit, max_possible + 1):
        # Boolean: does position j match position j + unit_len?
        matches = (seq_bytes[:-unit_len] == seq_bytes[unit_len:])
        matches &= not_n[:-unit_len]
        matches &= not_n[unit_len:]

        # Find runs of consecutive True values
        padded = np.empty(len(matches) + 2, dtype=bool)
        padded[0] = False
        padded[-1] = False
        padded[1:-1] = matches

        diff = np.diff(padded.view(np.uint8).astype(np.int8))
        run_starts = np.where(diff == 1)[0]
        run_ends = np.where(diff == -1)[0]

        for rs, re in zip(run_starts, run_ends):
            run_len = re - rs
            copies = (run_len + unit_len) // unit_len
            if copies < min_copies:
                continue

            rep_start = rs
            rep_end = rs + copies * unit_len
            if rep_end > n:
                rep_end = n
                copies = (rep_end - rep_start) // unit_len
                if copies < min_copies:
                    continue

            unit = seq[rep_start:rep_start + unit_len]
            if _primitive_period(unit) == unit_len:
                repeats.append((rep_start, rep_end, unit, copies))

    return repeats


def _primitive_period(s: str) -> int:
    """Find the smallest period k such that s == (s[:k]) repeated."""
    n = len(s)
    for k in range(1, n):
        if n % k == 0 and s[:k] * (n // k) == s:
            return k
    return n


def _find_nearest_gene(start: int, end: int,
                       genes: list) -> tuple[str, int, str]:
    """Find the nearest gene to an interval. Returns (name, distance, direction)."""
    if not genes:
        return "none", -1, "none"

    mid = (start + end) // 2
    best_dist = float("inf")
    best_gene = "none"
    best_dir = "none"

    for g_start, g_end, g_name, g_strand in genes:
        if start < g_end and end > g_start:
            return g_name, 0, "overlapping"
        if mid < g_start:
            dist = g_start - mid
            direction = "upstream" if g_strand == "+" else "downstream"
        else:
            dist = mid - g_end
            direction = "downstream" if g_strand == "+" else "upstream"
        if dist < best_dist:
            best_dist = dist
            best_gene = g_name
            best_dir = direction

    return best_gene, int(best_dist), best_dir


def _plot_repeat_summary(repeat_df: pd.DataFrame, outdir: str):
    """Repeat unit length distribution and type breakdown."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].hist(repeat_df["unit_length"], bins=50, color="steelblue",
                 edgecolor="none", alpha=0.8)
    axes[0].set_xlabel("Repeat Unit Length (bp)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Repeat Unit Length Distribution")
    axes[0].set_yscale("log")

    type_counts = repeat_df["repeat_type"].value_counts()
    colors = {"microsatellite": "#2196F3", "minisatellite": "#FF9800",
              "satellite": "#4CAF50"}
    bar_colors = [colors.get(t, "gray") for t in type_counts.index]
    axes[1].bar(type_counts.index, type_counts.values, color=bar_colors)
    axes[1].set_ylabel("Count")
    axes[1].set_title("Repeat Types")

    plt.suptitle("Tandem Repeat Catalog — C. elegans", fontsize=13)
    plt.tight_layout()
    path = os.path.join(outdir, "repeat_summary.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    log_msg(f"  Wrote {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="frequency_analyzer",
        description=(
            "Statistical Fingerprinting Pipeline for Genomic Analysis.\n\n"
            "Applies cryptanalytic frequency analysis to genomic data,\n"
            "producing codon usage tables, k-mer distributions, information\n"
            "density maps, and repeat element catalogs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- codon ---
    p_codon = subparsers.add_parser(
        "codon",
        help="Codon frequency table with RSCU and chi-squared analysis",
        description=cmd_codon.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_codon.add_argument("--cds", required=True,
                         help="Path to CDS FASTA file")
    p_codon.add_argument("-o", "--outdir", default=".",
                         help="Output directory (default: current dir)")
    p_codon.set_defaults(func=cmd_codon)

    # --- kmer ---
    p_kmer = subparsers.add_parser(
        "kmer",
        help="K-mer distribution across genomic region types",
        description=cmd_kmer.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_kmer.add_argument("--genome", required=True,
                        help="Path to genome FASTA file")
    p_kmer.add_argument("--gff", required=True,
                        help="Path to GFF3 annotation file")
    p_kmer.add_argument("--kmin", type=int, default=2,
                        help="Minimum k-mer length (default: 2)")
    p_kmer.add_argument("--kmax", type=int, default=8,
                        help="Maximum k-mer length (default: 8)")
    p_kmer.add_argument("-o", "--outdir", default=".",
                        help="Output directory (default: current dir)")
    p_kmer.set_defaults(func=cmd_kmer)

    # --- entropy ---
    p_entropy = subparsers.add_parser(
        "entropy",
        help="Sliding-window information density map",
        description=cmd_entropy.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_entropy.add_argument("--genome", required=True,
                           help="Path to genome FASTA file")
    p_entropy.add_argument("--window", type=int, default=1000,
                           help="Window size in bp (default: 1000)")
    p_entropy.add_argument("--step", type=int, default=100,
                           help="Step size in bp (default: 100)")
    p_entropy.add_argument("-o", "--outdir", default=".",
                           help="Output directory (default: current dir)")
    p_entropy.set_defaults(func=cmd_entropy)

    # --- repeats ---
    p_repeats = subparsers.add_parser(
        "repeats",
        help="Tandem repeat element catalog",
        description=cmd_repeats.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_repeats.add_argument("--genome", required=True,
                           help="Path to genome FASTA file")
    p_repeats.add_argument("--gff", required=True,
                           help="Path to GFF3 annotation file")
    p_repeats.add_argument("--min-unit", type=int, default=1,
                           help="Minimum repeat unit size (default: 1)")
    p_repeats.add_argument("--max-unit", type=int, default=50,
                           help="Maximum repeat unit size (default: 50)")
    p_repeats.add_argument("--min-copies", type=int, default=3,
                           help="Minimum copy count (default: 3)")
    p_repeats.add_argument("-o", "--outdir", default=".",
                           help="Output directory (default: current dir)")
    p_repeats.set_defaults(func=cmd_repeats)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
