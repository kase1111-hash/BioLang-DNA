# GeneLang Phase 1 — C. elegans Cryptanalytic Fingerprinting

## Context

We're treating DNA as an undeciphered programming language and applying cryptanalytic methodology to decode the regulatory grammar. Full project spec is attached: `dna-decode-spec.md`. Read it first — it contains the compiler architecture model, organism ladder, and methodology.

This session starts Phase 1: statistical fingerprinting of the *C. elegans* genome, followed by formal specification of the vulval induction circuit.

## Task 1: Repository Setup

Create a GitHub repo called `genelang` with the directory structure from the spec. Initialize with the spec as `SPEC.md` and a README that summarizes the project in 2-3 paragraphs.

## Task 2: C. elegans Genome Acquisition

Download the *C. elegans* reference genome (WBcel235) from NCBI:
- Full genome FASTA: GCF_000002985.6
- Gene annotations: GFF3 format
- Protein sequences: predicted proteome FASTA
- Store in `organisms/c-elegans/data/`

Also pull the coding sequence (CDS) FASTA separately — we need coding regions isolated from non-coding for comparative frequency analysis.

## Task 3: Statistical Fingerprinting Pipeline

Build `tools/frequency_analyzer.py` that takes a genome FASTA + GFF3 annotation and produces:

### 3a. Codon Frequency Table
- Count all 64 codons across all CDS regions
- Compute relative synonymous codon usage (RSCU) — this reveals "dialect" preferences
- Compare observed vs expected frequencies (chi-squared test for each amino acid's codon family)
- Output: CSV + summary stats + visualization

### 3b. K-mer Distribution Analysis
- Compute k-mer frequencies for k=2 through k=8
- Do this SEPARATELY for three region types (use GFF3 to classify):
  - Coding regions (exons within CDS)
  - Regulatory regions (2kb upstream of each gene start = putative promoter)
  - Intergenic regions (everything else)
- For each region type and each k, compute:
  - Frequency distribution
  - Shannon entropy (bits per base pair)
  - Kullback-Leibler divergence from uniform distribution
  - Top 20 overrepresented and underrepresented k-mers
- Output: CSVs per region type per k, plus comparative heatmap showing entropy across region types and k values

### 3c. Information Density Map
- Sliding window entropy analysis across each chromosome
- Window size: 1000bp, step size: 100bp
- Shannon entropy per window using dinucleotide frequencies
- Flag regions with entropy significantly above or below genome-wide mean (>2 SD)
- High-entropy regions = information-dense (potential regulatory complexity)
- Low-entropy regions = repetitive/simple (tandem repeats, structural)
- Output: BED file of flagged regions + whole-genome entropy profile plot per chromosome

### 3d. Repeat Element Catalog
- Identify tandem repeats (exact and approximate) using a simple algorithm or TRF-style approach
- Catalog microsatellites, minisatellites, and larger repeat blocks
- Cross-reference repeat locations with gene annotations (are repeats near genes? in introns? in regulatory regions?)
- Output: repeat catalog CSV with columns: chr, start, end, repeat_unit, copy_count, nearest_gene, distance_to_gene, region_type

## Task 4: Regulatory Motif Scan (Initial)

Build `tools/motif_scanner.py`:
- Extract all sequences 2kb upstream of annotated gene starts (putative promoter regions)
- Search for TATA box variants (TATAAA and known variants within -25 to -35 of TSS)
- Search for CAAT box variants
- Search for GC box / Sp1 binding sites
- Count: what percentage of C. elegans genes have each of these "standard" promoter elements?
- For genes WITHOUT a TATA box, flag them as "alternative syntax" candidates
- Output: CSV mapping each gene to its detected promoter elements + summary statistics

## Task 5: Vulval Induction Gene Network Data Pull

This is prep for the formal specification work. Pull all available data on these genes from WormBase (use their REST API at https://wormbase.org/rest or scrape if needed):

| Gene | Role in Circuit |
|---|---|
| lin-3 | EGF ligand (function call from anchor cell) |
| let-23 | EGF receptor (event listener on VPCs) |
| sem-5 | Adaptor protein (signal routing) |
| let-60 | RAS GTPase (signal amplification) |
| lin-45 | RAF kinase (cascade step) |
| mek-2 | MEK kinase (cascade step) |
| mpk-1 | MAPK (cascade terminal output) |
| lin-12 | Notch receptor (lateral inhibition / mutex) |
| lin-1 | ETS transcription factor (conditional execution) |
| lin-31 | Winged-helix TF (conditional execution) |
| lin-39 | Hox gene (scope declaration for vulval competence) |
| lst-1 through lst-4 | Lateral signal targets |
| lip-1 | MAPK phosphatase (signal dampener / negative feedback) |

For each gene, collect:
- Sequence (genomic + CDS)
- Known regulatory elements (promoter, enhancers if annotated)
- Expression pattern (which cells, which developmental stages)
- Loss-of-function phenotype (what breaks when you knock it out)
- Gain-of-function phenotype (what happens with overexpression)
- Known protein-protein interactions
- Known regulatory targets (what does it activate/repress)

Store in `organisms/c-elegans/circuits/vulval-induction/gene-data/` as individual JSON files per gene plus a combined `network.json` that captures the interaction graph.

## Task 6: Differential Expression Scaffold

Build `tools/expression_differ.py` — a scaffold that can take two gene expression datasets (e.g., "VPC in 1° fate" vs "VPC in 3° fate") and compute:
- Differentially expressed genes (fold change + statistical test)
- Gene set enrichment against known pathways
- Regulatory module detection (co-expressed gene clusters)

This doesn't need real data yet — build it with a clean interface that accepts expression matrices (gene × sample) and outputs ranked gene lists with statistics. We'll feed it real GEO data in the next session.

## Constraints

- Use Python 3.10+ with standard scientific stack (numpy, scipy, pandas, matplotlib, seaborn, biopython)
- All outputs should be both CSV (for analysis) and visual (PNG plots)
- Every tool should be callable from CLI with `--help`
- Write docstrings explaining the cryptanalytic analogy for each analysis (e.g., "Codon frequency analysis is analogous to letter frequency analysis in substitution cipher cracking")
- Use the repo structure from the spec
- Commit incrementally with meaningful messages

## Priority Order

If you need to triage: Task 1 → Task 2 → Task 3a → Task 3b → Task 4 → Task 5 → Task 3c → Task 3d → Task 6

Get the codon frequencies and k-mer distributions running first — those are the "frequency analysis" foundation everything else builds on.
