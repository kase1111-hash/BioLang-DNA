# C. elegans Genome Data Manifest

## Assembly

- **Species:** *Caenorhabditis elegans*
- **Assembly:** WBcel235
- **NCBI Accession:** GCF_000002985.6
- **Genome size:** ~100 Mb (6 chromosomes + mitochondrial)

## Files

| File | Description | Source | Format |
|------|-------------|--------|--------|
| `genome.fna` | Full reference genome | NCBI RefSeq | FASTA |
| `annotations.gff3` | Gene models, exons, CDS, UTRs | NCBI RefSeq | GFF3 |
| `proteome.faa` | Predicted protein sequences | NCBI RefSeq | FASTA |
| `cds.fna` | Coding sequences extracted from genome | NCBI RefSeq | FASTA |

## Acquisition

Run the download script from the repository root:

```bash
./organisms/c-elegans/download_genome.sh
```

Or download manually from:
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/
```

## Why These Files

- **genome.fna** — The complete "ciphertext." Contains all coding, regulatory, and intergenic sequences. Needed for whole-genome entropy analysis and repeat cataloging.
- **annotations.gff3** — The partial "codebook." Maps coordinates to known features (genes, exons, UTRs, regulatory elements). Required to classify regions for separate frequency analysis.
- **proteome.faa** — "Known plaintext." The protein outputs whose sequences are the solved portion of the code (codon→amino acid = simple substitution cipher, cracked in the 1960s).
- **cds.fna** — Coding sequences isolated from the genome. Needed for codon frequency analysis and RSCU computation without contamination from non-coding regions.
