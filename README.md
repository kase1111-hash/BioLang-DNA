# GeneLang — Decoding the Programming Language of DNA

GeneLang applies cryptanalytic methodology to genomic data, treating DNA as an undeciphered programming language rather than a chemistry problem. DNA independently arrived at the same core computational operations found in every programming language — conditionals, loops, variables, error handling, modularity, and inheritance — because information processing demands them. Despite 70+ years of molecular biology, no one has produced a formal, readable grammar of the regulatory code that turns raw sequence into living behavior. This project aims to crack that code using the same techniques cryptographers use to break ciphers: frequency analysis, pattern matching, known-plaintext attacks, and differential analysis.

The project follows a four-tier organism ladder, starting with *C. elegans* (959 cells, fully mapped connectome, deterministic development) and progressing through *Drosophila melanogaster*, environmental-switching organisms (*Daphnia pulex*, *Pristionchus pacificus*), and ultimately *Apis mellifera* — where identical genomes produce queens or workers based on regulatory "compiler flags." Each tier builds the pattern libraries and analytical skills required for the next. The first target circuit is the *C. elegans* vulval induction pathway (RAS/MAPK signaling), the most completely mapped developmental decision in any animal.

The end goal is a human-readable formal specification of genetic regulatory programs — a "GeneLang" grammar — validated against known mutant phenotypes and generalizable across organisms. Success means the specification can predict what happens when you knock out a gene the same way a compiler spec predicts what happens when you delete a function. The full project spec is in [SPEC.md](SPEC.md), and the Phase 1 implementation plan is in [claude-code-kickoff.md](claude-code-kickoff.md).

## Repository Structure

```
genelang/
├── README.md
├── SPEC.md                          # Full project specification
├── grammar/                         # Extracted grammar rules and command mappings
├── organisms/
│   ├── c-elegans/                   # Tier 1: roundworm
│   │   ├── data/                    # Genome data (FASTA, GFF3, proteome)
│   │   ├── fingerprint/             # Statistical analysis outputs
│   │   ├── circuits/                # Formal circuit specifications
│   │   │   └── vulval-induction/    # Phase 1 target circuit
│   │   └── analysis/
│   ├── drosophila/                  # Tier 2: fruit fly
│   ├── daphnia/                     # Tier 3a: water flea
│   ├── pristionchus/                # Tier 3b: predatory roundworm
│   └── apis-mellifera/              # Tier 4: honeybee
├── tools/                           # Analysis pipelines
├── analysis/                        # Cross-organism comparative work
└── papers/                          # Publication drafts
```

## License

TBD
