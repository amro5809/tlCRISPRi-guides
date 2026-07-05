# tlCRISPRi-guides

A tool for designing transcription-level CRISPRi (tlCRISPRi) guide RNA spacers
from an annotated GenBank genome and a CSV list of target genes.

---

## Usage

```bash
pip install biopython

python design_guides.py <gene_csv> <genome.gbk> -o <output.csv>
```

- `<gene_csv>` — CSV with columns `Gene` and `Locus tag`
- `<genome.gbk>` — Annotated GenBank file (`.gbk` / `.gbff`); **plain FASTA is not supported**
- `-o` — Output CSV path (default: `guides_output.csv`)

The output CSV extends the input table with four columns:
`RBS`, `CDS`, `RBS spacer`, and `CDS spacer`.

---

## Testing

A ready-to-run test scaffold is provided in the [`test/`](test/) folder, including
a 3-gene CSV for *Rhodopseudomonas palustris* CGA009 and a helper shell script.

```bash
# Run from repo root (supply path to your annotated GenBank file)
./test/run_test.sh /path/to/GCF_000013785.1_genomic.gbff
```

See [`test/README.md`](test/README.md) for full instructions.