# test/

This folder contains test scaffolding for `design_guides.py`.

## Contents

| File | Description |
|------|-------------|
| `crRNA_test_generation.csv` | Three *R. palustris* CGA009 genes used for testing |
| `run_test.sh` | Shell script that invokes `design_guides.py` with the test CSV |
| `guides_output.csv` | Generated output (created after a successful run; not committed) |

---

## Requirements

### Python dependencies

```bash
pip install biopython
```

### Genome annotation file

`design_guides.py` requires an **annotated GenBank** file (`.gbk` or `.gbff`), **not** a plain FASTA (`.fna`).  
The script needs CDS feature annotations to locate the ribosome-binding site (RBS) and early coding-sequence (CDS) windows for each locus tag.

For *R. palustris* CGA009, download the annotated genome from NCBI:  
<https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000013785.1/>

The downloaded archive typically contains a file named `GCF_000013785.1_genomic.gbff`.

---

## Running the test

### Option A — direct Python call (from repo root)

```bash
python design_guides.py \
    test/crRNA_test_generation.csv \
    /path/to/GCF_000013785.1_genomic.gbff \
    -o test/guides_output.csv
```

Replace `/path/to/GCF_000013785.1_genomic.gbff` with the actual path to
your downloaded GenBank file.

### Option B — shell script (from repo root)

```bash
./test/run_test.sh /path/to/GCF_000013785.1_genomic.gbff
```

### Option C — shell script (from inside `test/`)

```bash
cd test
./run_test.sh /path/to/GCF_000013785.1_genomic.gbff
```

Running `run_test.sh` without arguments prints a usage message and exits:

```
ERROR: No GenBank file path provided.

Usage:
  ./test/run_test.sh /path/to/genome.gbk
...
```

---

## Expected output

A file `test/guides_output.csv` is created with the original columns from
`crRNA_test_generation.csv` plus four appended columns:

| Column | Description |
|--------|-------------|
| `RBS` | 50-nt window upstream of the CDS start codon (coding-strand orientation) |
| `CDS` | 50-nt window from the CDS start codon (coding-strand orientation) |
| `RBS spacer` | 20-nt protospacer targeting the RBS window (first NGG PAM site found) |
| `CDS spacer` | 20-nt protospacer targeting the early CDS window (first NGG PAM site found) |

Example header row:

```
Gene,Locus tag,RBS,CDS,RBS spacer,CDS spacer
```

---

## Test genes

| Gene | Locus tag | Organism |
|------|-----------|----------|
| phaR | RPA1795 | *Rhodopseudomonas palustris* CGA009 |
| ftsZ | RPA3522 | *Rhodopseudomonas palustris* CGA009 |
| phaZ | RPA1786 | *Rhodopseudomonas palustris* CGA009 |
