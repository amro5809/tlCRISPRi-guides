# tlCRISPRi-guides

A Python script for designing transcriptional-level CRISPRi (tlCRISPRi) guide
sequences (crRNA spacers) that target the ribosome-binding site (RBS) and
coding sequence (CDS) of bacterial genes.

---

## Requirements

- Python ≥ 3.8
- [BioPython](https://biopython.org/) ≥ 1.80
- [pandas](https://pandas.pydata.org/) ≥ 1.5

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

```
python design_guides.py <gene_list.csv> <genome.gbk> [options]
```

### Positional arguments

| Argument | Description |
|---|---|
| `gene_list.csv` | CSV file with gene targets and their locus tags |
| `genome.gbk` | Annotated GenBank file (`.gb`, `.gbk`, or `.gbff`) |

### Optional arguments

| Flag | Default | Description |
|---|---|---|
| `-o FILE` / `--output FILE` | stdout | Write output CSV to FILE |
| `--locus-tag-col COLUMN` | `locus_tag` | Column name for locus tags in the input CSV |
| `--upstream BP` | `100` | Upstream window (bp) to search for the RBS / Shine-Dalgarno sequence |
| `--spacer-length BP` | `28` | Length of each crRNA spacer (bp) |

### Example

```bash
python design_guides.py targets.csv genome.gbk -o guides_output.csv
```

---

## Input format

**`gene_list.csv`** – must contain at least a column with locus tags (default
column name: `locus_tag`).  Any additional columns are preserved in the output.

```
gene,locus_tag
lacZ,b0344
araC,b0064
```

**`genome.gbk`** – a fully annotated GenBank file for the organism.  Plain
FASTA files are not supported because they lack CDS feature annotations.

---

## Output columns

Four columns are appended to the input CSV:

| Column | Description |
|---|---|
| `RBS` | Suspected RBS region (5'→3', sense strand). A window of `--spacer-length` bp within `--upstream` bp of the start codon that contains the Shine-Dalgarno sequence (`AGGAGG` or a similar purine-rich sequence). |
| `CDS` | Full annotated coding sequence (5'→3', sense strand), including the ATG start codon. |
| `RBS spacer` | crRNA spacer complementary to the RBS region (reverse complement of the `RBS` column). |
| `CDS spacer` | crRNA spacer complementary to the first `--spacer-length` bp of the CDS, including the ATG start codon (reverse complement of the first 28 bp of the `CDS` column). |

---

## How the RBS region is identified

1. The script retrieves up to `--upstream` bp (default 100 bp) immediately 5'
   of the annotated start codon, presented in the gene's 5'→3' direction
   (strand-aware).
2. It searches for the Shine-Dalgarno consensus sequence `AGGAGG` and a series
   of common variants (`AAGGAG`, `GAGGAG`, `AGGAAG`, `AGGAG`, …).  The
   occurrence closest to the start codon is used.
3. A `--spacer-length`-bp (default 28 bp) window that covers the SD sequence
   and extends toward the start codon is selected as the RBS region.
4. If no SD-like sequence is found, the `--spacer-length` bp immediately
   upstream of the start codon are used as a fallback.

---

## crRNA spacer design

- **RBS spacer** = reverse complement of the RBS region.
  The spacer directs dCas9 to the RBS, sterically blocking ribosome assembly
  and thereby repressing translation.
- **CDS spacer** = reverse complement of the first 28 bp of the CDS (ATG…).
  The spacer directs dCas9 to the start of the coding sequence, blocking
  transcription elongation and/or translation initiation.

Both spacers are reported in the 5'→3' direction as they would be cloned
into the crRNA expression vector.

---

## Example calculation

Given an RBS region:

```
5'-AATTCATTAAAGAGGAGAAAGGTACC-3'
```

The RBS spacer (reverse complement) is:

```
5'-GGTACCTTTCTCCTCTTTAATGAATT-3'
```

Given a CDS starting with:

```
5'-ATGAAACGCATTAGCACCACCATTACC...-3'
```

The CDS spacer is the reverse complement of the first 28 bp of the CDS.
