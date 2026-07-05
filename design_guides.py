#!/usr/bin/env python3
"""
tlCRISPRi guide design script
==============================
Given a CSV file listing gene targets with their locus tags and an annotated
GenBank genome file, this script:

  1. Locates the annotated CDS for each locus tag.
  2. Identifies the suspected RBS region (a 28-bp window within ~100 bp upstream
     of the start codon that contains the Shine-Dalgarno sequence AGGAGG or a
     similar purine-rich sequence).
  3. Generates two 28-bp crRNA spacer sequences:
       - RBS spacer : reverse complement of the RBS region.
       - CDS spacer : reverse complement of the first 28 bp of the CDS
                      (starting at the ATG start codon).

Output columns appended to the input CSV:
  RBS        – suspected RBS region (5'→3', sense strand)
  CDS        – full annotated coding sequence (5'→3', sense strand)
  RBS spacer – 28-bp crRNA spacer targeting the RBS
  CDS spacer – 28-bp crRNA spacer targeting the start of the CDS

Usage
-----
    python design_guides.py gene_list.csv genome.gbk
    python design_guides.py gene_list.csv genome.gbk -o output.csv
    python design_guides.py gene_list.csv genome.gbk --locus-tag-col locus_tag
"""

import argparse
import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SPACER_LENGTH = 28       # length of each crRNA spacer (bp)
UPSTREAM_WINDOW = 100    # bp upstream of start codon to search for RBS

# Shine-Dalgarno consensus patterns ordered by decreasing similarity to the
# canonical sequence AGGAGG.  rfind() is used so that the instance closest to
# the start codon (i.e. the rightmost one) is preferred.
SD_PATTERNS = [
    "AGGAGG",  # canonical
    "AAGGAG",  # common variant
    "GAGGAG",  # common variant
    "AGGAAG",  # common variant
    "AGGAG",   # 5-mer core
    "AAGGA",
    "GGAGG",
    "AGGA",
    "GGAG",
    "AAGG",
    "AGGG",
]

# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence string."""
    return str(Seq(seq).reverse_complement())


def find_sd_position(upstream_seq: str):
    """
    Search *upstream_seq* for a Shine-Dalgarno sequence.

    Tries each pattern in SD_PATTERNS in order and returns the (start, end)
    coordinates of the rightmost (closest-to-start-codon) hit.  If no pattern
    matches, returns (None, None).
    """
    seq_upper = upstream_seq.upper()
    for pattern in SD_PATTERNS:
        pos = seq_upper.rfind(pattern)
        if pos != -1:
            return pos, pos + len(pattern)
    return None, None


# ---------------------------------------------------------------------------
# RBS extraction
# ---------------------------------------------------------------------------


def extract_rbs_region(record, feature, upstream_window: int, spacer_length: int) -> str:
    """
    Identify and return the *spacer_length*-bp RBS region for *feature*.

    Strategy
    --------
    1. Fetch up to *upstream_window* bp upstream of the start codon (accounting
       for strand), presented in the gene's 5'→3' direction.
    2. Search for a Shine-Dalgarno sequence.  Use the instance closest to the
       start codon.
    3. Build a *spacer_length*-bp window that covers the SD and extends toward
       the start codon.
    4. If no SD is found, use the *spacer_length* bp immediately upstream of the
       start codon.
    """
    strand = feature.location.strand

    if strand == 1:
        gene_start = int(feature.location.start)
        up_start = max(0, gene_start - upstream_window)
        upstream_seq = str(record.seq[up_start:gene_start])
    else:
        gene_end = int(feature.location.end)
        up_end = min(len(record.seq), gene_end + upstream_window)
        # Reverse-complement so the sequence is presented 5'→3' relative to the gene
        upstream_seq = reverse_complement(str(record.seq[gene_end:up_end]))

    sd_start, sd_end = find_sd_position(upstream_seq)

    if sd_start is not None:
        # Extend the window a few bases past the SD toward the start codon, then
        # take spacer_length bp ending at that point.
        window_end = min(len(upstream_seq), sd_end + 10)
        window_start = max(0, window_end - spacer_length)
        rbs_region = upstream_seq[window_start:window_end]
    else:
        # No SD found – use the spacer_length bp closest to the start codon
        rbs_region = (
            upstream_seq[-spacer_length:]
            if len(upstream_seq) >= spacer_length
            else upstream_seq
        )

    return rbs_region


# ---------------------------------------------------------------------------
# CDS extraction
# ---------------------------------------------------------------------------


def extract_cds_sequence(record, feature) -> str:
    """Return the full CDS sequence (5'→3') for *feature*.

    BioPython's :meth:`~Bio.SeqFeature.SeqFeature.extract` handles strand
    automatically, so the returned string always starts with the ATG start
    codon (when the annotation is correct).
    """
    return str(feature.extract(record.seq))


# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------


def build_locus_index(genome_records: list) -> dict:
    """
    Build a mapping of locus_tag → (record, feature) for every CDS feature
    in *genome_records*.
    """
    index = {}
    for record in genome_records:
        for feature in record.features:
            if feature.type == "CDS":
                for lt in feature.qualifiers.get("locus_tag", []):
                    index[lt] = (record, feature)
    return index


def process_gene_list(
    gene_df: pd.DataFrame,
    locus_col: str,
    genome_records: list,
    upstream_window: int,
    spacer_length: int,
) -> list:
    """
    For every row in *gene_df*, look up the CDS by locus tag and compute the
    RBS region, CDS sequence, and two crRNA spacers.

    Returns a list of dicts with keys: RBS, CDS, ``RBS spacer``, ``CDS spacer``.
    """
    locus_index = build_locus_index(genome_records)
    rows = []

    for _, row in gene_df.iterrows():
        locus_tag = str(row[locus_col]).strip()
        entry = {
            "RBS": None,
            "CDS": None,
            "RBS spacer": None,
            "CDS spacer": None,
        }

        if locus_tag in locus_index:
            record, feature = locus_index[locus_tag]

            rbs_seq = extract_rbs_region(record, feature, upstream_window, spacer_length)
            cds_seq = extract_cds_sequence(record, feature)

            entry["RBS"] = rbs_seq
            entry["CDS"] = cds_seq

            if rbs_seq:
                entry["RBS spacer"] = reverse_complement(rbs_seq)

            if cds_seq:
                first_n = cds_seq[:spacer_length]
                entry["CDS spacer"] = reverse_complement(first_n)
        else:
            print(
                f"Warning: locus tag '{locus_tag}' not found in genome.",
                file=sys.stderr,
            )

        rows.append(entry)

    return rows


# ---------------------------------------------------------------------------
# File-format detection
# ---------------------------------------------------------------------------


def detect_format(filepath: str):
    """
    Infer whether *filepath* is a GenBank or FASTA file.

    Checks the file extension first; falls back to inspecting the first line.
    Returns ``'genbank'``, ``'fasta'``, or ``None`` if undetermined.
    """
    lower = filepath.lower()
    if any(lower.endswith(ext) for ext in (".gb", ".gbk", ".genbank", ".gbff")):
        return "genbank"
    if any(lower.endswith(ext) for ext in (".fa", ".fasta", ".fna", ".ffn")):
        return "fasta"

    try:
        with open(filepath) as fh:
            first_line = fh.readline()
    except OSError:
        return None

    if first_line.startswith("LOCUS") or first_line.startswith("ID "):
        return "genbank"
    if first_line.startswith(">"):
        return "fasta"

    return None


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Design CRISPRi crRNA spacer sequences targeting the RBS and CDS of "
            "gene targets.  Reads a CSV file of gene targets with locus tags and "
            "an annotated GenBank genome file, then appends four columns: RBS, "
            "CDS, 'RBS spacer', and 'CDS spacer'."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "gene_list",
        help="CSV file with gene targets and locus tags.",
    )
    parser.add_argument(
        "genome",
        help="Annotated GenBank (.gb / .gbk / .gbff) genome file.",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="FILE",
        help="Write the output CSV to FILE instead of stdout.",
    )
    parser.add_argument(
        "--locus-tag-col",
        default="locus_tag",
        metavar="COLUMN",
        help="Name of the locus-tag column in the input CSV.",
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=UPSTREAM_WINDOW,
        metavar="BP",
        help="Size of the upstream window (bp) searched for the RBS / SD sequence.",
    )
    parser.add_argument(
        "--spacer-length",
        type=int,
        default=SPACER_LENGTH,
        metavar="BP",
        help="Length of each crRNA spacer sequence (bp).",
    )

    args = parser.parse_args()

    # ---- Load gene list ---------------------------------------------------
    try:
        gene_df = pd.read_csv(args.gene_list)
    except Exception as exc:
        sys.exit(f"Error reading gene list CSV '{args.gene_list}': {exc}")

    if args.locus_tag_col not in gene_df.columns:
        sys.exit(
            f"Error: column '{args.locus_tag_col}' not found in "
            f"'{args.gene_list}'.  Available columns: {list(gene_df.columns)}"
        )

    # ---- Load genome -------------------------------------------------------
    fmt = detect_format(args.genome)
    if fmt is None:
        sys.exit(
            f"Error: could not determine the format of '{args.genome}'.  "
            "Please use a file with a .gb / .gbk / .gbff extension."
        )
    if fmt == "fasta":
        sys.exit(
            f"Error: '{args.genome}' appears to be a plain FASTA file, which "
            "does not carry CDS annotations.  Please supply an annotated "
            "GenBank file (.gb / .gbk / .gbff)."
        )

    try:
        genome_records = list(SeqIO.parse(args.genome, fmt))
    except Exception as exc:
        sys.exit(f"Error reading genome file '{args.genome}': {exc}")

    if not genome_records:
        sys.exit(f"Error: no sequence records found in '{args.genome}'.")

    # ---- Process -----------------------------------------------------------
    results = process_gene_list(
        gene_df,
        locus_col=args.locus_tag_col,
        genome_records=genome_records,
        upstream_window=args.upstream,
        spacer_length=args.spacer_length,
    )

    # ---- Build output dataframe --------------------------------------------
    out_df = gene_df.copy()
    out_df["RBS"] = [r["RBS"] for r in results]
    out_df["CDS"] = [r["CDS"] for r in results]
    out_df["RBS spacer"] = [r["RBS spacer"] for r in results]
    out_df["CDS spacer"] = [r["CDS spacer"] for r in results]

    csv_text = out_df.to_csv(index=False)

    if args.output:
        try:
            with open(args.output, "w") as fh:
                fh.write(csv_text)
        except OSError as exc:
            sys.exit(f"Error writing output file '{args.output}': {exc}")
        print(f"Results written to '{args.output}'.")
    else:
        sys.stdout.write(csv_text)


if __name__ == "__main__":
    main()
