#!/usr/bin/env python3
"""
tlCRISPRi Guide RNA Designer
=============================
Designs 28 bp CRISPRi spacer sequences (crRNA) targeting the RBS and CDS of user-defined genes.

Usage
-----
    python design_guides.py targets.csv genome.gb -o results.csv

Positional arguments
--------------------
    csv     Input CSV with at least two columns: a gene name column and a locus tag column.
    genome  Annotated GenBank file (preferred) or FASTA file of the genome.

Optional arguments
------------------
    -o / --output     Output CSV file path (default: stdout).
    --locus-col       Exact column header for the locus tag (auto-detected if omitted).
    --spacer-len      Length of the crRNA spacer in bp (default: 28).
    --rbs-window      Upstream window in bp to search for the RBS/SD sequence (default: 100).

Output columns added to the CSV
---------------------------------
    RBS         The identified RBS region (contains the Shine-Dalgarno sequence).
    CDS         The full CDS nucleotide sequence (annotated region, 5'→3').
    RBS spacer  28 bp spacer = reverse complement of the RBS region.
    CDS spacer  28 bp spacer = reverse complement of the first 28 bp of the CDS.

Dependencies
------------
    pip install biopython
"""

import sys
import csv
import argparse

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "Error: Biopython is required.\n"
        "Install it with:  pip install biopython"
    )

# Ordered from strongest to weakest Shine-Dalgarno consensus matches.
SD_PATTERNS = [
    "AGGAGG",
    "AAGGAG",
    "GAGGAG",
    "AGGAAG",
    "AGAAGG",
    "AGGAG",
    "AAGGAG",
    "AGGAG",
    "GAGG",
    "AAGG",
]


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence string."""
    return str(Seq(seq).reverse_complement())


def find_rbs_region(upstream_seq: str, spacer_len: int) -> str:
    """
    Locate the Shine-Dalgarno sequence within *upstream_seq* (5'→3', ending at the
    start codon) and return a *spacer_len*-bp window centred on it.

    If no SD-like pattern is found the last *spacer_len* bp are returned as a
    conservative fallback.
    """
    best_pos = -1
    best_pattern = ""

    for pattern in SD_PATTERNS:
        pos = upstream_seq.rfind(pattern)
        # Prefer the most downstream (closest to start codon) hit.
        if pos != -1 and pos > best_pos:
            best_pos = pos
            best_pattern = pattern

    if best_pos != -1:
        sd_center = best_pos + len(best_pattern) // 2
        start = sd_center - spacer_len // 2
        end = start + spacer_len

        # Clamp to sequence boundaries.
        if start < 0:
            start = 0
            end = min(spacer_len, len(upstream_seq))
        if end > len(upstream_seq):
            end = len(upstream_seq)
            start = max(0, end - spacer_len)

        return upstream_seq[start:end]

    # Fallback: no SD motif detected — use the region immediately upstream.
    if len(upstream_seq) >= spacer_len:
        return upstream_seq[-spacer_len:]
    return upstream_seq


def find_feature_by_locus_tag(genome_records, locus_tag: str):
    """
    Search all records for a CDS feature whose ``locus_tag`` or ``gene``
    qualifier matches *locus_tag*.  Returns ``(record, feature)`` or
    ``(None, None)`` if not found.
    """
    for record in genome_records:
        for feature in record.features:
            if feature.type != "CDS":
                continue
            if locus_tag in feature.qualifiers.get("locus_tag", []):
                return record, feature
            if locus_tag in feature.qualifiers.get("gene", []):
                return record, feature
    return None, None


def get_upstream_sequence(record, feature, window: int) -> str:
    """
    Return the *window*-bp sequence immediately upstream of the CDS start codon,
    oriented 5'→3' relative to the gene.

    For plus-strand genes  : record.seq[ start-window : start ]
    For minus-strand genes : reverse complement of record.seq[ end : end+window ]
    """
    strand = feature.location.strand
    if strand == 1:
        start = int(feature.location.start)
        upstream_start = max(0, start - window)
        return str(record.seq[upstream_start:start])
    else:
        end = int(feature.location.end)
        downstream_end = min(len(record.seq), end + window)
        return str(record.seq[end:downstream_end].reverse_complement())


def detect_format(path: str) -> str:
    """Guess the file format ('genbank' or 'fasta') from the file extension."""
    lower = path.lower()
    if any(lower.endswith(ext) for ext in (".gb", ".gbk", ".genbank", ".gbff")):
        return "genbank"
    if any(lower.endswith(ext) for ext in (".fa", ".fna", ".fasta", ".fas")):
        return "fasta"
    # Try GenBank first; fall back to FASTA.
    return "genbank"


def parse_genome(path: str):
    """Parse the genome file and return a list of SeqRecord objects."""
    formats_to_try = [detect_format(path)]
    if formats_to_try[0] == "genbank":
        formats_to_try.append("fasta")
    else:
        formats_to_try.append("genbank")

    for fmt in formats_to_try:
        try:
            records = list(SeqIO.parse(path, fmt))
            if records:
                return records
        except Exception:
            continue

    sys.exit(f"Error: Could not parse genome file '{path}'. "
             "Supported formats: GenBank (.gb, .gbk, .gbff) and FASTA (.fa, .fna, .fasta).")


def detect_locus_col(fieldnames: list[str]) -> str:
    """Heuristically identify the locus-tag column from the CSV header."""
    for col in fieldnames:
        if any(kw in col.lower() for kw in ("locus", "tag", "id", "accession")):
            return col
    # Fall back to the second column (first column is often the gene name).
    return fieldnames[1] if len(fieldnames) > 1 else fieldnames[0]


def process(csv_path: str, genome_path: str, output_path: str | None,
            locus_col_arg: str | None, spacer_len: int, rbs_window: int) -> None:

    genome_records = parse_genome(genome_path)

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            sys.exit("Error: Input CSV has no header row.")
        rows = list(reader)

    locus_col = locus_col_arg or detect_locus_col(fieldnames)
    if locus_col not in fieldnames:
        sys.exit(f"Error: Locus-tag column '{locus_col}' not found in CSV. "
                 f"Available columns: {fieldnames}")

    new_cols = ["RBS", "CDS", "RBS spacer", "CDS spacer"]
    out_fields = fieldnames + [c for c in new_cols if c not in fieldnames]

    out_rows = []
    not_found = []

    for row in rows:
        locus_tag = row.get(locus_col, "").strip()

        if not locus_tag:
            row.update(dict.fromkeys(new_cols, ""))
            out_rows.append(row)
            continue

        record, feature = find_feature_by_locus_tag(genome_records, locus_tag)

        if record is None:
            not_found.append(locus_tag)
            row.update(dict.fromkeys(new_cols, "NOT FOUND"))
            out_rows.append(row)
            continue

        # Full CDS sequence (5'→3', including start codon).
        cds_seq = str(feature.extract(record.seq))

        # Upstream sequence for RBS identification.
        upstream_seq = get_upstream_sequence(record, feature, rbs_window)

        # RBS region (centred on SD sequence).
        rbs_seq = find_rbs_region(upstream_seq, spacer_len)

        # crRNA spacers: reverse complement of the target DNA.
        rbs_spacer = reverse_complement(rbs_seq)[:spacer_len]
        cds_spacer = reverse_complement(cds_seq[:spacer_len])

        row.update({
            "RBS": rbs_seq,
            "CDS": cds_seq,
            "RBS spacer": rbs_spacer,
            "CDS spacer": cds_spacer,
        })
        out_rows.append(row)

    if not_found:
        print(
            f"Warning: {len(not_found)} locus tag(s) not found in genome: "
            + ", ".join(not_found),
            file=sys.stderr,
        )

    out_fh = open(output_path, "w", newline="") if output_path else sys.stdout
    try:
        writer = csv.DictWriter(out_fh, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)
    finally:
        if output_path:
            out_fh.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Design 28 bp CRISPRi spacer sequences targeting the RBS and CDS "
            "of genes listed in an input CSV."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("csv", help="Input CSV file with gene targets and locus tags.")
    parser.add_argument(
        "genome",
        help="Annotated GenBank (preferred) or FASTA genome file.",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="FILE",
        help="Write results to FILE instead of stdout.",
    )
    parser.add_argument(
        "--locus-col",
        metavar="COLUMN",
        help="Exact CSV column header containing locus tags (auto-detected if omitted).",
    )
    parser.add_argument(
        "--spacer-len",
        type=int,
        default=28,
        metavar="N",
        help="Length of the crRNA spacer in bp (default: 28).",
    )
    parser.add_argument(
        "--rbs-window",
        type=int,
        default=100,
        metavar="N",
        help="Upstream window in bp to search for the Shine-Dalgarno sequence (default: 100).",
    )

    args = parser.parse_args()
    process(
        csv_path=args.csv,
        genome_path=args.genome,
        output_path=args.output,
        locus_col_arg=args.locus_col,
        spacer_len=args.spacer_len,
        rbs_window=args.rbs_window,
    )


if __name__ == "__main__":
    main()
