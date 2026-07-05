#!/usr/bin/env python3
"""
design_guides.py — tlCRISPRi guide design from annotated GenBank + gene CSV.

Usage:
    python design_guides.py <csv_file> <genbank_file> [-o <output_csv>]

Arguments:
    csv_file      Path to CSV with columns: Gene, Locus tag
    genbank_file  Path to annotated GenBank file (.gbk / .gbff)
    -o            Output CSV path (default: guides_output.csv)

Output columns added to the input table:
    RBS           Sequence of the ribosome-binding site window (−50 to −1 nt upstream of CDS start)
    CDS           Sequence of the early coding region window (+1 to +50 nt from CDS start)
    RBS spacer    20-nt spacer sequence targeting the RBS region (NGG PAM, non-template strand)
    CDS spacer    20-nt spacer sequence targeting the early CDS region (NGG PAM, non-template strand)
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "ERROR: Biopython is required.\n"
        "Install it with:  pip install biopython"
    )


# ---------------------------------------------------------------------------
# Guide-design constants
# ---------------------------------------------------------------------------
RBS_UPSTREAM = 50       # nt upstream of CDS start to include in RBS window
CDS_DOWNSTREAM = 50     # nt downstream of CDS start to include in CDS window
SPACER_LENGTH = 20      # length of the protospacer (excludes PAM)
PAM = "GG"              # NGG PAM — the G-dinucleotide checked after the 'N'


def find_pam_spacer(sequence: str, label: str) -> str:
    """
    Return the first 20-nt spacer from *sequence* whose 3′ end is followed
    by an NGG PAM on the same strand (i.e., search for ...SPACER[N]GG...).

    If no NGG PAM is found, return an empty string and emit a warning.
    """
    seq = sequence.upper()
    for i in range(len(seq) - SPACER_LENGTH - 2):
        spacer_end = i + SPACER_LENGTH          # exclusive index
        pam_start = spacer_end + 1              # skip the 'N' in NGG
        if seq[pam_start : pam_start + 2] == PAM:
            return seq[i:spacer_end]
    print(f"  WARNING: no NGG PAM found in {label} window; leaving spacer blank.")
    return ""


def design_guides_for_feature(feature, record_seq: str, locus_tag: str) -> dict:
    """
    Given a CDS SeqFeature and the full chromosome sequence (as a string),
    return a dict with keys: RBS, CDS, RBS spacer, CDS spacer.
    """
    strand = feature.location.strand
    start = int(feature.location.start)   # 0-based, first nt of CDS
    end = int(feature.location.end)       # 0-based, exclusive

    if strand == 1:
        # Forward strand — RBS is upstream (lower indices), CDS is downstream
        rbs_start = max(0, start - RBS_UPSTREAM)
        rbs_seq = record_seq[rbs_start:start]
        cds_seq = record_seq[start : start + CDS_DOWNSTREAM]
    elif strand == -1:
        # Reverse strand — complement & reverse; upstream means higher indices
        rbs_end = min(len(record_seq), end + RBS_UPSTREAM)
        rbs_raw = record_seq[end:rbs_end]
        cds_raw = record_seq[end - CDS_DOWNSTREAM : end]
        # Reverse-complement so the guide is designed on the coding strand
        rbs_seq = str(Seq(rbs_raw).reverse_complement())
        cds_seq = str(Seq(cds_raw).reverse_complement())
    else:
        sys.exit(f"ERROR: CDS feature for {locus_tag} has unknown strand: {strand!r}")
    rbs_spacer = find_pam_spacer(rbs_seq, f"{locus_tag} RBS")
    cds_spacer = find_pam_spacer(cds_seq, f"{locus_tag} CDS")

    return {
        "RBS": rbs_seq,
        "CDS": cds_seq,
        "RBS spacer": rbs_spacer,
        "CDS spacer": cds_spacer,
    }


def load_genbank(gbk_path: str) -> dict:
    """
    Parse *gbk_path* and return a mapping of locus_tag → (feature, record_seq).
    """
    tag_map = {}
    path = Path(gbk_path)
    if not path.exists():
        sys.exit(f"ERROR: GenBank file not found: {gbk_path}")

    for record in SeqIO.parse(str(path), "genbank"):
        seq_str = str(record.seq).upper()
        for feature in record.features:
            if feature.type != "CDS":
                continue
            tags = feature.qualifiers.get("locus_tag", [])
            for tag in tags:
                tag_map[tag] = (feature, seq_str)
    return tag_map


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Design tlCRISPRi guides from a gene CSV and annotated GenBank file."
    )
    parser.add_argument("csv_file", help="Input CSV (columns: Gene, Locus tag)")
    parser.add_argument("genbank_file", help="Annotated GenBank file (.gbk / .gbff)")
    parser.add_argument(
        "-o", "--output", default="guides_output.csv",
        help="Output CSV path (default: guides_output.csv)"
    )
    args = parser.parse_args()

    # Validate inputs
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        sys.exit(f"ERROR: CSV file not found: {args.csv_file}")

    print(f"Loading GenBank annotations from: {args.genbank_file}")
    tag_map = load_genbank(args.genbank_file)
    print(f"  → {len(tag_map)} CDS features indexed.")

    # Read input CSV
    rows = []
    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        input_fieldnames = reader.fieldnames or []
        required = {"Gene", "Locus tag"}
        missing_cols = required - set(input_fieldnames)
        if missing_cols:
            sys.exit(f"ERROR: Input CSV is missing required columns: {', '.join(sorted(missing_cols))}")

        for row in reader:
            rows.append(row)

    if not rows:
        sys.exit("ERROR: Input CSV contains no data rows.")

    fieldnames = input_fieldnames + ["RBS", "CDS", "RBS spacer", "CDS spacer"]

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    found = 0
    missing = []

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            locus_tag = row.get("Locus tag", "").strip()
            gene = row.get("Gene", "").strip()
            print(f"  Processing {gene} ({locus_tag}) …")

            if locus_tag in tag_map:
                feature, seq_str = tag_map[locus_tag]
                guide_data = design_guides_for_feature(feature, seq_str, locus_tag)
                row.update(guide_data)
                found += 1
            else:
                print(f"  WARNING: locus tag '{locus_tag}' not found in GenBank; skipping.")
                row.update({"RBS": "", "CDS": "", "RBS spacer": "", "CDS spacer": ""})
                missing.append(locus_tag)

            writer.writerow(row)

    print(f"\nDone. {found}/{len(rows)} genes processed successfully.")
    if missing:
        print(f"Missing locus tags: {', '.join(missing)}")
    print(f"Output written to: {out_path}")


if __name__ == "__main__":
    main()
