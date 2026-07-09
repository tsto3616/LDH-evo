# this is for the formatting of the fa files for the csv output - run this script after the fa files have been generated to get the csv summary of the accessions in the fa files 

import csv
import re
from pathlib import Path

def parse_fasta_headers_to_table(fasta_path, isoform_name, out_csv):
    rows = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()

                # Example header:
                # HOG0002747[0]_African Elephant_XP_010600002.1
                parts = header.split("_")

                tree_label = parts[0]            # HOG0002747[0]
                accession = parts[-1]            # XP_010600002.1
                species = "_".join(parts[1:-1])  # handles multi-word species

                rows.append([isoform_name, tree_label, species, accession])

    with open(out_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["LDH_iso", "Tree", "Species", "Accession"])
        writer.writerows(rows)

# Example usage:
fasta_path = Path("phylo_hierarchy_0/HOG0002747.fa")
hog_id = fasta_path.stem  # extracts "HOG0002747"

outdir = Path("phylo_hierarchy_0")
outdir.mkdir(exist_ok=True)

parse_fasta_headers_to_table(
    fasta_path=fasta_path,
    isoform_name="LDHA",
    out_csv=outdir / f"{hog_id}_summary.csv"
)
