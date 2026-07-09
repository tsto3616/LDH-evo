# everything is fully contained within this script, no other scripts are required
# fill in with the details you need (eg the path to the OrthoFinder output files, the column name for the platypus orthologs, and the list of RefSeq protein accessions for the platypus orthologs of the gene family of interest)

import re
import pandas as pd
from pathlib import Path
import subprocess
from Bio import Entrez

# this is the list of RefSeq protein accessions for the platypus orthologs of the gene family of interest
platypus_roots = [
    "NP_001121092.1",  
    "XP_028931082.1",  
    "XP_028931083.1",  
    "XP_028914549.1",  
    "XP_028914550.1"   
]

# proteome
platypus_col = "Platypus _GCF_004115215.2_mOrnAna1.pri.v4_translated_cds_ uncompressed"

# for retrieving sequences from NCBI Entrez, you must provide an email address
Entrez.email = "thomas.stocker@sydney.edu.au"

# EXECUTABLES
MAFFT_EXE = r"C:\Users\tsto3616\Downloads\mafft-win\mafft.bat"
IQTREE_EXE = r"C:\Users\tsto3616\Downloads\iqtree-3.1.3-Windows\bin\iqtree3.exe"

# now for parsing the data to be suitable for the NCBI retrieval
def extract_refseq_ids(cell):
    if pd.isna(cell):
        return []
    ids = []
    for entry in str(cell).split(","):
        m = re.search(r"(XP|NP|YP|WP)_[0-9]+\.[0-9]+", entry)
        if m:
            ids.append(m.group(0))
    return ids

# retrieve the rows of the dataframe that contain at least one of the platypus root accessions in the specified column
def filter_hogs_by_platypus(df, platypus_col, root_accessions):
    keep_rows = []
    for idx, row in df.iterrows():
        plat_ids = extract_refseq_ids(row[platypus_col])
        if any(pid in root_accessions for pid in plat_ids):
            keep_rows.append(idx)
    return df.loc[keep_rows].reset_index(drop=True)

# fetch the protein sequence in FASTA format from NCBI Entrez for a given accession
def fetch_protein_record(acc):
    h = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
    seq = h.read()
    h.close()
    return seq

# collect all ortholog sequences for a given HOG row, returning a list of tuples (species_col, protein_id, fasta_sequence)
def collect_ortholog_seqs_for_hog(row):
    species_cols = row.index[3:]  # after HOG, OG, Parent Clade
    seqs = []
    for col in species_cols:
        ids = extract_refseq_ids(row[col])
        for pid in ids:
            try:
                fasta = fetch_protein_record(pid)
            except Exception:
                continue
            seqs.append((col, pid, fasta))
    return seqs

# write the fasta sequences for a given HOG to a file, with headers formatted for IQ-TREE
def write_hog_fasta(hog_id, seqs, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = outdir / f"{hog_id}.fa"

    with open(fasta_path, "w") as f:
        for i, (species_col, pid, fasta) in enumerate(seqs):
            # IQ-TREE never strips bracketed tokens → guaranteed unique
            header = f">{hog_id}[{i}]_{species_col}_{pid}"
            seq_lines = [line.strip() for line in fasta.splitlines() if not line.startswith(">")]
            seq = "".join(seq_lines)
            f.write(header + "\n")
            f.write(seq + "\n")

    return fasta_path

# align through MAFFT
def run_mafft(fasta_path):
    aln_path = str(fasta_path) + ".aln.fa"
    cmd = [MAFFT_EXE, "--auto", str(fasta_path)]
    with open(aln_path, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    return aln_path

# count the number of sequences in a FASTA file (for analysis in R)
def count_sequences_in_fasta(fasta_path):
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

# perform the phylogeny
def run_iqtree(aln_path):
    cmd = [
        IQTREE_EXE,
        "-s", aln_path,
        "-m", "MFP",
        "-bb", "1000",
        "-nt", "AUTO"
    ]
    subprocess.run(cmd, check=True)


# This performs the entire process for a given hierarchy level (from OrthoFinder), filtering HOGs by platypus orthologs, fetching sequences, aligning them, and running IQ-TREE.
def process_hierarchy_for_phylogeny(path, platypus_col, root_accessions, outdir):
    df = pd.read_csv(path, sep="\t")
    filtered = filter_hogs_by_platypus(df, platypus_col, root_accessions)

    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    for _, row in filtered.iterrows():
        hog_id = row["HOG"]

        # Skip if already processed (treefile exists)
        treefile = outdir / f"{hog_id}.fa.aln.fa.treefile"
        if treefile.exists():
            print(f"Skipping {hog_id}: already processed.")
            continue

        seqs = collect_ortholog_seqs_for_hog(row)
        if len(seqs) < 3:
            print(f"Skipping {hog_id}: only {len(seqs)} sequences.")
            continue

        fasta_path = write_hog_fasta(hog_id, seqs, outdir)
        aln_path = run_mafft(fasta_path)

        nseq = count_sequences_in_fasta(aln_path)
        if nseq < 4:
            print(f"Skipping IQ-TREE for {hog_id}: only {nseq} sequences.")
            continue

        run_iqtree(aln_path)


# Perform for all levels and print to a directory named phylo_hierarchy_{level} for a tsv per level
for i in range(13):
    process_hierarchy_for_phylogeny(
        path=f"level-{i}_orthofinder.tsv",
        platypus_col=platypus_col,
        root_accessions=platypus_roots,
        outdir=f"phylo_hierarchy_{i}"
    )
