# this one requries the MACSE alignment to be done externally, and then the pairwise dN/dS is computed from that alignment
FASTA_FILE      = "LDHD_CDS_NT.fa"              # input CDS FASTA
MACSE_ALN_FILE  = "LDHD_CDS_NT.fa"    # MACSE nucleotide alignment output

import numpy as np
import plotly.express as px
import csv
from Bio import SeqIO, AlignIO
from Bio.Data import CodonTable
from Bio import Entrez

# load the fasta file and return a list of tuples (accession, sequence)
def load_cds_fasta(fasta_file):
    cds_records = []
    for rec in SeqIO.parse(fasta_file, "fasta"):
        cds_records.append((rec.id, str(rec.seq).upper()))
    return cds_records

def load_macse_alignment(aln_file):
    return AlignIO.read(aln_file, "fasta")

# query genbank for species name given an accession
Entrez.email = "thomas.stocker@sydney.edu.au"

def extract_real_accession(header):
    """
    Converts:
    lcl|XM_037837472.1_cds_XP_037693400.1_1
    → XM_037837472.1
    """
    if header.startswith("lcl|"):
        header = header.replace("lcl|", "")
    return header.split("_cds_")[0]

def get_species_from_genbank(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = handle.read()
        handle.close()

        for line in record.split("\n"):
            if line.strip().startswith("ORGANISM"):
                return line.replace("ORGANISM", "").strip()

        return "Unknown_species"

    except Exception as e:
        print(f"Failed to fetch {accession}: {e}")
        return "Unknown_species"

def build_species_lookup(cds_records):
    species_map = {}
    for acc, seq in cds_records:
        real_acc = extract_real_accession(acc)
        species_map[acc] = get_species_from_genbank(real_acc)
    return species_map

standard_table = CodonTable.unambiguous_dna_by_id[1]

# calculate dN/dS manually from two aligned sequences
def manual_dn_ds(aln_seq1, aln_seq2):
    codons1 = [aln_seq1[i:i+3] for i in range(0, len(aln_seq1), 3)]
    codons2 = [aln_seq2[i:i+3] for i in range(0, len(aln_seq2), 3)]

    nonsyn = syn = nonsyn_sites = syn_sites = 0

    for c1, c2 in zip(codons1, codons2):
        if "-" in c1 or "-" in c2:
            continue
        if len(c1) != 3 or len(c2) != 3:
            continue
        if c1 not in standard_table.forward_table or c2 not in standard_table.forward_table:
            continue

        aa1 = standard_table.forward_table[c1]
        aa2 = standard_table.forward_table[c2]

        syn_sites += 1 if aa1 == aa2 else 0
        nonsyn_sites += 1 if aa1 != aa2 else 0

        if c1 != c2:
            if aa1 == aa2:
                syn += 1
            else:
                nonsyn += 1

    dN = nonsyn / nonsyn_sites if nonsyn_sites else np.nan
    dS = syn / syn_sites if syn_sites else np.nan
    return dN, dS

# compute pairwise dN/dS for all sequences in the MSA
def compute_dn_ds_matrix(msa):
    dn_ds_results = {}
    n = len(msa)

    for i in range(n):
        for j in range(i+1, n):
            seq1 = str(msa[i].seq)
            seq2 = str(msa[j].seq)
            dN, dS = manual_dn_ds(seq1, seq2)
            omega = dN / dS if (dS and not np.isnan(dS)) else np.nan
            dn_ds_results[(msa[i].id, msa[j].id)] = omega

    return dn_ds_results

# compute mean dN/dS per CDS and include species name
def mean_omega_per_cds_with_species(cds_records, dn_ds_results, species_map):
    """
    Compute mean dN/dS per CDS and include species name.
    """
    omega_means = {}

    for acc, seq in cds_records:
        vals = []

        # collect all ω values involving this CDS
        for (a, b), w in dn_ds_results.items():
            if a == acc or b == acc:
                if w is not None and not np.isnan(w):
                    vals.append(w)

        mean_w = np.mean(vals) if vals else np.nan
        species = species_map.get(acc, "Unknown_species")

        omega_means[acc] = (species, mean_w)

    return omega_means

# run it - run it for each gene family too
if __name__ == "__main__":
    cds_records = load_cds_fasta(FASTA_FILE)

    # Query GenBank for species names
    species_map = build_species_lookup(cds_records)

    # Load MACSE alignment
    msa = load_macse_alignment(MACSE_ALN_FILE)

    # Compute pairwise dN/dS
    dn_ds_results = compute_dn_ds_matrix(msa)

    # Compute mean ω per CDS (with species)
    omega_means = mean_omega_per_cds_with_species(cds_records, dn_ds_results, species_map)

    # Save CSV
    with open("LDHD_mean_dnds_per_CDS.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["CDS_ID", "Species", "Mean_dN_dS"])
        for acc, (species, w) in omega_means.items():
            writer.writerow([acc, species, w])
