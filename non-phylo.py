# this is for the CDS retrieval from the LDH protein accessions - convert the tsv to csv and run this script to get the CDS sequences in FASTA format
# all is contained within this script
import pandas as pd

df = pd.read_csv("phylo_hierarchy_0/N0.HOG0009223_annotations.csv")

accessions = df["accession"].dropna().astype(str).unique().tolist()
print("Accessions collected:", accessions)

###############################################
# 2. Retrieve CDS from GenBank
###############################################

from Bio import Entrez, SeqIO
Entrez.email = "thomas.stocker@sydney.edu.au"

def fetch_cds_from_protein(protein_acc):
    try:
        handle = Entrez.efetch(db="protein", id=protein_acc,
                               rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
    except Exception as e:
        print(f"Error fetching {protein_acc}: {e}")
        return None, None

    for feature in record.features:
        if feature.type == "CDS":
            coded_by = feature.qualifiers.get("coded_by", [""])[0]
            nuc_acc = coded_by.split(":")[0]

            try:
                h2 = Entrez.efetch(db="nuccore", id=nuc_acc,
                                   rettype="fasta_cds_na", retmode="text")
                cds_text = h2.read()
                h2.close()
                return nuc_acc, cds_text
            except Exception as e:
                print(f"Error fetching CDS for {protein_acc}: {e}")
                return None, None

    print(f"No CDS feature found for {protein_acc}")
    return None, None

###############################################
# 3. Retrieve all CDS sequences
###############################################

cds_records = []

for acc in accessions:
    nuc_acc, cds = fetch_cds_from_protein(acc)
    if cds:
        cds_records.append((acc, nuc_acc, cds))
        print(f"Retrieved CDS for {acc} → {nuc_acc}")
    else:
        print(f"No CDS found for {acc}")

###############################################
# 4. Save CDS FASTA
###############################################

with open("LDHD_CDS.fa", "w") as f:
    for protein_acc, nuc_acc, cds in cds_records:
        f.write(cds)

