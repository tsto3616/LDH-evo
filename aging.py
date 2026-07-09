# load the metadata file and compute the number of derived SNPs per lineage
from collections import defaultdict
def load_metadata(metadata_file):
    """
    metadata.txt format (TAB-separated):
    Individual    Latitude    Longitude    Elevation    Lin    site
    Ignore any non-tabbed lines (Abbreviations block).
    """
    mapping = {}

    with open(metadata_file) as f:
        header = f.readline()  # skip header

        for line in f:
            line = line.strip()
            if not line:
                continue

            # Skip non-tabbed lines (Abbreviations block)
            if "\t" not in line:
                continue

            parts = line.split("\t")

            # Real metadata rows have EXACTLY 6 columns
            if len(parts) != 6:
                continue

            sample = parts[0]      # KG505b
            lin = parts[4]         # SN, NRM, CR

            mapping[sample] = lin

    return mapping

def compute_linage(vcf_file, metadata_file):
    """
    Compute phylogenetic SNP age per lineage.
    LinAge(Lin) = number of SNPs where ANY sample in Lin has GT != 0/0
    """

    sample_to_lin = load_metadata(metadata_file)

    # Parse VCF header
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                raw_samples = header[9:]
                samples = [s.split("_")[0] for s in raw_samples]
                break

    lin_age = defaultdict(int)

    # Count derived SNPs per lineage
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            genotypes = fields[9:]

            for raw_sample, gt in zip(raw_samples, genotypes):
                sample = raw_sample.split("_")[0]
                lin = sample_to_lin.get(sample, "Unknown")
                gt = gt.split(":")[0]

                if gt in ("0/1", "1/1"):
                    lin_age[lin] += 1

    return lin_age

# write as a table
def write_linage_table(vcf_file, metadata_file, output_file):
    lin_age = compute_linage(vcf_file, metadata_file)

    with open(output_file, "w") as out:
        out.write("Lin\tLinAge\n")
        for lin, age in lin_age.items():
            out.write(f"{lin}\t{age}\n")

    print("LinAge per lineage written to", output_file)

# run the script
if __name__ == "__main__":
    write_linage_table(
        vcf_file="363_29709.vcf",
        metadata_file="indv_info_update.txt",
        output_file="LinAge_phylo.tsv"
    )
