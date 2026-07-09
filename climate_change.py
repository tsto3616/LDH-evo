#!/usr/bin/env python3
import csv

# search for the LDH genes in the GFF and extract their genomic coordinates
ldh_gene_names = ["Lactate Dehydrogenase A", "Lactate Dehydrogenase B", "Lactate Dehydrogenase D"]

# retrieve the LDH genes
def load_ldh_regions(gff_file, ldh_gene_names):
    """
    Parse GFF and return genomic intervals for LDH genes.
    ldh_gene_names = ["LDHA", "LDHB", "LDHC", "LDHD", "LDHAL6A", ...]
    """
    regions = []

    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attributes = fields

            # Only keep gene or mRNA features
            if feature not in ["gene", "mRNA", "exon", "CDS"]:
                continue

            # Check if LDH gene name appears in attributes
            for gene in ldh_gene_names:
                if gene.lower() in attributes.lower():
                    regions.append({
                        "gene": gene,
                        "chrom": chrom,
                        "start": int(start),
                        "end": int(end)
                    })

    return regions


# isolate SNPs that fall within the LDH gene intervals
def extract_ldh_snps(vcf_file, ldh_regions, output_file):
    """
    Scan VCF and extract SNPs that fall inside LDH gene intervals.
    """

    out = open(output_file, "w")
    writer = csv.writer(out, delimiter="\t")

    # Write header
    writer.writerow([
        "CHROM", "POS", "REF", "ALT", "GENE", "SAMPLES"
    ])

    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            samples = fields[9:]  # all genotype columns

            # Check each LDH region
            for region in ldh_regions:
                if chrom == region["chrom"] and region["start"] <= pos <= region["end"]:
                    writer.writerow([
                        chrom, pos, ref, alt, region["gene"], ",".join(samples)
                    ])
                    break

    out.close()
    
import csv

# reformat sample names to match metadata
def strip_name(name):
    return name.split("_")[0]

def load_metadata(metadata_file):
    meta = {}
    with open(metadata_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                raw = parts[0]
                clean = strip_name(raw)
                meta[clean] = " ".join(parts[1:])
    return meta

# extract altenrative alleles and metadata for a specific SNP
def extract_non_ref(tsv_file, vcf_file, metadata_file, chrom, pos, output_txt):

    # --- Load metadata ---
    metadata = load_metadata(metadata_file)

    # --- Load VCF sample names ---
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                raw_samples = header[9:]
                vcf_samples = [strip_name(s) for s in raw_samples]
                break

    # --- Load TSV row ---
    with open(tsv_file) as tsv, open(output_txt, "w") as out:
        reader = csv.DictReader(tsv, delimiter="\t")

        for row in reader:
            if row["CHROM"] == chrom and int(row["POS"]) == pos:

                genotypes = row["SAMPLES"].split(",")

                out.write(f"SNP: {chrom}:{pos}\n")
                out.write(f"REF: {row['REF']}\n")
                out.write(f"ALT: {row['ALT']}\n")
                out.write(f"GENE: {row['GENE']}\n\n")
                out.write("Non-reference samples:\n\n")

                found = False

                for name, gt in zip(vcf_samples, genotypes):
                    gt = gt.split(":")[0]

                    if gt != "0/0":
                        found = True
                        meta = metadata.get(name, "NoMetadata")
                        out.write(f"{name}\t{gt}\t{meta}\n")

                if not found:
                    out.write("No non-reference samples found.\n")

                break

    print(f"Done. Output written to {output_txt}")

# now run the main function to execute the script
if __name__ == "__main__":

    gff_file = "GCF_014633375.1_OchPri4.0_genomic.gff"
    vcf_file = "363_29709.vcf"
    output_file = "LDH_SNPs.tsv"

    print("Loading LDH genomic regions from GFF...")
    ldh_regions = load_ldh_regions(gff_file, ldh_gene_names)

    print(f"Found {len(ldh_regions)} LDH genomic intervals.")

    # After loading ldh_regions and before scanning the VCF
    ldh_chroms = sorted({r["chrom"] for r in ldh_regions})
    print("LDH chromosomes:", ldh_chroms)

    vcf_chroms = set()
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            vcf_chroms.add(chrom)
            if len(vcf_chroms) > 50:
                break

    print("VCF chromosomes (sample):", sorted(vcf_chroms))

    for r in ldh_regions[:10]:
        print(r)

    hits = 0

    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            samples = fields[9:]

            for region in ldh_regions:
                if chrom != region["chrom"]:
                    continue
                if region["start"] <= pos <= region["end"]:
                    hits += 1
                    print("HIT:", chrom, pos, ref, alt, region["gene"])
                    break

    print("Total LDH SNP hits:", hits)

    print("Extracting LDH SNPs from VCF...")
    extract_ldh_snps(vcf_file, ldh_regions, output_file)

    print(f"Done. LDH SNPs written to: {output_file}")

    extract_non_ref(
        tsv_file="LDH_SNPs.tsv",
        vcf_file="363_29709.vcf",
        metadata_file="indv_info_update.txt",
        chrom="NC_050564.1",
        pos=12726230,
        output_txt="LDH_SNP_12726230_with_metadata.txt"
    )

    