ldh_gene_names = {"LDHD": "comp131078_c0_seq1", # this is from the transcriptome assembly, not the genome - manually searched eggNOG for LDH genes and found these transcripts
             "LDHC": "comp57240_c0_seq1",
             "LDHB": "comp88044_c0_seq1"}

import csv

def extract_ldh_snps_from_vcf(vcf_file, LDH_genes, output_file):
    """
    Extract SNPs from a transcript-mapped VCF.
    CHROM = transcript ID (e.g., comp131078_c0_seq1)
    POS   = CDS position
    """

    out = open(output_file, "w")
    writer = csv.writer(out, delimiter="\t")

    writer.writerow(["CHROM", "POS", "REF", "ALT", "GENE", "SAMPLES"])

    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]          # transcript ID
            pos = int(fields[1])       # CDS position
            ref = fields[3]
            alt = fields[4]
            samples = fields[9:]

            # Check if CHROM matches any LDH transcript
            for gene_name, transcript_id in LDH_genes.items():
                if chrom == transcript_id:
                    writer.writerow([
                        chrom, pos, ref, alt, gene_name, ",".join(samples)
                    ])
                    break

    out.close()
    print(f"LDH SNPs written to {output_file}")

extract_ldh_snps_from_vcf("end-pleistocene.vcf", ldh_gene_names, "LDH_END_PILO.tsv")

def extract_all_aberrant(tsv_file, vcf_file, output_txt):

    # --- Load VCF sample names (raw) ---
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                vcf_samples = header[9:]   # use raw names
                break

    # --- Process all TSV rows ---
    with open(tsv_file) as tsv, open(output_txt, "w") as out:
        reader = csv.DictReader(tsv, delimiter="\t")

        for row in reader:
            chrom = row["CHROM"]
            pos   = int(row["POS"])
            ref   = row["REF"]
            alt   = row["ALT"]
            gene  = row["GENE"]
            genotypes = row["SAMPLES"].split(",")

            out.write(f"SNP: {chrom}:{pos}\n")
            out.write(f"REF: {ref}\n")
            out.write(f"ALT: {alt}\n")
            out.write(f"GENE: {gene}\n\n")
            out.write("Aberrant genotypes (non-reference):\n")

            found = False

            for sample_name, gt in zip(vcf_samples, genotypes):
                gt_clean = gt.split(":")[0]

                # Only print aberrant genotypes
                if gt_clean in ("0/1", "1/1"):
                    found = True
                    out.write(f"{sample_name}\t{gt_clean}\n")

            if not found:
                out.write("None\n")

            out.write("\n" + "-"*50 + "\n\n")

    print(f"Done. Output written to {output_txt}")

extract_all_aberrant("LDH_END_PILO.tsv", "end-pleistocene.vcf", "LDH_nonref.txt")

from Bio import SeqIO

def load_transcript_lengths(fasta_file):
    lengths = {}
    for rec in SeqIO.parse(fasta_file, "fasta"):
        lengths[rec.id] = len(rec.seq)
    return lengths

def snps_per_transcript(vcf_file):
    counts = {}
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            chrom = line.split("\t")[0]
            counts[chrom] = counts.get(chrom, 0) + 1
    return counts

def normalize_snp_density(snp_counts, lengths):
    density = {}
    for tr, count in snp_counts.items():
        if tr in lengths:
            density[tr] = count / lengths[tr]
        else:
            density[tr] = None
    return density

def compute_percentile_ranks(density, LDH_genes, output_txt):
    
    # All transcript densities (ignore None)
    all_vals = sorted([v for v in density.values() if v is not None])
    total = len(all_vals)

    def percentile(x):
        # count how many transcripts have density <= x
        count = sum(1 for v in all_vals if v <= x)
        return (count / total) * 100

    with open(output_txt, "w") as out:
        out.write("Gene\tTranscript\tDensity\tPercentile\n")

        for gene, tr in LDH_genes.items():
            d = density.get(tr, None)
            if d is None:
                out.write(f"{gene}\t{tr}\tNA\tNA\n")
            else:
                p = percentile(d)
                out.write(f"{gene}\t{tr}\t{d:.6f}\t{p:.2f}\n")

    print(f"LDH percentile ranks written to {output_txt}")

vcf_file = "end-pleistocene.vcf"
fasta_file = "transcripts_end-pleistocene.fa" 

# 1. SNPs per transcript
snp_counts = snps_per_transcript(vcf_file)

# 2. Transcript lengths
lengths = load_transcript_lengths(fasta_file)

# 3. SNP density = SNPs / length
density = normalize_snp_density(snp_counts, lengths)

# 4. Percentile rank for LDH genes
compute_percentile_ranks(density, ldh_gene_names, "LDH_percentiles_end_plei.txt")
