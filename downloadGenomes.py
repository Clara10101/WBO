from Bio import Entrez
from Bio import SeqIO

Entrez.email = "*****@gmail.com"

#identyfikatory wczesniej wybranych genomow
genomeAccessions = ['NC_009925.1','NZ_CP011120.1','NZ_CP011389.1','NZ_LN609302.1','NC_013209.1']

handle = Entrez.efetch(db="nucleotide", id=genomeAccessions, rettype="gbwithparts")
records = SeqIO.parse(handle, "gb")

for i,record in enumerate(records):
    output_file = open(record.id + ".fasta", "w")

    for feature in record.features:
        if feature.type == "CDS":
            if "protein_id" in feature.qualifiers and "translation" in feature.qualifiers:
                output_file.write(">" + feature.qualifiers["protein_id"][0]+"\n")
                output_file.write(feature.qualifiers["translation"][0]+"\n")
