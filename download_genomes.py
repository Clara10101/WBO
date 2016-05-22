from Bio import Entrez
from Bio import SeqIO

Entrez.email = "*****@gmail.com"

#identyfikatory wczesniej wybranych genomow
genomeAccessions = ['NC_009925.1','NZ_CP011120.1','NZ_CP011389.1','NZ_LN609302.1','NC_013209.1']
"""['CP004420.1','CP004624.1','CP006341.2','CP004670.1','CP004764.1','CP004900.1','CP005204.1','CP006129.1','CP004992.1','CP005124.1','CP005301.1','CP006445.1','CP005402.1','CP005504.1','CP005598.1','CP006221.1','CP006484.1','CP004510.1',
'CM003553.1','CM003554.1','CM003556.1','CM003555.1','CM003557.1','CM003558.1','CM003559.1','CM003560.1','CM003561.1','CM003562.1','CM003563.1','CM003564.1','CM003565.1','CM003566.1','CM003567.1','CM003568.1','CM003569.1','CM003570.1',
'CM003579.1','CM003588.1','CM003590.1','CM003576.1','CM003586.1','CM003584.1','CM003587.1','CM003589.1','CM003577.1','CM003585.1','CM003582.1','CM003591.1','CM003578.1','CM003581.1','CM003580.1','CM003583.1','CM003593.1',
'CM001563.1','CM001564.1','CM001565.1','CM001566.1','CM001567.1','CM001568.1','CM001569.1','CM001570.1','CM001571.1','CM001572.1','CM001573.1','CM001574.1','CM001575.1','CM001576.1','CM001577.1','CM001578.1','CM001579.1',
'NC_001133.9','NC_001134.8','NC_001135.5','NC_001136.10','NC_001137.3','NC_001138.5','NC_001139.9','NC_001140.6','NC_001141.2','NC_001142.9','NC_001143.9','NC_001144.5','NC_001145.3','NC_001146.8','NC_001147.6','NC_001148.4','NC_001224.1',
'NC_015518.1','NC_015435.1','NZ_BBBY01000001.1','NC_003106.2','NC_007181.1',
'NC_021353.1', 'NC_020913.1','NZ_CAJE01000012.1','NC_002578.1','NC_005877.1']
"""

handle = Entrez.efetch(db="nucleotide", id=genomeAccessions, rettype="gbwithparts")
records = SeqIO.parse(handle, "gb")

for i,record in enumerate(records):
    output_file = open(record.id + ".fasta", "w")

    for feature in record.features:
        if feature.type == "CDS":
            if "protein_id" in feature.qualifiers and "translation" in feature.qualifiers:
                output_file.write(">" + feature.qualifiers["protein_id"][0]+"\n")
                output_file.write(feature.qualifiers["translation"][0]+"\n")