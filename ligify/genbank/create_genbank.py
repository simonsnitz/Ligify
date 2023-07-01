from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from annotations import get_plasmid_components, get_annotations


regulator_name = "VprR"
ligand_name = "4-ethylphenol"
promoter_seq = "GTCAATTCCTCCAAACTGATTTGTTATCTAATGAGCATTTGACGGCTAAATTCACCATTACTATAATGAAGTAACCGCTATCACACCCTCATCATAACGGGGGTATGGTGGTTAGCAAACAGACAATTTGGTTCGGCAGGCGAATCAATCGTGCGGATGGAGTTGAAAAA"
regulator_seq = "ATGGACCTTGATCGCTTACAGACCTTCCTTAAAGTCGTACAATACGGATCGTTTCAACAAGTCGCAGCGCGTGAGTACCGCAGTCAACGTACAGTCTCTAAGCAGATGACACAGCTGGAAAATGAGTTAAAAGTAACTTTGTTCGACCGCGGTCAGAACCGCATCCAACTTACGCCGCAGGGCCGTTTATTCTGGGCCTCCGCTCAAGACATCGTGAATAACTACACTACGGCTCTTACGGAATTACGTCAATTTAATGTGCCAACCGACCAAATTCTTCGCGTTGGCTATTTTTCGGCCTTCGAGCAACGTTTATTATTGCCCGCGTTATATGATTTAAAGCAACAACACTCCGAGTTACAGTTAGTAGTGCGCCAGGGCTCTAACGAGCATCTTGCACAGCAAGTTGCCGATGGCTCCCTTGATTTAGCATTATCTATTAACTATGGACGCCCAGCTGTAACCCCTGAATCCCAACTGACCGCTGTACCTATTTATCATAACCAGATGGTGATCGGGGTCTCACGCCTTAACCCATTATCGCGTCTTTCACAGTTGCCGCCATCGGCGTTAGCAACGCTGCCCATTTTATATTACTCCCCTGAGTCTTCAACTTTCTTGCTGGAATCCTTCCTTGCCTCAGCTCCTTTCATCCAGGATTATGAGCAAATTCGCCGTGTTTCTTCCGCCGAGCAAATGCACCTTTTAGTGGCCTTGAATCAAGCTTTAGCCTTCTATCCGGCTGGCTTGGTGCCAACCAAACACGATGAGCAAGTTGCGTATCTTCCGATTACCGACGCCGCTCAGCAGGGCTACGACATTGTCGCCTTGTTGAAATCTAACCGTAGTCGCCCTCTGATCGCCAAGTTGGTCCAGCGCTTGAAGGCTAATGCCAAATCCGATGAATAA"
regulator_CDS = "MDLDRLQTFLKVVQYGSFQQVAAREYRSQRTVSKQMTQLENELKVTLFDRGQNRIQLTPQGRLFWASAQDIVNNYTTALTELRQFNVPTDQILRVGYFSAFEQRLLLPALYDLKQQHSELQLVVRQGSNEHLAQQVADGSLDLALSINYGRPAVTPESQLTAVPIYHNQMVIGVSRLNPLSRLSQLPPSALATLPILYYSPESSTFLLESFLASAPFIQDYEQIRRVSSAEQMHLLVALNQALAFYPAGLVPTKHDEQVAYLPITDAAQQGYDIVALLKSNRSRPLIAKLVQRLKANAKSDE*"

plasmid_components = get_plasmid_components()
seq = plasmid_components["before_promoter"] + \
        promoter_seq + \
        plasmid_components["before_regulator"] + \
        regulator_seq + \
        plasmid_components["after_regulator"]

# Create the plasmid sequence
plasmid_sequence = Seq(seq)
record = SeqRecord(plasmid_sequence,
                   id='123456789', # random accession number
                   name="pLigify_"+regulator_name,
                   description='This is a genetic circuit designed by Ligify to express GFP in response to the ligand '+ligand_name+' using the regulator '+regulator_name,
                   annotations={"molecule_type": "DNA"})

# Set topology as circular
record.annotations["topology"] = "circular"

# Add annotations
annotations = get_annotations(promoter_seq, regulator_seq, regulator_name, regulator_CDS)
for annotation in annotations:
    if "translation" in annotation.keys():
        record.features.append(\
            SeqFeature(FeatureLocation(\
                start = annotation["start"], \
                end = annotation["end"], \
                strand = annotation["strand"]), \
                type = annotation["type"], \
                qualifiers={ \
                    "ApEinfo_fwdcolor":[annotation["color"]],\
                    "translation":[annotation["translation"]],\
                    "label": annotation["label"]}))
    else:
        record.features.append(\
            SeqFeature(FeatureLocation(\
                start = annotation["start"], \
                end = annotation["end"], \
                strand = annotation["strand"]), \
                type = annotation["type"], \
                qualifiers={ \
                    "ApEinfo_fwdcolor":[annotation["color"]],\
                    "label": annotation["label"]}))


# Save as GenBank file
output_file = open('pLigify_'+regulator_name+'.gb', 'w')
SeqIO.write(record, output_file, 'genbank')
