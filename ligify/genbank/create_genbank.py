from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from tempfile import NamedTemporaryFile

# from annotations import get_plasmid_components, get_annotations
# from codon_optimize import codon_opt

from ligify.genbank.annotations import get_plasmid_components, get_annotations
from ligify.genbank.codon_optimize import codon_opt

regulator_name = "VprR"
ligand_name = "4-ethylphenol"
promoter_seq = "GTCAATTCCTCCAAACTGATTTGTTATCTAATGAGCATTTGACGGCTAAATTCACCATTACTATAATGAAGTAACCGCTATCACACCCTCATCATAACGGGGGTATGGTGGTTAGCAAACAGACAATTTGGTTCGGCAGGCGAATCAATCGTGCGGATGGAGTTGAAAAA"
regulator_seq = "ATGGACCTTGATCGCTTACAGACCTTCCTTAAAGTCGTACAATACGGATCGTTTCAACAAGTCGCAGCGCGTGAGTACCGCAGTCAACGTACAGTCTCTAAGCAGATGACACAGCTGGAAAATGAGTTAAAAGTAACTTTGTTCGACCGCGGTCAGAACCGCATCCAACTTACGCCGCAGGGCCGTTTATTCTGGGCCTCCGCTCAAGACATCGTGAATAACTACACTACGGCTCTTACGGAATTACGTCAATTTAATGTGCCAACCGACCAAATTCTTCGCGTTGGCTATTTTTCGGCCTTCGAGCAACGTTTATTATTGCCCGCGTTATATGATTTAAAGCAACAACACTCCGAGTTACAGTTAGTAGTGCGCCAGGGCTCTAACGAGCATCTTGCACAGCAAGTTGCCGATGGCTCCCTTGATTTAGCATTATCTATTAACTATGGACGCCCAGCTGTAACCCCTGAATCCCAACTGACCGCTGTACCTATTTATCATAACCAGATGGTGATCGGGGTCTCACGCCTTAACCCATTATCGCGTCTTTCACAGTTGCCGCCATCGGCGTTAGCAACGCTGCCCATTTTATATTACTCCCCTGAGTCTTCAACTTTCTTGCTGGAATCCTTCCTTGCCTCAGCTCCTTTCATCCAGGATTATGAGCAAATTCGCCGTGTTTCTTCCGCCGAGCAAATGCACCTTTTAGTGGCCTTGAATCAAGCTTTAGCCTTCTATCCGGCTGGCTTGGTGCCAACCAAACACGATGAGCAAGTTGCGTATCTTCCGATTACCGACGCCGCTCAGCAGGGCTACGACATTGTCGCCTTGTTGAAATCTAACCGTAGTCGCCCTCTGATCGCCAAGTTGGTCCAGCGCTTGAAGGCTAATGCCAAATCCGATGAATAA"



def create_genbank(regulator_name, ligand_name, promoter_seq, regulator_protein_seq):
    
    # Codon optimize the natural sequence
    opt_regulator_seq = codon_opt(regulator_protein_seq)

    # Create the plasmid sequence record
    plasmid_components = get_plasmid_components()
    seq = plasmid_components["before_promoter"] + \
            promoter_seq + \
            plasmid_components["before_regulator"] + \
            opt_regulator_seq + \
            plasmid_components["after_regulator"]
    plasmid_sequence = Seq(seq)
    record = SeqRecord(plasmid_sequence,
                    id=str(regulator_name),
                    name="pLigify_"+regulator_name,
                    description='This is a genetic circuit designed by Ligify to express GFP in response to the ligand '+ligand_name+' using the regulator '+regulator_name,
                    annotations={"molecule_type": "DNA"})

    # Set topology as circular
    record.annotations["topology"] = "circular"

    # Add annotations
    annotations = get_annotations(promoter_seq, opt_regulator_seq, regulator_name, regulator_protein_seq)
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


    # Converts gbk into straight text
    def get_gbk(record):
        outfileloc=NamedTemporaryFile()
        with open(outfileloc.name, "w") as handle:
            SeqIO.write(record, handle, "genbank")
        with open(outfileloc.name) as handle:
            record=handle.read()
        outfileloc.close()

        return record

    record = get_gbk(record)
    
    return record

if __name__ == "__main__":

    regulator_seq = "ATGGACCTTGATCGCTTACAGACCTTCCTTAAAGTCGTACAATACGGATCGTTTCAACAAGTCGCAGCGCGTGAGTACCGCAGTCAACGTACAGTCTCTAAGCAGATGACACAGCTGGAAAATGAGTTAAAAGTAACTTTGTTCGACCGCGGTCAGAACCGCATCCAACTTACGCCGCAGGGCCGTTTATTCTGGGCCTCCGCTCAAGACATCGTGAATAACTACACTACGGCTCTTACGGAATTACGTCAATTTAATGTGCCAACCGACCAAATTCTTCGCGTTGGCTATTTTTCGGCCTTCGAGCAACGTTTATTATTGCCCGCGTTATATGATTTAAAGCAACAACACTCCGAGTTACAGTTAGTAGTGCGCCAGGGCTCTAACGAGCATCTTGCACAGCAAGTTGCCGATGGCTCCCTTGATTTAGCATTATCTATTAACTATGGACGCCCAGCTGTAACCCCTGAATCCCAACTGACCGCTGTACCTATTTATCATAACCAGATGGTGATCGGGGTCTCACGCCTTAACCCATTATCGCGTCTTTCACAGTTGCCGCCATCGGCGTTAGCAACGCTGCCCATTTTATATTACTCCCCTGAGTCTTCAACTTTCTTGCTGGAATCCTTCCTTGCCTCAGCTCCTTTCATCCAGGATTATGAGCAAATTCGCCGTGTTTCTTCCGCCGAGCAAATGCACCTTTTAGTGGCCTTGAATCAAGCTTTAGCCTTCTATCCGGCTGGCTTGGTGCCAACCAAACACGATGAGCAAGTTGCGTATCTTCCGATTACCGACGCCGCTCAGCAGGGCTACGACATTGTCGCCTTGTTGAAATCTAACCGTAGTCGCCCTCTGATCGCCAAGTTGGTCCAGCGCTTGAAGGCTAATGCCAAATCCGATGAATAA"

    protein = Seq(regulator_seq).translate()

    print(protein)
