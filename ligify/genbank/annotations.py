before_promoter = "ttatttgtatagttcatccatgccatgtgtaatcccagcagctgttacaaactcaagaaggaccatgtggtctctcttttcgttgggatctttcgaaagggcagattgtgtggacaggtaatggttgtctggtaaaaggacagggccatcgccaattggagtattttgttgataatggtctgctagttgaacgcttccatcttcaatgttgtgtctaattttgaagttaactttgattccattcttttgtttgtctgccatgatgtatacattgtgtgagttatagttgtattccaatttgtgtccaagaatgtttccatcttctttaaaatcaataccttttaactcgattctattaacaagggtatcaccttcaaatttgacttcagcacgtgtcttgtagttcccgtcatctttgaaaaatatagttctttcctgtacataaccttcgggcatggcactcttgaaaaagtcatgctgtttcatatgatctgggtatctcgcaaagcattgaagaccatacgcgaaagtagtgacaagtgttggccatggaacaggtagttttccagtagtgcaaataaatttaagggtaagttttccgtatgttgcatcaccttcaccctctccactgacagaaaatttgtgcccattaacatcaccatctaattcaacaagaattgggacaactccagtgaaaagttcttctcctttactcatATATACCCCCTTATTCTCCCGTATTAAACAAAATTATTTGTAGAGGCCCCATTTCGTCCTTTTGGACTCATCAGGGGTGGTACACACCACCCTATGGGGCT"
promoter_seq = "GTCAATTCCTCCAAACTGATTTGTTATCTAATGAGCATTTGACGGCTAAATTCACCATTACTATAATGAAGTAACCGCTATCACACCCTCATCATAACGGGGGTATGGTGGTTAGCAAACAGACAATTTGGTTCGGCAGGCGAATCAATCGTGCGGATGGAGTTGAAAAA"
before_regulator = "tttcttgcaaaaaaaaaccccgccgaagcggggtttttttttttacaggggAAACGaaaatatatttttcaaaagtatcgTTTACCgctagctcagtcctaggTACAATgagcacagTGGCagctgtcaccggatgtgctttccggtctgatgagtccgtgaggacgaaacagcctctaCAAATAATTTTGTTTAAGGGCCCAAGTTCACTTAAAAAGGAGATCAACAATGAAAGCAATTTTCGTACTGAAACATCTTAATCATGCGCCGGTGTTTTAAA"
regulator = "ATGGACCTTGATCGCTTACAGACCTTCCTTAAAGTCGTACAATACGGATCGTTTCAACAAGTCGCAGCGCGTGAGTACCGCAGTCAACGTACAGTCTCTAAGCAGATGACACAGCTGGAAAATGAGTTAAAAGTAACTTTGTTCGACCGCGGTCAGAACCGCATCCAACTTACGCCGCAGGGCCGTTTATTCTGGGCCTCCGCTCAAGACATCGTGAATAACTACACTACGGCTCTTACGGAATTACGTCAATTTAATGTGCCAACCGACCAAATTCTTCGCGTTGGCTATTTTTCGGCCTTCGAGCAACGTTTATTATTGCCCGCGTTATATGATTTAAAGCAACAACACTCCGAGTTACAGTTAGTAGTGCGCCAGGGCTCTAACGAGCATCTTGCACAGCAAGTTGCCGATGGCTCCCTTGATTTAGCATTATCTATTAACTATGGACGCCCAGCTGTAACCCCTGAATCCCAACTGACCGCTGTACCTATTTATCATAACCAGATGGTGATCGGGGTCTCACGCCTTAACCCATTATCGCGTCTTTCACAGTTGCCGCCATCGGCGTTAGCAACGCTGCCCATTTTATATTACTCCCCTGAGTCTTCAACTTTCTTGCTGGAATCCTTCCTTGCCTCAGCTCCTTTCATCCAGGATTATGAGCAAATTCGCCGTGTTTCTTCCGCCGAGCAAATGCACCTTTTAGTGGCCTTGAATCAAGCTTTAGCCTTCTATCCGGCTGGCTTGGTGCCAACCAAACACGATGAGCAAGTTGCGTATCTTCCGATTACCGACGCCGCTCAGCAGGGCTACGACATTGTCGCCTTGTTGAAATCTAACCGTAGTCGCCCTCTGATCGCCAAGTTGGTCCAGCGCTTGAAGGCTAATGCCAAATCCGATGAATAA"
after_regulator = "TAATCACTTTCAGCCAAAAAACTTAAGACCGCCGGTCTTGTCCACTACCTTGCAGTAATGCGGTGGACAGGATCGGCGGTTTTCTTTTCTCTTCTCAAaggctaggtggaggctcagtgatgataagtctgcgatggtggatgcatgtgtcatggtcatagctgtttcctgtgtgaaattgttatccgctcagagggcacaatcctattccgcgctatccgacaatctccaagacattaggtggagttcagttcggcgagcggaaatggcttacgaacggggcggagatttcctggaagatgccaggaagatacttaacagggaagtgagagggccgcggcaaagccgtttttccataggctccgcccccctgacaagcatcacgaaatctgacgctcaaatcagtggtggcgaaacccgacaggactataaagataccaggcgtttccccctggcggctccctcgtgcgctctcctgttcctgcctttcggtttaccggtgtcattccgctgttatggccgcgtttgtctcattccacgcctgacactcagttccgggtaggcagttcgctccaagctggactgtatgcacgaaccccccgttcagtccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggaaagacatgcaaaagcaccactggcagcagccactggtaattgatttagaggagttagtcttgaagtcatgcgccggttaaggctaaactgaaaggacaagttttggtgactgcgctcctccaagccagttacctcggttcaaagagttggtagctcagagaaccttcgaaaaaccgccctgcaaggcggttttttcgttttcagagcaagagattacgcgcagaccaaaacgatctcaagaagatcatcttattaagtctgacgctctattcaacaaagccgccgtccatgggtagggggcttcaaatcgtccgctctgccagtgttacaaccaattaacaaattctgattagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagcttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagcgagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgttttcccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaatcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctcgagcaagacgtttcccgttgaatatggctcataacaccccttgtattactgtttatgtaagcagacagttttattgttcatgatgatatatttttatcttgtgcaatgtaacatcagagattttgagacacaacgtggctttcccccgccgctctagaactagtggatccaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcgcattatacgagacgtccaggttgggatacctgaaacaaaacccatcgtacggccaaggaagtctccaataactgtgatccaccacaagcgccagggttttcccagtcacgacgttgtaaaacgacggccagtcatgcataatccgcacgcatctggaataaggaagtgccattccgcctgacctggaccaaaacgaaaaaaggggagcggtttcccgctcccctcttttctggaatttggtaccgagGAATGAAgCAGGATta"
seq = before_promoter+promoter_seq+before_regulator+regulator+after_regulator

def get_plasmid_components():
    plasmid_components = {
        "before_promoter": before_promoter,
        "before_regulator": before_regulator,
        "after_regulator": after_regulator,
    }
    return plasmid_components


def get_annotations(promoter_seq:str, regulator_seq:str, regulator_name:str, regulator_CDS:str):
    
    prom_len = len(promoter_seq)
    reg_len = len(regulator_seq)

    annotations = [
        {"type": "misc_feature",
        "label": "GFP_mut2",
        "color": "#95ff7d",
        "start": 0,
        "end": 717,
        "strand": -1},

        {"type": "misc_feature",
        "label": "GFP_RBS8",
        "color": "#ff7c6b",
        "start": 717,
        "end": 740,
        "strand": -1},

        {"type": "misc_feature",
        "label": "ElvJ",
        "color": "#cd8cff", 
        "start": 740,
        "end": 818,
        "strand": -1},

        {"type": "misc_feature",
        "label": regulator_name+"_promoter",
        "color": "#faff66",
        "start": 818,
        "end": 818+prom_len,
        "strand": -1},

        {"type": "misc_feature",
        "label": "Bidirectional_terminator",
        "color": "#9e9e9e",
        "start": 818+prom_len,
        "end": 818+prom_len+49,
        "strand": -1},

        {"type": "misc_feature",
        "label": "P250_promoter",
        "color": "#ffae2b",
        "start": 818+prom_len+52,
        "end": 818+prom_len+56+65,
        "strand": 1},

        {"type": "misc_feature",
        "label": "RiboJ",
        "color": "#cd8cff",
        "start": 818+prom_len+56+65,
        "end": 818+prom_len+56+65+75,
        "strand": 1},

        {"type": "misc_feature",
        "label": "5'-UTR",
        "color": "#ffffff",
        "start": 818+prom_len+56+65+75,
        "end": 818+prom_len+56+65+75+32,
        "strand": 1},

        {"type": "misc_feature",
        "label": "4.6k_Leader_peptide",
        "color": "#8589ff",
        "start": 818+prom_len+56+65+75+32,
        "end": 818+prom_len+56+65+75+32+51,
        "strand": 1},

        {"type": "misc_feature",
        "label": regulator_name,
        "color": "#f9abff",
        "start": 818+prom_len+56+65+75+32+52,
        "end": 818+prom_len+56+65+75+32+52+reg_len,
        "strand": 1},

        {"type": "CDS",
        "label": regulator_name+" protein",
        "color": "#f9abff",
        "start": 818+prom_len+56+65+75+32+52,
        "end": 818+prom_len+56+65+75+32+52+reg_len,
        "strand": 1,
        "translation": regulator_CDS},

        {"type": "misc_feature",
        "label": "ECK12_spy_terminator",
        "color": "#9e9e9e",
        "start": 818+prom_len+56+65+75+32+52+reg_len+8,
        "end": 818+prom_len+56+65+75+32+52+reg_len+8+90,
        "strand": 1},

        {"type": "misc_feature",
        "label": "p15A_origin",
        "color": "#abfffe",
        "start": 818+prom_len+56+65+75+32+52+reg_len+8+90+252,
        "end": 818+prom_len+56+65+75+32+52+reg_len+9+90+252+546,
        "strand": -1},

        {"type": "misc_feature",
        "label": "Kanamycin_resistance",
        "color": "#7bad89",
        "start": 818+prom_len+56+65+75+32+52+reg_len+9+90+252+546+110,
        "end": 818+prom_len+56+65+75+32+52+reg_len+9+90+252+546+110+816,
        "strand": -1},

        {"type": "misc_feature",
        "label": "L3S2P00_terminator",
        "color": "#9e9e9e",
        "start": 818+prom_len+56+65+75+32+52+reg_len+9+90+252+546+110+816+381,
        "end": 818+prom_len+56+65+75+32+52+reg_len+9+90+252+546+110+816+381+63,
        "strand": -1},
    ]
    return annotations