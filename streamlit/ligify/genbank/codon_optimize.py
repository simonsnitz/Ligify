from dnachisel import *


def codon_opt(protein_seq:str):

    # Create a random DNA seq given the protein seq. Append a stop codon.
    protein_dna_seq = reverse_translate(protein_seq+"*")

    # DEFINE THE OPTIMIZATION PROBLEM
    problem = DnaOptimizationProblem(
        sequence=protein_dna_seq,
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceGCContent(mini=0.35, maxi=0.65, window=50),
            EnforceTranslation(location=(0, len(protein_dna_seq)))
        ],
        objectives=[CodonOptimize(species='e_coli', location=(0, len(protein_dna_seq)))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

    final_sequence = problem.sequence  # string
    return final_sequence


if __name__ == "__main__":
    codon_opt("MTTIRWRRMSIHSERITLADSPLHWAHTLNGSMRTHFEVQRLERGRGAYLARSRFGAGELYSAIAPSQVLRHFNDQRNANEAEHSYLIQIRSGALGVASGGRKVILANGDCSIVDSRQDFTLSSNSSTQGVVIRFPVSWLGAWVSNPEDLIARRVDAEIGWGRALSASVSNLDPLRIDDLGSNVNSIAEHVAMLISLASSAVSSEDGGVALRKMREVKRVLEQSFADANLEPESVSSQLGISKRYLHYVFAACGTTFGRELLEIRLGKAYRMLCATSGSGAVLKVAMSSGFSDSSHFSKKFKERYGVSPVSLVRQA")

