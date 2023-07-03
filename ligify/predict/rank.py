import json
import re


def calculate_rank(r):

    operon = r["protein"]["context"]["operon"]

    # Calculate the total number of genes in the operon
    total_genes = len(operon)

    # Get the regulator position within the operon
    reg_index = 0
    for gene in operon:
        if gene["accession"] == r["refseq"]:
            break
        else:
            reg_index += 1
    # Get the enzyme position within the operon
    enz_index = r['protein']['context']['enzyme_index']

    # Calculate the distance between the regulator and the ligand-associated enzyme
    enz_reg_distance = (max(enz_index, reg_index)-min(enz_index, reg_index))

    # Calculate how many other regulators are in the operon
    regulator_re = re.compile(r"regulator|repressor|activator")
    total_regs = 0
    for gene in operon:
        if regulator_re.search(gene["description"]):
            total_regs += 1

    # Calculate final rank score
    rank = 100
        # Subtract 10 points for every gene between the enzyme and the regulator
    rank = rank - 10*(enz_reg_distance-1)
        # Subtract 15 points for every additional regulator in the operon
    rank = rank - 15*(total_regs-1)
        # Subtract 5 points for every additional gene in the operon
    rank = rank - 5*(total_genes - 2)

    # Color-code based on rank score
    if rank >= 70:
        color = "#02a602"  #Green
    elif rank >= 50  and rank < 70:
        color = "#d4d400"  # Yellow
    elif rank >= 30 and rank < 50:
        color = "#d48302"  # Orange
    else:
        color = "#f50b02"  # Red


    return {"rank": rank, "total_genes": total_genes, "enz_reg_distance": enz_reg_distance, "total_regs": total_regs, "color": color}




if __name__ == "__main__":

    with open("regulator.json", "r") as f:
        data = json.load(f)

        print(calculate_rank(data))