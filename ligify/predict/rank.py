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
    total_regs = -1
    for gene in operon:
        if "description" in gene.keys():
            if regulator_re.search(gene["description"]):
                total_regs += 1

    # Calculate final rank score
    rank = 100
        # Subtract 10 points for every gene between the enzyme and the regulator
    distance_deduction = -10*(enz_reg_distance-1)
    rank = rank + distance_deduction

        # Subtract 15 points for every additional regulator in the operon
    extra_reg_deduction = -15*(total_regs)
    rank = rank + extra_reg_deduction
        # Subtract 5 points for every additional gene in the operon
    total_genes_deduction = -5*(total_genes - 2)
    rank = rank + total_genes_deduction

    # Color-code based on rank score
    if rank >= 70:
        color = "#02a602"  #Green
    elif rank >= 50  and rank < 70:
        color = "#d4d400"  # Yellow
    elif rank >= 30 and rank < 50:
        color = "#d48302"  # Orange
    else:
        color = "#f50b02"  # Red


    return {"rank": rank, "color": color, "metrics": { \
            "Genes within operon": {"Value":total_genes, "Deduction": total_genes_deduction}, \
                "Enzyme-regulator distance": {"Value": enz_reg_distance, "Deduction": distance_deduction}, \
                    "Additional regulators": {"Value": total_regs, "Deduction": extra_reg_deduction}, \
                        }
            }




if __name__ == "__main__":

    with open("regulator.json", "r") as f:
        data = json.load(f)

        print(calculate_rank(data))