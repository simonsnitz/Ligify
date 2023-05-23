import requests
from urllib.error import HTTPError


    # Input must be KEGG compound ID. (from TFBMiner)
    ### NOT using this. EC numbers returned by KEGG are loose-fitting for the input chemical
        # For example, for acrylate, EC numbers corresponding to a broad enzymes class acting on a carboxylic acid is returned.
def fetch_kegg_reaction(compound):

    # Removes coefficients that obfuscate compound IDs.
    if " " in compound:
        compound = compound.split(" ")[1]
    search_term = "cpd:" + str(compound)

    data = requests.get(f"http://rest.kegg.jp/get/{search_term}")

    if data is None:
        return []
    else:
        try:
            # Retrieved HTML data is read line by line.
            current_section = None
            for line in data.text.strip().split("\n"):
                section = line[:12].strip()
                if not section == "":
                    current_section = section
                    # Finds section where reactions are listed
                    # and extracts their IDs from each line.
                    if current_section == "ENZYME":
                        index = data.text.strip().split("\n").index(line)
                        reactions = line[12:].split(" ")
                        for line in data.text.strip().split("\n")[int(index)+1:]:
                            # If the current line is not part of the REACTION 
                            # section then the processing ends.
                            if line[:12].strip() != (""):
                                break
                            reactions_ = line[12:].split(" ")
                            reactions.extend(reactions_)
            # Eliminates erroneous IDs resulting from whitespaces.
            reactions = [reaction for reaction in reactions if reaction != ""]
            
            return reactions
        except HTTPError:
            return None
