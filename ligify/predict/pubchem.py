import requests

def get_inchiKey(input, prop):

    if prop == "name":
        URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+str(input)+"/property/InChiKey/TXT"

        response = requests.get(URL)
        if response.ok:
                # get the first entry
            out = response.text.split("\n")[0]
            return out

    elif prop == "smiles":
        URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"+str(input)+"/property/InChiKey/TXT"

        response = requests.get(URL)
        if response.ok:
                # get the first entry
            out = response.text.split("\n")[0]
            return out
     

def get_smiles(input):

    URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"+str(input)+"/property/CanonicalSMILES/TXT"

    response = requests.get(URL)
    if response.ok:
            # get the first entry
        out = response.text.split("\n")[0]
        return out


def check_url(url):

    response = requests.get(url)
    if response.ok:
        return True
    else:
        return False

if __name__ == "__main__":
    # out = get_inchiKey("isovalerate", "name")
    # print(out)

    # out = get_inchiKey("C=CC(=O)[O-]", "smiles")
    # print(out)

    # out = get_smiles("isovalerate")
    # print(out)

    out = check_url("http://hulab.rxnfinder.org/smi2img/3")
    print(out)