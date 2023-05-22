import json
from pprint import pprint

with open("dopamine.json", "r") as f:
    data = json.load(f)["rxn_data"][0]["proteins"][0]

    pprint(", ".join(data["organism"])[:-2])