import json

with open("all_rhea_chemicals.json", "r") as f:
    data = json.load(f)

    names = "".join(i["name"]+"\n" for i in data)

    with open("names", "w+") as r:
        r.write(names)