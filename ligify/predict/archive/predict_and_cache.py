

import sys 
import pandas as pd
import os
import json
from pprint import pprint

from pubchem import get_inchiKey
from chemical2enzymes import chem2enzymes
from enzymes2operons import append_operons, pull_regulators


with open("temp/all.json", "r") as f:
    all_chemicals = json.load(f)

    for chemical_name in all_chemicals:
        InChiKey = get_inchiKey(str(chemical_name), "name")
        domain_filter = "Bacteria"
        lineage_filter_name = "Family"
        reviewed = True

        if os.path.exists("../temp/"+str(chemical_name)+".json"):
            print("data for "+str(chemical_name)+" already exists")

        else:

            data = chem2enzymes(InChiKey = InChiKey,
                domain_filter = domain_filter,
                lineage_filter_name = lineage_filter_name, 
                reviewed_bool = reviewed)
            print("enzyme data collected for "+str(chemical_name))

            data = append_operons(data, chemical_name)

            if data == None:
                print("No regulators found")
            
            else:

                regulators = pull_regulators(data, chemical_name)
                
                with open("../temp/"+str(chemical_name)+".json", "w+") as entry:
                    entry.write(json.dumps(regulators))
                    print("cached regulator data")

                if regulators == None or len(regulators) == 0:
                    print("No regulators found")
                else:
                    print("regulators found for "+str(chemical_name))