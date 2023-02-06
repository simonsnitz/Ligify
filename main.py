from pprint import pprint

from src.chemical2enzymes import chem2enzymes
from src.enzymes2operons import append_operons, pull_regulators, protein2chemicals
from src.make_chemDB import make_chemDB
from src.query_chemDB import query_chemDB


# InChiKey = "NIXOWILDQLNWCW-UHFFFAOYSA-M"    # Found AcuR
# chemical_name = "acrylate"

# InChiKey = "KWIUHFFTVRNATP-UHFFFAOYSA-N"     # Found BetI
# chemical_name = "choline"

# InChiKey = "ZHDLAGPONFNQMZ-UHFFFAOYSA-M"     # Found CymR
# chemical_name = "cumate"

# InChiKey = "NWXMGUDVXFXRIG-WESIUVDSSA-N"     # Nothing returned ...
# chemical_name = "tetracycline"

# InChiKey = "ISAKRJDGNUQOIC-UHFFFAOYSA-N"     # Found RutR
# chemical_name = "uracil"

# InChiKey = "DSSYKIVIOFKYAU-XCBNKYQSSA-N"     # Did NOT find CamR, but found CamA & CamC.
# chemical_name = "d-camphor"

# InChiKey = "WKOLLVMJNQIZCI-UHFFFAOYSA-N"     # Did NOT find DesR or VanR
# chemical_name = "vanillic_acid"

# InChiKey = "WTDRDQBEARUVNC-LURJTMIESA-N"     # Nothing returned
# chemical_name = "l-dopa"

# InChiKey = "VYFYYTLLBUKUHU-UHFFFAOYSA-N"     # Found two really interesting candidates! RDC20414.1 & RDC20415.1
# chemical_name = "dopamine"

# InChiKey = "LRFVTYWOQMYALW-UHFFFAOYSA-N"     # Didn't find VvPecS, but found several others
# chemical_name = "xanthine"

# InChiKey = "CZMRCDWAGMRECN-UGDNZRGBSA-N"     # Didn't find SghR, but found other candidates.
# chemical_name = "sucrose"

# InChiKey = "LOIYMIARKYCTBW-OWOJBTEDSA-N"        # Found 85% identical homolog of known regulator (HutC)
# chemical_name = "urconate"

#InChiKey = "CPJRRXSHAYUTGL-UHFFFAOYSA-N"        # Found 
chemical_name = "buten-2-one"
InChiKey = "FUSUHKVFWTUUBE-UHFFFAOYSA-N"
chemical_smiles = "CC(=O)C=C"


    # Rhea enzyme output filter criteria
domain_filter = "Bacteria"
lineage_filter_name = "Family"
reviewed= "false"


    # file where I'm storing info on all chemicals included in Rhea
rhea_chemDB = "data/all_rhea_chemicals.json"
make_chemDB(rhea_chemDB)
scores = query_chemDB(rhea_chemDB, chemical_smiles)
pprint(scores[0:5])


chem2enzymes(InChiKey = InChiKey,
                chemical_name = chemical_name, 
                domain_filter = domain_filter,
                lineage_filter_name = lineage_filter_name, 
                reviewed = reviewed)

append_operons(chemical_name)

pull_regulators(chemical_name)

# protein2chemicals('WP_010889412.1')