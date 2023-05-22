import json
import requests
from Bio.Blast.Applications import NcbiblastpCommandline


def create_fasta():
    with open("groovSequences.json", "r") as f:
        data = json.load(f)


    fasta = ""
    for i in data:
        header = ">"+str(i["name"]+"\n")
        seq = str(i["sequence"])
        protein = header + seq + "\n\n"
        fasta += protein

    with open("groovSeq.fasta", "w+") as out:
        out.write(fasta)



    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        return None





    # Input protein sequence. Output cached blast results
def blast(acc: str):

    print("NOTE: Starting BLAST")
    blast_db = "blast/groovSeq.fasta"
    num_aligns = 5

    seq = accID2sequence(acc)

    if seq != None:
            # Must have BLAST+ executables in PATH to run this
        blast_cline = NcbiblastpCommandline(db=blast_db, outfmt="6 sseqid pident qcovs", \
            num_alignments=num_aligns)

        results, err = blast_cline(stdin=seq)

        results = results.split("\n")[:-1]

        print(results)
        
        # homologs = [{"accession": r.split("|")[1], \
        #         "identity": r.split("|")[2].split("\t")[1], \
        #         "coverage": r.split("|")[2].split("\t")[2].strip()} \
        #         for r in results ]

            # filter out homologs that don't meet the percent identity cutoff
        # homologs = [h for h in homologs if float(h["identity"]) >= pident_cutoff]

    #     alignment = json.dumps(homologs, indent=4)
    #     return alignment

    # else:
    #     return None


blast("KJF19170.1")
