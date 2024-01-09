

from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO






def blast_align(sequence1, sequence2):

    # Create two sequence files
    seq1 = SeqRecord(Seq(sequence1),
                    id="seq1")
    seq2 = SeqRecord(Seq(sequence2),
                    id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse the output as XML
    output = NcbiblastpCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]
    blast_result_record = NCBIXML.read(StringIO(output))

    # Print some information on the result
    data = []
    for alignment in blast_result_record.alignments:

        for hsp in alignment.hsps:

            identity = round((hsp.identities/hsp.align_length)*100, 2)

            coverage = round((hsp.align_length/blast_result_record.query_length)*100, 2)

            entry = {
                    "e value": hsp.expect,
                    "identity": identity,
                    "coverage": coverage}
            data.append(entry)

    try:
        best_e = min([i["e value"] for i in data])
        best_align = [i for i in data if i["e value"] == best_e][0]
        return best_align

    except:
        entry = {
        "e value": 1,
        "identity": 0,
        "coverage": 0}
        return entry


    
if __name__ == "__main__":

    seq1 = "MNKGQRHIKIREIITANEIETQDELVDILKKDGYNVTQATVSRDIKELHLVKVPTNNGSYKYSLPADQRFNPLSKLKRSLMDAFVKIDSASHLIVLKTMPGNAQAIGALMDNLDWEEIMGTICGDDTILIICRTHDDTKVVQKKILELL"
    seq2 = "MSNLEKSLRHEVILDIIESKCICKQEELIIELKECGINVTQATLSRDLHEMNIIRKSFENHEHRYIVTKEDKFIKKFKNIFKSSVRSISMQEYFISVNTLDGMASLIGEFIDRLGDGRIAGTVAKKNHILVLCRSVNATQQIFKELDSLRV"

    print(blast_align(seq1, seq2))