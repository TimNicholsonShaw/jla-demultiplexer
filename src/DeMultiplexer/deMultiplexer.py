from Tailer import LocalAligner as la
from Bio import SeqIO

def parseFASTQ(read1, read2):
    with open(read1, 'r') as r1, open(read2, 'r') as r2:
        r1 = SeqIO.parse(r1, "fastq")
        r2 = SeqIO.parse(r2, "fastq")
        x = zip(list(r1), list(r2))
        return x
    

def getBarcode(entry, ranmerlen=0, trimlen=0):
    pass

def trimEntry(Entry, length=0, loc="5"):
    pass






if __name__ == "__main__":
    x = parseFASTQ("tests/503701_S3_L001_R1_001.fastq", "tests/503701_S3_L001_R2_001.fastq")

    targ = "TGTGCCATCCTTTTCTTGGGGTTGC"

    for pair in x:
        if str(pair[1].seq.reverse_complement())[:len(targ)] == targ:
            print("woo")