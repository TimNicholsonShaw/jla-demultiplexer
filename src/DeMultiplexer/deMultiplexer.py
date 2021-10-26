from Bio import SeqIO, SeqRecord
from Affirmations import affirm
from distance import hamming
import argparse, tempfile, os, subprocess
import pandas as pd
from Tailer import LocalAligner as la

class ReadPair():
    """Object to hold readPairs and associated methods"""

    def __init__(self, r1:SeqRecord.SeqRecord, r2:SeqRecord.SeqRecord):
        self.r1 = r1
        self.r2 = r2

    def getBarcode(self, barLen, revComp=False, read=2):
        """Looks for barcodes at 5' ends of reads. Defaults to reverse complement of read 2"""
        if read==1:
            if revComp:
                return self.r1RevComp()[:barLen].seq
            else:
                return self.r1[:barLen].seq

        elif read==2:
            if revComp:
                return self.r2RevComp()[:barLen].seq
            else:
                return self.r2[:barLen].seq

    def trim(self, length, read, end):
        if read == 1: # read 1 is reversed in illumina sequencing
            if end == 3: # So we trim from the front for 3'
                return ReadPair(self.r1[length:], self.r2)
            elif end == 5: # And from the back half here
                return ReadPair(self.r1[:-length], self.r2)
        elif read == 2: # read 2 is in correct orientation
            if end == 3: # So we trim from the front for 3'
                return ReadPair(self.r1, self.r2[:-length])
            elif end == 5: # And from the back half here
                return ReadPair(self.r1, self.r2[length:])
    def getRanMer(self, ranMerLen):
        return self.r1RevComp()[-ranMerLen:].seq
    def getAnalysisSeq(self, ranMerLen, additionalLigation):
        return self.r1RevComp()[-(ranMerLen+additionalLigation):].seq
    def r1RevComp(self):
        return self.r1.reverse_complement()
    def r2RevComp(self):
        return self.r2.reverse_complement()
class Experiment():
    """Singular experiment consisting of matched r1 and r2"""
    def __init__(self, *args, **kwargs):
        if type(args[0]) == str:
            self.experimentReadPairs = []
            with open(args[0], 'r') as r1, open(args[1], "r") as r2:
                r1 = list(SeqIO.parse(r1, "fastq"))
                r2 = list(SeqIO.parse(r2, "fastq"))
                self.dupRemovedFlag = False

                assert len(r1) == len(r2) # everything should have a pair

                for i in range(len(r1)):
                    self.experimentReadPairs.append(ReadPair(r1[i], r2[i]))
        else: # Allows for creation with just the experiment readpairs
            self.experimentReadPairs = args[0]

    def __iter__(self):
        for read in self.experimentReadPairs:
            yield read

    def __len__(self):
        return(len(self.experimentReadPairs))

    def filterBarcodedPairs(self, barcode, revComp=False, read=2):
        for pair in self.experimentReadPairs:
            if pair.getBarcode(len(barcode), revComp=revComp, read=read) == barcode:
                yield pair

    def filterAndRemovePCRDuplicatesandTrim(self, ranmerLen, barcode, maxHam=1):
        """Very specific to our sequencing methodology.
        Tims adapter"""
        dupDict ={}
        dupRemovedReadPairs = []
        for pair in self.filterBarcodedPairs(barcode):
            seq = pair.getAnalysisSeq(ranmerLen, 2)
            ranmer = pair.getRanMer(ranmerLen)

            previouslySeenRanMers = dupDict.get(seq, [])

            if not previouslySeenRanMers:
                pass
            elif min([hamming(ranmer, x) for x in previouslySeenRanMers]) <= maxHam: 
                continue


            previouslySeenRanMers.append(ranmer)
            trimmed = pair.trim(len(barcode), read=2, end=5)
            trimmed = trimmed.trim(ranmerLen+2, read=1, end=3) # For the AG
            trimmed = trimmed.trim(len(barcode), read=1, end=5) # In case reads to barcode, probably too much
            dupRemovedReadPairs.append(trimmed)
            dupDict[seq] = previouslySeenRanMers
        
        return Experiment(dupRemovedReadPairs)
    def toCSV(self, r1_out, r2_out):
        with open(r1_out, 'w') as r1out, open(r2_out, 'w') as r2out:
            for pair in self:
                r1out.write(pair.r1.format('fastq'))
                r2out.write(pair.r2.format('fastq'))
    def writeRead1(self, outfile):
        with open(outfile, 'w') as file:
            for pair in self:
                file.write(pair.r1.format('fastq'))

    def trimAndDeDuplicate(self, ranmerLen, maxHam=1):
        """Very specific to our sequencing methodology.
        Tims adapter"""
        dupDict ={}
        dupRemovedReadPairs = []
        for pair in self:
            seq = pair.getAnalysisSeq(ranmerLen, 2)
            ranmer = pair.getRanMer(ranmerLen)

            previouslySeenRanMers = dupDict.get(seq, [])

            if not previouslySeenRanMers:
                pass
            elif min([hamming(ranmer, x) for x in previouslySeenRanMers]) <= maxHam: 
                continue


            previouslySeenRanMers.append(ranmer)
            trimmed = pair.trim(ranmerLen+2, read=1, end=3) # For the AG
            dupRemovedReadPairs.append(trimmed)
            dupDict[seq] = previouslySeenRanMers
        
        return Experiment(dupRemovedReadPairs)

class ManifestEntry():
    def __init__(self, csvline):
        "Uses a csv/xls as an input to log the different experimental conditions"
        self.ID = csvline[0].strip()
        self.barcode= csvline[2].strip() + csvline[3].strip() # combines barcode with gene-specific portion 
        self.ranmer = int(csvline[4])
class Manifest():
    def __init__(self, file):
        self.entries =[]
        with open(file, 'r') as infile:
            try:
                if file.endswith(".csv"):
                    next(infile) # skips header
                    infile = [x.split(",") for x in infile]
                else:  # excel file handler
                    infile = pd.ExcelFile(file)
                    infile = infile.parse()
                    infile = [row.values for _ , row in infile.iterrows()]
            except:
                raise("file type not supported")
            
            for line in infile: # works with either csv or excel after above
                if not line[0]: continue # handles emptyp lines
                self.entries.append(ManifestEntry(line)) # creates a new manifest object
    def __iter__(self):
        for entry in self.entries:
            yield entry


@affirm(0.1)       
def fastqBreakDown():
    parser = argparse.ArgumentParser(description="Specific to JLA gene specific sequencing. Removes PCR duplicates, trims barcodes, outputs fastq")
    
    parser.add_argument("-r1", "--read1", help="Location to read1, fastq format, no compression", required=True)
    parser.add_argument("-r2", "--read2", help="Location to read2, fastq format, no compression", required=True)
    parser.add_argument("-b", "--barcode", help="barcode plus gene specific portion, no common adapter", required=True)
    parser.add_argument("-r", "--ranmerlen", type=int, help="Length of random-mer. 10 or 11", required=True)

    args=parser.parse_args()

    print("Reading in Files...")
    experiment = Experiment(args.read1, args.read2)
    print("Filtering/Trimming...")
    experiment = experiment.filterAndRemovePCRDuplicatesandTrim(args.ranmerlen, args.barcode)
    print("Writing to file...")
    experiment.toCSV(args.read1+".processed.fastq", args.read2+".processed.fastq")
    print("Done")

@affirm(0.3)
def manifestAlign():
    parser = argparse.ArgumentParser(description="")
    group = parser.add_mutually_exclusive_group(required=True)

    parser.add_argument("-r1", "--read1", help="Location to read1, fastq format, no compression", required=True)
    parser.add_argument("-r2", "--read2", help="Location to read1, fastq format, no compression", required=True)
    parser.add_argument("-m", "--manifest", help="Location of manifest file", required=True)
    parser.add_argument("-t", "--trimmed", action="store_true", help="Set if you've already trimmed the randommer off")

    group.add_argument("-e", "--ensids", help="Ensembl IDs for database generation comma separated, this or -f is required")
    group.add_argument("-f", "--fasta", help="FASTA to use as reference database. This or -e is required")


    args = parser.parse_args()

    tempDir = tempfile.TemporaryDirectory() # temporary file holder
    queryFile = tempDir.name + "/query.fasta"
    dbFile = tempDir.name + "/db.fa"
    preQueryFile = tempDir.name + "/preQuery.fasta"
    XML_file = tempDir.name + "/temp.xml"

    readPairs = Experiment(args.read1, args.read2) # read fastq files into handler
    manifest = Manifest(args.manifest) # get handled manifest

    args.ensids = args.ensids.split(",")


    pre, ext = os.path.splitext(args.read1) # for giving good file names

    print("Making Database")
    if args.ensids:
        la.buildDBFromEID(args.ensids, dbFile)
    else:
        subprocess.call(["makeblastdb", "-in", args.fasta, "-dbtype", "nucl"])
        dbFile = args.fasta

    for entry in manifest:
        tempReadPairs = readPairs.filterAndRemovePCRDuplicatesandTrim(entry.ranmer, entry.barcode)
        if len(tempReadPairs) == 0:
            print("Nothing found for ", entry.ID)
            continue
        tempReadPairs.writeRead1(preQueryFile)

        print("Parsing...", entry.ID)
        tempReadPairs = la.parseFASTQ(preQueryFile)
        la.queryFormatter(tempReadPairs, queryFile, rev_comp=True)
        la.alignBlastDB(queryFile, dbFile, XML_file)
        tempReadPairs = la.BlastResultsParser(XML_file, tempReadPairs)

        la.tailbuildr(tempReadPairs, pre+entry.ID+"_tails.csv")

@affirm(0.5)
def trimDeDup():
    parser = argparse.ArgumentParser(description="Specific to JLA gene specific sequencing. Removes PCR duplicates, trims randomers, outputs fastq")
    
    parser.add_argument("-r1", "--read1", help="Location to read1, fastq format, no compression", required=True)
    parser.add_argument("-r2", "--read2", help="Location to read2, fastq format, no compression", required=True)
    parser.add_argument("-r", "--ranmerlen", type=int, help="Length of random-mer. 10 or 11", required=True)

    args=parser.parse_args()

    print("Reading in Files...")
    experiment = Experiment(args.read1, args.read2)
    print("Filtering/Trimming...")
    experiment = experiment.trimAndDeDuplicate(args.ranmerlen)
    experiment.toCSV(args.read1+".trimmed.fastq", args.read2+".trimmed.fastq")
    print("Done")



if __name__ == "__main__":
    trimDeDup()


