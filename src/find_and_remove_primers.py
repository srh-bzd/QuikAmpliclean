#!/usr/bin/env python3

import os, argparse, sys

"""
USAGE
    ./find_and_remove_primers.py [-h] <fastq_file> <cutadapt_info_file> <output_file> [-r <yes|no>] [-i <output_log>] [-c <output_counts>]
DESCRIPTION
    find_and_remove_primers.py is a script to keep sequences where the correct combinations of primers have been found and removed.
PREREQUISITE
    - python3
    - Cutadapt info-file
INFORMATION
    This script is a part of QuikClean
AUTHOR :
    Sarah BOUZIDI
    Engineer in bioinformatics
    Centre National de la Recherche Scientifique (CNRS)
    Laboratory MIVEGEC, IRD, Montpellier
"""


def parse_fastq(input_file):
    """
    Parse input FASTQ file and keep it into a dict like : {seqId : {"seq" : , "qual" : }}
    """
    fastqDict = dict()
    seqId = None
    for line in input_file:
        if line.startswith("@") and "?" not in line: # Rq : some quality lines start with @ this is why there is the second part of the if
            seqId = line.strip().replace("@","").split('\t')[0]
            complementary_lines = []
            fastqDict[seqId] = {}
            fastqDict[seqId]["seq"] = ""
            fastqDict[seqId]["qual"] = ""
        else:
            complementary_lines.append(line.strip())
            if len(complementary_lines)==3:
                fastqDict[seqId]["seq"] += complementary_lines[0] 
                fastqDict[seqId]["qual"] += complementary_lines[2]
    return fastqDict


def parse_cutadapt_log(cutadapt_file):
    """
    Parse input Cutadapt info-file and keep seqId of sequences where primers was find into a dict : {seqId : [primers]}
    """
    seqWrongPrimersDict = {"without":[], "one":[], "wrong":[]}
    seqGoodPrimersDict = dict()
    seqPrimerDict = dict()
    seqId = None
    primer = None
    for line in cutadapt_file:
        seqId = line.strip().split('\t')[0]
        if (int(line.strip().split('\t')[1]) == -1):
            seqWrongPrimersDict["without"].append(seqId)
        else : 
            primer = line.strip().split('\t')[7]
            if (seqId not in seqPrimerDict.keys()) :
                seqPrimerDict[seqId] = [primer]
            else:
                seqPrimerDict[seqId].append(primer)
    for seqId in seqPrimerDict:
        if (len(seqPrimerDict[seqId]) == 2) and (seqPrimerDict[seqId][0].endswith("_R") and seqPrimerDict[seqId][1].endswith("_F_RC") or seqPrimerDict[seqId][0].endswith("_F") and seqPrimerDict[seqId][1].endswith("_R_RC")):
            seqGoodPrimersDict[seqId] = seqPrimerDict[seqId]
        if (len(seqPrimerDict[seqId]) == 2) and not (seqPrimerDict[seqId][0].endswith("_R") and seqPrimerDict[seqId][1].endswith("_F_RC") or seqPrimerDict[seqId][0].endswith("_F") and seqPrimerDict[seqId][1].endswith("_R_RC")):
            seqWrongPrimersDict["wrong"].append(seqId)
        if (len(seqPrimerDict[seqId]) == 1):
            seqWrongPrimersDict["one"].append(seqId)
    return seqGoodPrimersDict, seqWrongPrimersDict


def write_output_files(fastqDict, seqGoodPrimersDict, seqWrongPrimersDict, output_file, reverse_complement, output_log, output_counts):
    """
    Write the new FASTQ file with sequences saved and log files : 1/ with id of sequences not saved, 2/ with counts of sequences not saved
    """
    # File with saved sequences
    for seqId in seqGoodPrimersDict.keys():
        if (seqGoodPrimersDict[seqId][0].endswith("_R") and reverse_complement == "yes"):
            comp = fastqDict[seqId]["seq"].maketrans('ATGC', 'TACG')
            print("@"+seqId, file=output_file)
            print(fastqDict[seqId]["seq"].translate(comp)[::-1], file=output_file)
            print("+", file=output_file)
            print(fastqDict[seqId]["qual"][::-1], file=output_file)
        else:
            print("@"+seqId, file=output_file)
            print(fastqDict[seqId]["seq"], file=output_file)
            print("+", file=output_file)
            print(fastqDict[seqId]["qual"], file=output_file)
    # Log files
    # id
    print("*---------- Sequences where there isn't any primer ("+str(len(seqWrongPrimersDict["without"]))+") : ", file=output_log)
    print("\n".join(seqWrongPrimersDict["without"]), file=output_log)
    print("\n*---------- Sequences where there only one primer ("+str(len(seqWrongPrimersDict["one"]))+") : ", file=output_log)
    print("\n".join(seqWrongPrimersDict["one"]), file=output_log)
    print("\n*---------- Sequences where there wrong primers ("+str(len(seqWrongPrimersDict["wrong"]))+") : ", file=output_log)
    print("\n".join(seqWrongPrimersDict["wrong"]), file=output_log)
    # counts
    outputName = os.path.split(output_file.name)[1] 
    outputName = os.path.splitext(outputName)[0]
    outputName = outputName.replace("_with_primers","")
    print(outputName+"\t"+str(len(fastqDict))+"\t"+str(len(seqGoodPrimersDict))+"\t"+str(len(seqWrongPrimersDict["without"]))+"\t"+str(len(seqWrongPrimersDict["one"]))+"\t"+str(len(seqWrongPrimersDict["wrong"])), file=output_counts)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="find_and_remove_primers.py is a script to keep sequences where the correct combinations of primers have been found and removed.")
    parser.add_argument('fastq_file', help="FASTQ file, not compressed, with sequences;", type=argparse.FileType('r'))
    parser.add_argument('cutadapt_info_file', help="Cutadapt info-file obtained with the argument --info-file;", type=argparse.FileType('r'))
    parser.add_argument('output_file', help="Output file with saved sequences;", type=argparse.FileType('w'))
    parser.add_argument('--reverse_complement', '-r', required=False, help="Reverse-complement sequences that begin with a reverse primer. yes or no;", type=str)
    parser.add_argument('--output_log', '-i', required=False, help="Output log file with id of sequences not saved;", type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--output_counts', '-c', required=False, help="Output log file with counts of sequences by category;", type=argparse.FileType('a'), default=sys.stdout)
    args = parser.parse_args()
    
    file_dict = parse_fastq(args.fastq_file)
    seq_good_primers, seq_wrong_primers = parse_cutadapt_log(args.cutadapt_info_file)
    write_output_files(file_dict, seq_good_primers, seq_wrong_primers, args.output_file, args.reverse_complement, args.output_log, args.output_counts)
