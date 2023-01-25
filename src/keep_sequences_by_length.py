#!/usr/bin/env python3


import os, argparse, sys

"""
USAGE
    ./keep_sequences_by_length.py [-h] <fastq_file> <min_length> <max_length> <output_file> <output_file2> [-i <output_log>] [-c <output_counts>]
DESCRIPTION
    keep_sequences_by_length.py is a script to keep sequences that have an expected length.
PREREQUISITE
    - python3
INFO
    This script is a part of QuikAmpliclean
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
            seqId = line.strip().replace("@","").split()[0]
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


def discriminate_by_length(fastqDict, min_length, max_length):
    """
    Detects reads that are too long or too short
    """
    seqWrongLength = {"too_long":[], "too_short":[]}
    seqGoodLength = []
    for seqId in fastqDict.keys():
        if len(fastqDict[seqId]["seq"]) < min_length:
            seqWrongLength["too_short"].append(seqId)
        elif len(fastqDict[seqId]["seq"]) > max_length:
            seqWrongLength["too_long"].append(seqId)
        else:
            seqGoodLength.append(seqId)
    return seqGoodLength, seqWrongLength 


def write_output_files(seqGoodLength, fastqDict, seqWrongLength, output_file, output_file2, output_log, output_counts):
    """
    Write the new FASTQ file without the sequences discriminated, a FASTQ file with sequences unkeep and the log file with all id of reads not keeped
    """
    # File with saved sequences
    for seqId in seqGoodLength:
        print("@"+seqId, file=output_file)
        print(fastqDict[seqId]["seq"], file=output_file)
        print("+", file=output_file)
        print(fastqDict[seqId]["qual"], file=output_file)
    # File with unsaved sequences
    for typeErrorLength in seqWrongLength.keys():
        for seqId in seqWrongLength[typeErrorLength]:
            print("@"+seqId, file=output_file2)
            print(fastqDict[seqId]["seq"], file=output_file2)
            print("+", file=output_file2)
            print(fastqDict[seqId]["qual"], file=output_file2)
    # Log files
    # id
    print("\n*---------- Sequences too short removed ("+str(len(seqWrongLength["too_short"]))+") : ", file=output_log)
    output_log.write("\n".join(seqWrongLength["too_short"]))
    print("\n*---------- Sequences too long removed ("+str(len(seqWrongLength["too_long"]))+") : ", file=output_log)
    output_log.write("\n".join(seqWrongLength["too_long"]))
    # counts
    outputName = os.path.split(output_file.name)[1] 
    outputName = os.path.splitext(outputName)[0]
    print(outputName+"\t"+str(len(seqGoodLength))+"\t"+str(len(seqWrongLength["too_short"]))+"\t"+str(len(seqWrongLength["too_long"])), file=output_counts)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="keep_sequences_by_length.py is a script to keep sequences that have an expected length.")
    parser.add_argument('fastq_file', help="FASTQ file, not compressed, with sequences;", type=argparse.FileType('r'))
    parser.add_argument('min_length', help="Minimum length of a sequence to be keep;", type=int)
    parser.add_argument('max_length', help="Maximum length of a sequence to be keep;", type=int)
    parser.add_argument('output_file', help="Output file with saved sequences;", type=argparse.FileType('w'))
    parser.add_argument('output_file2', help="Output file with unsaved sequences;", type=argparse.FileType('w'))
    parser.add_argument('--output_log', '-i', required=False, help="Output log file with id of sequences not saved;", type=argparse.FileType('a'), default=sys.stdout)
    parser.add_argument('--output_counts', '-c', required=False, help="Output log file with counts of sequences by category;", type=argparse.FileType('a'), default=sys.stdout)
    args = parser.parse_args()
    
    file_dict = parse_fastq(args.fastq_file)
    seq_good_length, seq_wrong_length = discriminate_by_length(file_dict, args.min_length, args.max_length)
    write_output_files(seq_good_length, file_dict, seq_wrong_length, args.output_file, args.output_file2, args.output_log, args.output_counts)