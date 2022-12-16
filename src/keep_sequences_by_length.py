#!/usr/bin/env python3


import os, argparse, sys

"""
USAGE
    ./keep_sequences_by_length.py [-h] <fastq_file> <min_length> <max_length> <output_file> [-i <output_log>] [-c <output_counts>]
DESCRIPTION
    keep_sequences_by_length.py is a script to keep sequences that have an expected length.
PREREQUISITE
    - python3
INFO
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
    Detects reads that are too long or too short, keep seqId into lists and create a dict without them
    """
    idSeqTooShort = list()
    idSeqTooLong = list()
    newfastqDict = dict()
    for seqId in fastqDict.keys():
        if len(fastqDict[seqId]["seq"]) < min_length:
            idSeqTooShort.append(seqId)
        elif len(fastqDict[seqId]["seq"]) > max_length:
            idSeqTooLong.append(seqId)
        else:
            newfastqDict[seqId] = fastqDict[seqId]
    return newfastqDict, idSeqTooShort, idSeqTooLong 


def write_output_files(newfastqDict, fastqDict, idSeqTooShort, idSeqTooLong, output_file, output_log, output_counts):
    """
    Write the new FASTQ file without the reads discriminated and the log file with all id of reads not keeped
    """
    # File with saved sequences
    for seqId in newfastqDict.keys():
        print("@"+seqId, file=output_file)
        print(newfastqDict[seqId]["seq"], file=output_file)
        print("+", file=output_file)
        print(newfastqDict[seqId]["qual"], file=output_file)
    # Log files
    # id
    print("\n*---------- Sequences too short removed ("+str(len(idSeqTooShort))+") : ", file=output_log)
    output_log.write("\n".join(idSeqTooShort))
    print("\n*---------- Sequences too long removed ("+str(len(idSeqTooLong))+") : ", file=output_log)
    output_log.write("\n".join(idSeqTooLong))
    # counts
    outputName = os.path.split(output_file.name)[1] 
    outputName = os.path.splitext(outputName)[0]
    print(outputName+"\t"+str(len(newfastqDict))+"\t"+str(len(idSeqTooShort))+"\t"+str(len(idSeqTooLong)), file=output_counts)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="keep_sequences_by_length.py is a script to keep sequences that have an expected length.")
    parser.add_argument('fastq_file', help="FASTQ file, not compressed, with sequences;", type=argparse.FileType('r'))
    parser.add_argument('min_length', help="Minimum length of a sequence to be keep;", type=int)
    parser.add_argument('max_length', help="Maximum length of a sequence to be keep;", type=int)
    parser.add_argument('output_file', help="Output file with saved sequences;", type=argparse.FileType('w'))
    parser.add_argument('--output_log', '-i', required=False, help="Output log file with id of sequences not saved;", type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--output_counts', '-c', required=False, help="Output log file with counts of sequences by category;", type=argparse.FileType('a'), default=sys.stdout)
    args = parser.parse_args()
    
    file_dict = parse_fastq(args.fastq_file)
    new_file_dict, id_seq_too_short, id_seq_too_long = discriminate_by_length(file_dict, args.min_length, args.max_length)
    write_output_files(new_file_dict, file_dict, id_seq_too_short, id_seq_too_long, args.output_file, args.output_log, args.output_counts)