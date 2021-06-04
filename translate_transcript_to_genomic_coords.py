#!/usr/bin/env python3
"""
A script that reads in 2 required files as command line arguments and translates the input transcript coordinates to
genomic coordinates using a position and CIGAR string. The (0-based) position is input as CHR1:3,
and the CIGAR string is 8M7D6M2I2M11D7M.

Input args:
--input_file1: file path
A four column (tab-separated) file containing the transcripts. The first column is the transcript
name, and the remaining three columns indicate it’s genomic mapping: chromosome name,
0-based starting position on the chromosome, and CIGAR string indicating the mapping.

--input_file2: file path
A two column (tab-separated) file indicating a set of queries. The first column is a transcript
name, and the second column is a 0-based transcript coordinate.

@Author: Ankit Jambusaria
@Date: 5/27/2021
"""

import os
import sys
import re
from argparse import ArgumentParser


# Function for reading in the genomic mapping file (input_file1).
# Create a dictionary which stores the chr, start coordinate, and cigar string for each transcript-chr mapping
def create_transcript_genomic_dict(genome_mapping_file):
    transcript_to_genomic_dict = {}
    genome_map = open(genome_mapping_file, "r")
    for row in genome_map:
        row = row.rstrip()
        rows = re.split(r'\s+', row)

        # check input genome_map format
        check_genome_map_line_format(rows)

        transcript_id = rows[0]
        chr_num = rows[1]
        start_coord = rows[2]
        cigar = rows[3]

        transcript_to_genomic_dict[transcript_id] = {}
        transcript_to_genomic_dict[transcript_id][chr_num] = {}
        transcript_to_genomic_dict[transcript_id][chr_num]["start_coord"] = int(start_coord)
        transcript_to_genomic_dict[transcript_id][chr_num]["cigar"] = process_cigar_string(cigar)

    genome_map.close()
    return transcript_to_genomic_dict


# Function that verifies the input genome mapping file format
# verifies that every line in the genome mapping file contains 4 fields
# verifies that 3rd column in the genome mapping file is an integer
# verifies that the cigar string contains valid characters
def check_genome_map_line_format(genome_map_line_arr):
    # checks number of columns
    if not len(genome_map_line_arr) == 4:
        print("Error! Genome mapping file must have 4 columns")
        sys.exit()

    # checks if the third column contains an integer
    if not genome_map_line_arr[2].isdigit():
        print("Error! Reference coordinate is invalid")
        sys.exit()

    if not (re.search(r'\d+[MIDNSHP=X]', genome_map_line_arr[3], re.I)):
        print("Error! CIGAR string contains invalid characters")
        sys.exit()


# Function that processes cigar string and returns a list where each element of list
# contains an integer and cigar char
def process_cigar_string(cigar_string):
    cigar_list = list()
    cigar_int = ''
    for cigar_char in cigar_string:
        if not (cigar_char.isalpha() or cigar_char == '='):
            cigar_int = cigar_int + cigar_char
        else:
            if cigar_int == '':
                print("Error! Wrong CIGAR string format")
                sys.exit()
            cigar_int = int(cigar_int)
            cigar_list.append([cigar_int, cigar_char])
            cigar_int = ''
    return cigar_list


# Function to establish coordinate mappings using the cigar list
# internally calls the "process_cigar_list" function that returns a dictionary containing coordinate correspondence
def map_coordinates(transcript_to_genomic_dict):
    coord_map = {}
    for tr in transcript_to_genomic_dict:
        coord_map[tr] = {}
        for ch in transcript_to_genomic_dict[tr]:
            coord_map[tr][ch] = {}
            start_coord = transcript_to_genomic_dict[tr][ch]["start_coord"]
            cigar_vec = transcript_to_genomic_dict[tr][ch]["cigar"]
            genomic_dict = generate_genomic_dict(start_coord, cigar_vec)
            coord_map[tr][ch] = genomic_dict
    return coord_map


# Function that returns a dictionary containing the genomic coordinate corresponding
# to the transcript coord and chromosome
# Assumptions:
# 1. Every transcript coordinate is mapped to a unique genomic coordinate
# 2. Cigar string contains only 3 chars (M,D,I)
def generate_genomic_dict(chr_start_coord, cigar_arr):
    genomic_dict = {}
    tr_idx = 0
    chr_idx = chr_start_coord

    for entry in cigar_arr:
        cigar_int = entry[0]
        cigar_char = entry[1]

        # input: query and reference
        if re.match(r'[M=X]', cigar_char, re.I):
            for idx in range(cigar_int):
                genomic_dict[tr_idx] = chr_idx
                tr_idx += 1
                chr_idx += 1
        # input: reference
        elif re.match(r'[DN]', cigar_char, re.I):
            for idx in range(cigar_int):
                chr_idx += 1
        # input: query
        elif re.match(r'[IS]', cigar_char, re.I):
            for idx in range(cigar_int):
                tr_idx += 1

    return genomic_dict


# Function to read in transcript processing file and coord_map from the map_coordinates function
# results are written to --output
def merge_transcript_file(transcript_processing_filename, coord_map, output):
    transcript_processing_file = open(transcript_processing_filename, "r")
    output_file = open(output, "w")
    for query in transcript_processing_file:
        query = query.rstrip()
        queries = re.split(r'\s+', query)

        # check format
        check_transcript_line_format(queries)

        tr_id = queries[0]
        tr_coord = queries[1]

        # checks if the transcript id is known
        if not (tr_id in coord_map):
            print("Transcript", tr_id, "does not exist in genome mapping file")
            continue

        # get the corresponding genomic coordinate for the given transcript coordinate
        for ch in coord_map[tr_id]:
            # checks if transcript coordinate is defined
            if not (int(tr_coord) in coord_map[tr_id][ch]):
                print("Transcript coordinate", tr_coord, "for transcript", tr_id, "does not exist")
                continue

            output_file.write(tr_id + "\t" + tr_coord + "\t" + ch + "\t" + str(
                coord_map[tr_id][ch][int(tr_coord)]) + "\n")

    transcript_processing_file.close()
    output_file.close()


# Function to verify transcript processing file format
def check_transcript_line_format(transcript_line_arr):
    # check number of columns
    if not len(transcript_line_arr) == 2:
        print("Error! Transcript Processing File must contain 2 columns")
        sys.exit()

    # check second column for integers only
    if not transcript_line_arr[1].isdigit():
        print("Error! Transcript coordinates must be integers")
        sys.exit()


# Function to validate the input files and return custom error message for invalid inputs
def validate_input_args(input_args):
    if not os.path.isfile(input_args.genome_mapping_file):
        return False, "Genome mapping file does not exist"
    if not os.path.isfile(input_args.transcript_processing_file):
        return False, "Transcript processing file does not exist"
    if not os.path.isdir(os.path.dirname(os.path.abspath(input_args.output_file))):
        return False, "Output file location does not exist"

    return True, ""


# Function that executes the entire workflow
def transcript_to_genomic_coordinates(genome_mapping_file, transcript_processing_file, output):
    transcript_genomic_alignment = create_transcript_genomic_dict(genome_mapping_file)
    genomic_coords = map_coordinates(transcript_genomic_alignment)
    merge_transcript_file(transcript_processing_file, genomic_coords, output)


######################################################################
######################################################################
#######  MAIN  #######################################################
######################################################################
######################################################################

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--genome-mapping-file", required=True, dest="genome_mapping_file",
                        help="file containing the transcripts (e.g., input_file1.txt)")
    parser.add_argument("--transcript-processing-file", required=True, dest="transcript_processing_file",
                        help="file containing a set of queries (e.g., input_file2.txt")
    parser.add_argument("--output", required=False, dest="output_file", default='output.txt',
                        help="Filename for output file. Default: output.txt)")
    args = parser.parse_args()

    (is_input_valid, msg) = validate_input_args(args)
    if not is_input_valid:
        sys.stderr.write(msg + "\n")
        sys.exit(-1)

transcript_to_genomic_coordinates(args.genome_mapping_file,
                                  args.transcript_processing_file,
                                  args.output_file)
