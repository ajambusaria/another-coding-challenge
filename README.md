The program can be executed using the following command:
python3 translate_transcript_to_genomic_coords.py --genome-mapping-file input_file1.txt --transcript-processing-file input_file2.txt --output output.txt

translate_transcript_to_genomic_coords.py is a script that reads in a genome mapping file and a transcript processing file which are given as command line arguments and translates input transcript coordinates to genomic coordinates using the position and CIGAR string. The (0-based) position is input as CHR1:3, and the CIGAR string is 8M7D6M2I2M11D7M.

Assumptions: 
•	The transcript is always mapped from genomic 5’ to 3’.
•	The transcript and genomic coordinates are 0-based.
•	The transcript processing file contains transcript ids where each transcript maps to a unique location on at most one chromosome. 

