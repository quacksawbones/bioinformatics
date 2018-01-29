# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:25:58 2018

@author: Darren Cullerne
"""
import argparse
import collections
from Bio.Blast import NCBIXML


#Only need input and output. Can expand later for number of alignments etc...
parser = argparse.ArgumentParser(description='Extract query and hit description from a BLAST XML file, save as tab separated text file')

parser.add_argument('-i','--input', dest='INPUT', metavar='FILE', help='Input BLASTXML file', required=True, type=argparse.FileType('r'))
parser.add_argument('-o','--output', dest='OUTPUTFILENAME', help='Output file (tab delimited)', required=True)

args = parser.parse_args()

#Create new Dictionaly type
alignment_dict = collections.defaultdict(list)

#Parse XML file
blast_records = NCBIXML.parse(args.INPUT)
blast_records = list(blast_records)

#Iterate through, if no hits, return "No Hits", otherwise, return hit descriptions
for blast_record in blast_records:
    temp_query = blast_record.query
    if (blast_record.alignments):
        for alignment in blast_record.alignments:
            alignment_dict[temp_query].append(str(alignment.hit_def))        
    else:    
        alignment_dict[temp_query].append("No Hits")

#Create filename string
outputfile="./%s" % args.OUTPUTFILENAME

#Open file, output key and each value as a list of "pairs"
with open(outputfile, "w") as output:
    for x in alignment_dict.keys():
        for y in alignment_dict[x]:
            #print "Key=%s, value=%s" % (x,y)
            output.write("%s\t%s\n" % (x,y))
            
