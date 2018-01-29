# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:25:58 2018

@author: Darren Cullerne
"""
import argparse
import collections
from Bio.Blast import NCBIXML



parser = argparse.ArgumentParser(description='Extract query and hit description from a BLAST XML file, save as tab separated text file')

parser.add_argument('-i','--input', dest='INPUT', metavar='FILE', help='Input BLASTXML file', required=True, type=argparse.FileType('r'))
parser.add_argument('-o','--output', dest='OUTPUTFILENAME', help='Output file (tab delimited)', required=True)


args = parser.parse_args()


alignment_dict = collections.defaultdict(list)



#blast_records = NCBIXML.parse(open("/media/u4473753/DC_UNI_ARCHIVE1/PhD - Darren/Andrew/txOnCTg.nearDartMarkersx70_1kbp-updown_ctS317g_c0007401_6.xml"),'r')
blast_records = NCBIXML.parse(args.INPUT)
blast_records = list(blast_records)

for blast_record in blast_records:
    # for header in blast_records:
    # print header.reference,
    temp_query = blast_record.query
    #print(blast_record.query)
    if (blast_record.alignments):
        for alignment in blast_record.alignments:
        #temp_hit_def = alignment.hit_def
        #alignment_array.append({temp_query,alignment.hit_def})
            #print(alignment.hit_def)
            alignment_dict[temp_query].append(str(alignment.hit_def))
        
    else:    
        #print("No Hits")
        alignment_dict[temp_query].append("No Hits")
        #NB: In here, you will need to work out how to make a key and index or sanme key, different value sort of thing


#d[k].append(v)


#print(sorted(alignment_dict.items()))

outputfile="./%s" % args.OUTPUTFILENAME

with open(outputfile, "w") as output:
    for k,v in alignment_dict.iteritems():
        #print "%s\t%s" % (str(k), str(v)[2:-2])
        #writer.writerows(sorted(alignment_dict.items()))
        output.write("%s\t%s\n" % (str(k), str(v)[2:-2]))
        
