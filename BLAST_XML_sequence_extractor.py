import argparse
import operator
import math
from Bio.Blast import NCBIXML
from Bio import SeqIO

alignment_array = []
region = []
temp_region = []
candidate = []


parser = argparse.ArgumentParser(description='Parse BLAST XML files for finding stacked alignments on contigs/scaffolds (nucleotide only) v0.1')
parser.add_argument('-i','--input', dest='BLASTXMLINPUT', metavar='FILE', help='Input file as Blast XML', required=True, type=file)
parser.add_argument('-I','--dbinput', dest='BLASTDBFASTA', metavar='FILE', help='Input blast db FASTA file', required=True, type=file)
parser.add_argument('-o','--output', dest='CANDIDATEOUTPUTDIR', help='Output directory (default is current directory)', required=False, default='./')
parser.add_argument('-t','--type', dest='gene_type', help='Gene family BLASTed (default is \'gene\'', required=False, default='gene')
parser.add_argument('-s','--spec', dest='spec', help='species BLASTed (default is \'gene\'', required=False, default='gene')
parser.add_argument('-l','--length', dest='LENGTH', help='Length of alignment (default = 500bp)', required=False, default='500', type=int)
parser.add_argument('-g','--gap', dest='GAP', help='Allowed length of gap (default = 250bp)', required=False, default='250', type=int)
args = parser.parse_args()

blast_records = NCBIXML.parse(args.BLASTXMLINPUT)
blast_records = list(blast_records)

fasta_records = list(SeqIO.parse(args.BLASTDBFASTA, "fasta"))

#blast_record = blast_records.next()
for blast_record in blast_records:
	# for header in blast_records:
		# print header.reference,
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			alignment_array.append({"Hit":blast_record.query, "Start":hsp.query_start, "End":hsp.query_end})

alignment_array.sort(key=operator.itemgetter("Hit","Start","End"))



for x in range(0,(len(alignment_array))):
	if x == 0:
		temp_region = [alignment_array[0]["Hit"], alignment_array[0]["Start"], alignment_array[0]["End"]]
		continue
		
	elif alignment_array[x]["Hit"] == alignment_array[x-1]["Hit"]:
		if alignment_array[x]["Start"] == temp_region[1]:
			if alignment_array[x]["End"] > temp_region[2]:
				temp_region[2] = alignment_array[x]["End"]
				continue
		elif (alignment_array[x]["Start"] - temp_region[2] < args.LENGTH):	
			if (alignment_array[x]["End"] > temp_region[2]):
				temp_region[2] = alignment_array[x]["End"]
			
		elif ((alignment_array[x]["Start"] - temp_region[1]) > args.LENGTH):
			region.append(temp_region)
			temp_region = [alignment_array[x]["Hit"], alignment_array[x]["Start"], alignment_array[x]["End"]]
		
	elif alignment_array[x]["Hit"] != alignment_array[x-1]:
		region.append(temp_region)
		temp_region = [alignment_array[x]["Hit"], alignment_array[x]["Start"], alignment_array[x]["End"]]		
		

		
for y in range(0,len(region)):
	try:
		if (region[y][2] - region[y][1]) > args.GAP:
			candidate.append(region[y])
	except IndexError:
		continue


		
for z in range(0,len(candidate)):

	for f in range(0,len(fasta_records)):
		
		if candidate[z][0] == fasta_records[f].id:
			if candidate[z][1] - 200 <= 0:
				candidate[z][1] = 0
			else:
				candidate[z][1] = candidate[z][1] - 200

			if candidate[z][2] + 200 > len(fasta_records[f].seq):
				candidate[z][2] = len(fasta_records[f].seq)
			else:
				candidate[z][2] = candidate[z][2] + 200

			subseq = "".join(fasta_records[f].seq[(candidate[z][1]):(candidate[z][2])])

			candidate[z].append(subseq)


			# output as FASTA format	

for r in range(0,len(candidate)):

	# filename = (args.spec+'_'+args.gene_type+'_'+candidate[z][0]+'_'+str(candidate[z][1] - 200)+'_'+str(candidate[z][2] + 200)+'_+-200bp_nt.fasta')
	# output = open(args.CANDIDATEOUTPUTDIR+"/"+filename, 'w')
	# output.write('>'+args.spec+'_'+args.gene_type+'_'+candidate[z][0]+'_'+str(candidate[z][1] - 200)+'_'+str(candidate[z][2] + 200)+'_+-200bp_nt\n')
	# output.write(candidate[z][3])
	# output.close()	
		
	filename = (args.spec+'_'+args.gene_type+'_'+candidate[r][0]+'_'+str(candidate[r][1])+'_'+str(candidate[r][2])+'_200bp_nt.fasta')
	output = open(args.CANDIDATEOUTPUTDIR+"/"+filename, 'w')
	output.write('>'+args.spec+'_'+args.gene_type+'_'+candidate[r][0]+'_'+str(candidate[r][1])+'_'+str(candidate[r][2])+'_200bp_nt\n')
	output.write(candidate[r][3])
	output.close()	
		
		
		
