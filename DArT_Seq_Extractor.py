# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 08:42:13 2017

@author: Darren Cullerne
"""

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


#TODO:
#Get ARGSPARSE working
#Parse file parameters to location to save (Default  = Location of input, current = location of Python file)


#import argparse
#
#parser = argparse.ArgumentParser(description='Extract DArT Markers from file and convert to FASTA format')
#
#parser.add_argument('-i','--input', dest='INPUT', metavar='FILE', help='Input DArT file', required=True, type=file)
#parser.add_argument('-o','--output', dest='OUTPUTPREFIX', help='Output file prefix', required=True)
#parser.add_argument('-t','--type', dest='extracttype', help='Which DArT marker types to extract (as single flag): 1=SilicoDArT, 2=SNP1Row, 4=SNP2Row', type=int, choices=range(1,8))
#parser.add_argument('-s','--species', dest='spec', help='Organism which DArT sequences are obtained from \(default = \'species\'\)', required=False, default='species')
#
#args = parser.parse_args()
#
#
#
#parseSilicoDArT = format(args.extracttype,"03b")[0]
#parseSNP1Row = format(args.extracttype,"03b")[1]
#parseSNP2Row = format(args.extracttype,"03b")[2]

parseSilicoDArT = 1
parseSNP1Row = 1
parseSNP2Row = 1


organismPrefix = "CarTin"

if parseSilicoDArT or parseSNP1Row or parseSNP2Row:
    DArT_Import = pd.ExcelFile("/media/u4473753/DC_UNI_ARCHIVE/PhD - Darren/crossing_x/DArT/X395-Report-DSaf16-2055/Report-DSaf16-2055.xlsx")
    #DArT_Import = pd.ExcelFile(args.INPUT)
    #print(DArT_Import.sheet_names)
    
else:
    print "No DArT marker type selected"

#FASTA_output = []

if parseSilicoDArT:
    silico_seq = []    
    SilicoDArT = DArT_Import.parse(sheetname="ScoringData_SilicoDArT", skiprows=6, usecols=['CloneID','AlleleSequence'])
    #f = open("silicodart.fa", 'w')
    
    for x in SilicoDArT.itertuples():
        record=SeqRecord(Seq(str(x[2].rstrip())), id="%s_%s_SilicoDart" % (x[1],organismPrefix), description = "", name = "")
        #output to FASTA
        silico_seq.append(record)

    outputSilico = "%s_%s_silicoDArT.fa" % ("Safflower", organismPrefix)
    SeqIO.write(silico_seq,outputSilico,"fasta")
            
        
    


if parseSNP1Row:
    onerow_seq = []    
    SNP1Row = DArT_Import.parse(sheetname="ScoringData_SNP_1row", skiprows=6, usecols=['CloneID','AlleleSequenceREF']) 
    
    for x in SNP1Row.itertuples():
        name = x[1].split("|")[0]
        SNP = x[1].split("|")[2].split(":")[2].split("-")[0].split(">")[1]
        pos = x[1].split("|")[2].split(":")[1].split("-")[1]
        record=SeqRecord(Seq(str(x[2].rstrip())), id="%s_%s_SNP_%s_%s" % (name,organismPrefix,pos,SNP))
        onerow_seq.append(record)
        
    output1Row = "%s_%s_DArTSNP_1row.fa" % ("Safflower", organismPrefix)
    SeqIO.write(onerow_seq,output1Row,"fasta")


if parseSNP2Row:
    tworow_seq = []
    SNP2Row = DArT_Import.parse(sheetname="ScoringData_SNP_2rows", skiprows=6, usecols=['CloneID','AlleleSequence'])   
    
    for x in SNP2Row.itertuples():
        if x[1].find("--") == -1:
            snptype = "ref"
            nucl = x[1].split("|")[2].split("-")[2].split(":")[1].split(">")[0]
        else:
            snptype = "snp"
            nucl = x[1].split("|")[2].split("-")[2].split(":")[1].split(">")[1]

        name = x[1].split("|")[0]
        pos = x[1].split("|")[2].split("-")[2].split(":")[0]
        
        record=SeqRecord(Seq(str(x[2].rstrip())), id="%s_%s_SNP_%s_%s_%s" % (name,organismPrefix, snptype,pos,nucl))
        tworow_seq.append(record)

    output2Row = "%s_%s_DArTSNP_2row.fa" % ("Safflower", organismPrefix)
    SeqIO.write(tworow_seq,output2Row,"fasta")



#Things to do:

#1) Read though specific tab in the DArT file and find the first marker line (probably a header row)

#2) Parse each line and extract the name and sequence

#3) Rearrange the data and output the file into FASTA format

#4) Dump file output
