#!/usr/bin/env python
'''Reads in nucleotide sequences from a fasta file as SeqRecords and sorts into a dictionary,
with unique sequences as keys and lists of SeqRecords sharing these sequevars as values. If the user
specified a fasta file containing sequences which have been empirically determined to amplify in the 
PCR (i.e. "safelist"), the sequevars are compared to this safelist prior to being output to fasta. Sequence
identifiers (usually Genbank accession numbers) of sequences in a sequevar are printed to a text file.

Usage with safelist: python sequevarsToDict.py -s safelist.fasta unmatching_sequences.fasta  Mar20_sequevars
Usage without safelist: python sequevars_to_dict.py unmatching_sequences.fasta  Mar20_sequevars
'''

'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory, March 2018'''

import sys,string,os, time, Bio, re, argparse
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData

#parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-s", dest = "safelist", help = "specify filename of previously vetted sequences",
    action = "store")
parser.add_argument("fastaToParse")
parser.add_argument("outFileHandle")
results = parser.parse_args()
print("ARGUMENTS:")
print(results)
print('safelist = ', results.safelist)

#construct output filenames
gb_accessionsHandle = results.outFileHandle + ".txt" #user-specified output file name
gb_accessions= open(gb_accessionsHandle,'w') #output file of GB access #'s and details re: which oligos NOT found
outputFastaHandle = results.outFileHandle + ".fasta" #create fasta file from  user specified name

def outputInformationFile(sequevars):
    '''Outputs text file with unique amplicon sequences and id's of sequences in each sequevar.'''
    localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
    gb_accessions.write("---------------------------------------------------------------------------\n")
    gb_accessions.write("RESULTS OF SEARCH FOR EXACT MATCHES TO OLIGONUCLEOTIDES IN QUERY SEQUENCES:\n")
    gb_accessions.write("---------------------------------------------------------------------------\n\n")
    gb_accessions.write("Analysis as of: " + localtime + "\n\n")
    gb_accessions.write("\n---------------------------------------------------------------------------\n")
    #if there are sequevars for further analysis, output them to the results .txt and fasta files
    if len(sequevars) > 0: #if there are sequevars to output, write to list with sequence id's in sequevar
        outputSequevarsToFasta(sequevars) #output fasta file of sequevars i.e.representative record
        gb_accessions.write("%i Different Amplicon Sequevars Found: %s \n\n" % (len(sequevars),results.fastaToParse))
        for sequevar in sequevars:
            print("\nSequevar: %s" % sequevar) #print sequevar to console
            recordList = sequevars[sequevar]
            gb_accessions.write("\n\n%i Records in Sequevar: %s" % (len(recordList),sequevar)) #write sequevar to file
            print("%i Record(s) have this sequence:" % (len(recordList)))
            for record in recordList:
                print("\t" + record.id) #print record id to console
                gb_accessions.write("\n\t" + record.id) #write records with this sequence to file
    else: #if not sequevars to output, print message to outfile to clarify this
        gb_accessions.write("No sequevars requiring further investigation.")
    gb_accessions.write("\n---------------------------------------------------------------------------")
    gb_accessions.write("\n\nEND OF RESULTS\n")
    #print output filepaths to console and to search results file
    gb_accessions.write("\nSearch Results (this file): " + gb_accessionsHandle)
    print("\nSearch Results: " + gb_accessionsHandle)
    if len(sequevars) > 0: #if there are sequevars to print to fasta, direct user to fasta file
        gb_accessions.write("\nFasta sequences for analysis: " + outputFastaHandle)
        print("Fasta sequences for analysis: " + outputFastaHandle)
    else:
        gb_accessions.write("\nFasta sequences for analysis: N/A")
        print("Fasta sequences for analysis: N/A\n")
    return

def outputSequevarsToFasta(sequevars):
    '''Output a fasta file with one representative sequence for each sequevar.'''
    output_fasta = open(outputFastaHandle, 'w') #create a writeable fasta file
    sequevars_list = []
    for sequevar in sequevars:
        first_record = sequevars[sequevar][0] #grab the first SeqRecord in the list sharing the sequevar
        sequevars_list.append(first_record)  #add it to the list of SeqRecords
    SeqIO.write(sequevars_list, output_fasta, "fasta") #write the sequevars to fasta
    return

with open(results.fastaToParse,'r') as inFile:
    #read nucleotide sequences from fasta file into SeqRecords, uppercases and adds to a list
    seqList = [rec.upper() for rec in list(SeqIO.parse(inFile, "fasta", alphabet=IUPAC.ambiguous_dna))]
    sequevars = {} # empty dictionary of sequevar: list<SeqRecord>
    safe_seqs = [] # empty list of nucleotide sequences

    if results.safelist: #if safelist specified by user, read these sequences into a list
        safe_seqs = [record.seq for record in SeqIO.parse(results.safelist, "fasta", alphabet=IUPAC.ambiguous_dna)]
        print("Safeseqs:")
        for seq in safe_seqs:
            print(seq) #print all the vetted unique sequences
        
    for record in seqList: #parse SeqRecord for exact match to amplicon or reverse complement
        sequence = str(record.seq)
        #check if the sequence has already been tested and vetted (i.e. in safelist)
        if (sequence in safe_seqs): #only add to dictionary if not previously tested
            print("Sequence from %s in safe list!" % (record.id) )#print mssg to console
        else:
            #if the sequence is a key in the dict, add SeqRecord to list
            if sequence in sequevars:
                sequevars[sequence].append(record)
            else:
                #add sequence as new key to dict, accessing a list with SeqRecord  as its first item
                sequevars[sequence] = [record]
    #get a sorted list of unique sequence keys
    sorted_unique_seq_keys = sorted(sequevars.keys())
    #process list of SeqRecords in each sequevar and print to results summary file for end-user        
    outputInformationFile(sequevars) #output information .txt file
    #outputSequevarsToFasta(sequevars) #output fasta file of sequevars i.e.representative record      

inFile.close()
gb_accessions.close()
