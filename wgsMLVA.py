#WHOLE GENOME SEQUENCE MULTIPLE LOCUS VARIABLE NUMBER TANDEM REPEAT SEARCH 
#last edit: 03/13/2018

#Imported Modules.
import sys
import re
import os
import argparse
from Bio import SeqIO

#------------------------------------------------------------------------------
#Program requires 2 arguments to be entered and will autormatically
#load in the VNTR chart data.

parser = argparse.ArgumentParser(prog = " MLVA in silico Program", \
								usage = 'Find MLVA code from closed \
										fasta sequence.')
parser.add_argument('-in', '--inputFile', required = True , \
					help = 'input fasta file. Required.')
parser.add_argument('-out', '--outputDir', required = True, \
					help = 'MLVA output directory. Required.')
parser.add_argument('-v', '--version', action = 'version', \
					version = '%(prog)s version 1.0')
parser.add_argument('-n', '--name', required = False, \
					help = 'Assign file name. If not assigned it will use \
							the file name of the fasta file.')
parser.add_argument('-u', '--unmask', action = "store_true", default= False, \
					help = 'Removes VNTR3 masking, gives actual VNTR3 \
							repeat values.')
parser.add_argument('-k', '--keepTemp', action = "store_true", default= False, \
					help = 'Keeps single-line temporary fasta file.')
args = parser.parse_args()
inFile = args.inputFile
outDir = args.outputDir
mlvaDir = os.path.dirname(os.path.realpath(sys.argv[0]))
sampleID = args.name
fix = args.unmask
keep = args.keepTemp

#------------------------------------------------------------------------------
#Check parsed arguments

color_red = '\33[31m'

if inFile.lower().endswith(".fasta"):
	pass
else:
	color_red = '\33[31m'
	print(color_red + "\n>>> Error: Please input fasta\n" + color_red)
	sys.exit()
#------------------------------------------------------------------------------
#Stores variables for program.
#Variable "name" changes the name of the files and directories created.
if sampleID:
	name = sampleID
	outFile = name + "fasta_out"
	outSummary = name + "summary"
else:	
	name = re.findall("(\w\d\d\d)",inFile)[0]
	outFile = name + "fasta_out"
	outSummary = name + "summary"
#------------------------------------------------------------------------------
#Create "MLVA" directory to hold results, creates directory in MLVA to hold
#each isolate.
if not os.path.exists("%s/MLVA" % outDir):
	os.mkdir("%s/MLVA" % outDir)
if not os.path.exists("%s/MLVA/summary.txt" % outDir):
	with open("%s/MLVA/summary.txt" % outDir, "a") as summary:
			summary.write(\
			"Isolate: MLVA_number [VNTR1, VNTR3a, VNTR3b, VNTR4, VNTR5, VNTR6]\n")	
if not os.path.exists("%s/MLVA/MLVA_data" % outDir):
	os.mkdir("%s/MLVA/MLVA_data" % outDir)		
if not os.path.exists("%s/MLVA/sequence_data" % outDir):
	os.mkdir("%s/MLVA/sequence_data" % outDir)
#------------------------------------------------------------------------------
print "\n>>> Opening", name, "fasta file..."

#Reads sequence file list and stores it as a string object. Close file.
#Splits string at the start of a line.
#First fasta in the file is split into an empty element and the first fasta.
#"del" removes this empty element.
with open(inFile,"r") as newFile:
	sequences = newFile.read()
	sequences = re.split("^>", sequences, flags = re.MULTILINE)
	del sequences[0]
	newFile.close()
#Conversts multiline fasta to single line. Writes new fasta to file: "X_out".
#1. Split each fasta into header and sequence.
#2. Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
#3. Replace newlines in sequence, remembering to add one to the end.
with open("%s/MLVA/MLVA_data/%s" % (outDir, outFile),"w") as newFasta:
	for fasta in sequences:
		try:
			header, sequence = fasta.split("\n", 1)
		except ValueError:
			print fasta
		header = ">" + header + "\n"
		sequence = sequence.replace("\n","") + "\n"
		newFasta.write(sequence)
	newFasta.close()
#------------------------------------------------------------------------------
#Each primer and reverse primer for genome search. 
#V(1 to 6)p(1 or 2)(start or end) = VNTR(1 to 6)primer(1 or 2)(start or end)
V1p1s = "AAAATTGCGGCATGTGGGCTGACTCTGA"
V1p1e = "CACCACCACGTCTCCCGCCGCCAGG"
V1p2s = "CCTGGCGGCGGGAGACGTGGTGGTG"
V1p2e = "TCAGAGTCAGCCCACATGCCGCAATTTT"
V3p1s = "GCCTCGGCGAAATTGCTGAAC"
V3p1e = "GGTCTCGGGCGTTTCCTCGCCCGC"
V3p2s = "GCGGGCGAGGAAACGCCCGAGACC"
V3p2e = "GTTCAGCAATTTCGCCGAGGC"
V4p1s = "GCCGCTGCTCGACGCCAGGGACAA"
V4p1e = "CAGGTCCAGGCGCAGGGCACG"
V4p2s = "CGTGCCCTGCGCCTGGACCTG"
V4p2e = "TTGTCCCTGGCGTCGAGCAGCGGC"
V5p1s = "TGCCGGGTTTCGGCATCTCGATGGGATACG"
V5p1e = "AAGAGCCTGGAGCTCGGGTGGGCCGGCTTC"
V5p2s = "GAAGCCGGCCCACCCGAGCTCCAGGCTCTT"
V5p2e = "CGTATCCCATCGAGATGCCGAAACCCGGCA"
V6p1s = "CCAACGGCGGTCTGCTGGGTGGTC"
V6p1e = "GGTAGCGGCGCAGCGGGCGGCG"
V6p2s = "CGCCGCCCGCTGCGCCGCTACC"
V6p2e = "GACCACCCAGCAGACCGCCGTTGG"

#VNTR repeats
VNTR1wt1 = "CTGCTTGGCGGGTTC"
VNTR1wt2 = "GAACCCGCCAAGCAG"
VNTR1mut1_13 = "CTGCTTGGCGGGCTC"
VNTR1mut2_13 = "GAGCCCGCCAAGCAG"
VNTR1mut1_1315 = "CTGCTTGGCGGGCTG"
VNTR1mut2_1315 = "CAGCCCGCCAAGCAG"
VNTR3_1 = "CTGGC"
VNTR3_2 = "GCCAG"
VNTR4_1 = "CAAGGACAAGGG"
#VNTR4mut1_1 = "TAAGGACAAGGG"
VNTR4_2 = "CCCTTGTCCTTG"
VNTR5_1 = "TGGTGC"
VNTR5_2 = "GCACCA"
VNTR6_1 = "GGCGGCTCG"
VNTR6_2 = "CGAGCCGCC"

#Order: [VNTR1, VNTR3_1, VNTR3b_2, VNTR4, VNTR5, VNTR6]
MLVA_code = []

#Creating and opening output file.
MLVAfile = name + "MLVA.txt"
q = open("%s/MLVA/MLVA_data/%s" % (outDir, MLVAfile), "w")
seq = open("%s/MLVA/sequence_data/%s" % (outDir, name), "a")

print >>q, name
with open("%s/MLVA/MLVA_data/%s" % (outDir, outFile), 'r') as fasta:
	all_dna = fasta.read()

print ">>> Starting MLVA Search..."

#------------------------------------------------------------------------------
#attempt at doing a fuzzy repeat search. Still WIP.
def fuzzy_search(search_key, text, strictness):
	lines  = text
	for line in lines:
		similarity = SequenceMatcher(None, line, search_key)
		if similarity.ratio() > strictness:
			print "yes"

#------------------------------------------------------------------------------
#all VNTR search scripts. 
def VNTR1_search(head, p1s, p1e, p2s, p2e, r1, r2, r3, r4, r5, r6):

	print>>q, "\n>>> " + head

	if p1s and p1e in all_dna:
		for x in re.finditer(p1s, all_dna):
			fp1 = int(x.start())
		for y in re.finditer(p1e, all_dna):
			fp2 = int(y.end())
		fDistance = fp2 - fp1
		fragment = all_dna[fp1:fp2]
		if r1 in fragment:
			#print fuzzy_search(r1, fragment, 0.85)
			#vhtr_count = fuzzy_search(r1, fragment, 0.85)
			#vntr_count = difflib.get_close_matches(r1, fragment, cutoff = 0.85)
			#print vntr_count
			vntr_count = fragment.count(r1) + fragment.count(r3) + fragment.count(r5)
			MLVA_code.append(vntr_count)	
		else:
			#print fuzzy_search(r2, fragment, 0.85)
			#vhtr_count = fuzzy_search(r2, fragment, 0.85)
			#vntr_count = difflib.get_close_matches(r2, fragment, cutoff = 0.85)
			vntr_count = fragment.count(r2) + fragment.count(r4) + fragment.count(r6)
			MLVA_code.append(vntr_count)
		print >>seq, ">%s" % (name+"_" +head), "\n" + fragment
		print >>q, "Fragment locus:", fp1,"-",fp2
		print >>q, "Fragment size:", fDistance	
		print >>q, "Number of repeats:", vntr_count
	else:
		pass

	if p2s and p2e in all_dna:
		for x in re.finditer(p2s, all_dna):
			rp1 = int(x.start())
		for y in re.finditer(p2e, all_dna):
			rp2 = int(y.end())
		rDistance = rp2 - rp1
		fragment = all_dna[rp1:rp2]
		replace = {"A": "T", "T": "A", "G": "C", "C": "G"}
		replace = dict((re.escape(k), v) for k, v in replace.iteritems())
		pattern = re.compile("|".join(replace.keys()))
		rev_text = pattern.sub(lambda m: replace[re.escape(m.group(0))], fragment)
		fragment_rev = rev_text[::-1]
		if r1 or r3 or r5 in fragment_rev:	
			vntr_count = fragment_rev.count(r1) + \
			fragment_rev.count(r3) + fragment_rev.count(r5)
			MLVA_code.append(vntr_count)
		else:	
			vntr_count = fragment_rev.count(r2) + \
			fragment_rev.count(r4) + fragment_rev.count(r6)
			MLVA_code.append(vntr_count)
		print >>seq, ">%s" % (name+"_" +head), "\n" + fragment_rev
		print >>q, "Fragment locus:", rp1, "-", rp2
		print >>q, "Fragment size:", rDistance	
		print >>q, "Number of repeats:", vntr_count	
	else:
		pass

def VNTR_search(head, p1s, p1e, p2s, p2e, r1, r2):
	print>>q, "\n>>> " + head
	if p1s and p1e in all_dna:
		for x in re.finditer(p1s, all_dna):
			fp1 = int(x.start())
		for y in re.finditer(p1e, all_dna):
			fp2 = int(y.end())
		fDistance = fp2 - fp1
		fragment = all_dna[fp1:fp2]	
		vntr_count = fragment.count(r1) or fragment_rev.count(r2)
		MLVA_code.append(vntr_count)
		print >>seq, ">%s" % (name+"_" +head), "\n" + fragment
		print >>q, "Fragment locus:", fp1, "-", fp2
		print >>q, "Fragment size:", fDistance
		print >>q, "Number of repeats:", vntr_count
	else:
		pass

	if p2s and p2e in all_dna:
		for x in re.finditer(p2s, all_dna):
			rp1 = int(x.start())
		for y in re.finditer(p2e, all_dna):
			rp2 = int(y.end())
		rDistance = rp2 - rp1
		fragment = all_dna[rp1:rp2]
		replace = {"A": "T", "T": "A", "G": "C", "C": "G"}
		replace = dict((re.escape(k), v) for k, v in replace.iteritems())
		pattern = re.compile("|".join(replace.keys()))
		rev_text = pattern.sub(lambda m: replace[re.escape(m.group(0))], fragment)
		fragment_rev = rev_text[::-1]
		vntr_count = fragment_rev.count(r1) or fragment_rev.count(r2)
		MLVA_code.append(vntr_count)
		print >>seq, ">%s" % (name+"_" +head), "\n" + fragment_rev
		print >>q, "Fragment locus:", rp1, "-", rp2
		print >>q, "Fragment size:", rDistance
		print >>q, "Number of repeats:", vntr_count
	else:
		pass

def VNTR3_search(p1s, p1e, p2s, p2e, r1, r2):

	V3p1_start = []
	V3p1_end = []
	V3p2_start = []
	V3p2_end = []

	if (p1s in all_dna) and (p2s not in all_dna):
		for x in re.finditer(p1s, all_dna):
			a = int(x.start())
			V3p1_start.append(a)
		for x in re.finditer(p1e, all_dna):
			b = int(x.end())
			V3p1_end.append(b)
		if len(V3p1_start) == 2:
			V3_start_1 = int(V3p1_start[0])
			V3_start_2 = int(V3p1_start[1])
			V3_end_1 = int(V3p1_end[0])
			V3_end_2 = int(V3p1_end[1])
			fDistance1 = V3_end_1 - V3_start_1
			fDistance2 = V3_end_2 - V3_start_2
			fragment1 = all_dna[V3_start_1:V3_end_1]
			fragment2 = all_dna[V3_start_2:V3_end_2]
			vntr_count1 = fragment1.count(r1) or fragment1.count(r2)
			vntr_count2 = fragment2.count(r1) or fragment2.count(r2)
			if vntr_count1 == vntr_count2 and fix == False:
				MLVA_code.append(vntr_count1)
				MLVA_code.append(0)
			elif vntr_count1 < vntr_count2:
				MLVA_code.append(vntr_count1)
				MLVA_code.append(vntr_count2)
			else:
				MLVA_code.append(vntr_count2)
				MLVA_code.append(vntr_count1)
			print >>seq, ">%s" % (name+"_VNTR3a"), "\n" + fragment1
			print >>seq, ">%s" % (name+"_VNTR3b"), "\n" + fragment2
			print >>q, "\n>>> VNTR3a"
			print >>q, "Fragment locus:", V3_start_1, "-", V3_end_1
			print >>q, "Fragment size:", fDistance1
			print >>q, "\n>>> VNTR3b"
			print >>q, "Fragment locus:", V3_start_2, "-", V3_end_2
			print >>q, "Fragment size:", fDistance2
			print >>q, "Number of repeats:", vntr_count

		elif (len(V3p1_start) == 1) and (p2s not in all_dna):
			V3_start_single_1 = int(V3p1_start[0])
			V3_end_single_1 = int(V3p1_end[0])
			fDistance = V3_end_single_1 - V3_start_single_1
			fragment = all_dna[V3_start_single_1:V3_end_single_1]
			vntr_count1 = fragment.count(r1) or fragment.count(r2)
			MLVA_code.append(vntr_count1)
			MLVA_code.append(0)
			print >>seq, ">%s" % (name+"_VNTR3a"), "\n" + fragment
			print >>q, "\n>>> VNTR3a"
			print >>q, "Fragment locus:", V3_start_single_1, "-",\
			V3_end_single_1
			print >>q, "Fragment size:", fDistance
			print >>q, "Number of repeats:", vntr_count
			print >>q, "\n>>> VNTR3b"
			print >>q, "Fragment locus: 0"
			print >>q, "Fragment size: 0"
			print >>q, "Number of repeats: 0"
		else:
			pass

	if (p2s in all_dna) and (p1s not in all_dna):
		for x in re.finditer(p2s, all_dna):
			a = int(x.start())
			V3p2_start.append(a)
		for x in re.finditer(p2e, all_dna):
			b = int(x.end())
			V3p2_end.append(b)
		if len(V3p2_start) == 2:
			V3_start_1 = int(V3p2_start[0])
			V3_start_2 = int(V3p2_start[1])
			V3_end_1 = int(V3p2_end[0])
			V3_end_2 = int(V3p2_end[1])

			fDistance1 = V3_end_1 - V3_start_1
			rDistance2 = V3_end_2 - V3_start_2
			fragment1 = all_dna[V3_start_1:V3_end_1]
			fragment2 = all_dna[V3_start_2:V3_end_2]

			replace = {"A": "T", "T": "A", "G": "C", "C": "G"}
			replace = dict((re.escape(k), v) for k, v in replace.iteritems())
			pattern = re.compile("|".join(replace.keys()))
			rev_text1 = pattern.sub(lambda m: replace[re.escape(m.group(0))], fragment1)
			rev_text2 = pattern.sub(lambda m: replace[re.escape(m.group(0))], fragment2)
			fragment1rev = rev_text1[::-1]
			fragment2rev = rev_text2[::-1]
			
			vntr_count1 = fragment1rev.count(r1) or fragment1rev.count(r2)
			vntr_count2 = fragment2rev.count(r2) or fragment2rev.count(r2)
			if vntr_count1 == vntr_count2 and fix == False:
				MLVA_code.append(vntr_count1)
				MLVA_code.append(0)
			elif vntr_count1 < vntr_count2:
				MLVA_code.append(vntr_count1)
				MLVA_code.append(vntr_count2)
			else:
				MLVA_code.append(vntr_count2)
				MLVA_code.append(vntr_count1)
			print >>seq, ">%s" % (name+"_VNTR3a"), "\n" + fragment1rev
			print >>seq, ">%s" % (name+"_VNTR3b"), "\n" + fragment2rev
			print >>q, "\n>>> VNTR3a"
			print >>q, "Fragment locus:", V3_start_1, "-", V3_end_1
			print >>q, "Fragment size:", fDistance1
			print >>q, "Number of repeats:", vntr_count
			print >>q, "\n>>> VNTR3b"
			print >>q, "Fragment locus:", V3_start_2, "-", V3_end_2
			print >>q, "Fragment size:", fDistance2	
			print >>q, "Number of repeats:", vntr_count
		elif (len(V3p2_start) == 1) and (p1s not in all_dna):
			V3_start_single_2 = int(V3p2_start[0])
			V3_end_single_2 = int(V3p2_end[0])
			rDistance = V3_end_single_2 - V3_start_single_2
			fragment = all_dna[V3_start_single_2:V3_end_single_2]
			replace = {"A": "T", "T": "A", "G": "C", "C": "G"}
			replace = dict((re.escape(k), v) for k, v in replace.iteritems())
			pattern = re.compile("|".join(replace.keys()))
			rev_text1 = pattern.sub(lambda m: replace[re.escape(m.group(0))], fragment)
			fragment_rev = rev_text1[::-1]
			vntr_count2 = fragment_rev.count(r1) or fragment_rev.count(r2)
			MLVA_code.append(vntr_count2)
			MLVA_code.append(0)
			print >>seq, ">%s" % (name+"_VNTR3a"), "\n" + fragment_rev
			print >>q, "\n>>> VNTR3a"
			print >>q, "Fragment locus:", V3_start_single_2, "-",\
			V3_end_single_2
			print >>q, "Fragment size:", rDistance
			print >>q, "Number of repeats:", vntr_count
			print >>q, "\n>>> VNTR3b"
			print >>q, "Fragment locus: 0"
			print >>q, "Fragment size: 0"
			print >>q, "Number of repeats: 0"
			
		else:
			pass
			
	if (p1s in all_dna) and (p2s in all_dna):
		for x in re.finditer(p1s, all_dna):
			a = int(x.start())
			V3p1_start.append(a)
		for x in re.finditer(p1e, all_dna):
			b = int(x.end())
			V3p1_end.append(b)
		for x in re.finditer(p2s, all_dna):
			a = int(x.start())
			V3p2_start.append(a)
		for x in re.finditer(p2e, all_dna):
			b = int(x.end())
			V3p2_end.append(b)
		V3_start_1 = int(V3p1_start[0])
		V3_end_1 = int(V3p1_end[0])
		V3_start_2 = int(V3p2_start[0])
		V3_end_2 = int(V3p2_end[0])
		fDistance = V3_end_1 - V3_start_1
		rDistance = V3_end_2 - V3_start_2
		fragment1 = all_dna[V3_start_1:V3_end_1]
		fragment2 = all_dna[V3_start_2:V3_end_2]

		replace = {"A": "T", "T": "A", "G": "C", "C": "G"}
		replace = dict((re.escape(k), v) for k, v in replace.iteritems())
		pattern = re.compile("|".join(replace.keys()))
		rev_text2 = pattern.sub(lambda m: replace[re.escape(m.group(0))], fragment2)
		fragment2rev = rev_text2[::-1]

		vntr_count1 = fragment1.count(r1) or fragment1.count(r2)
		vntr_count2 = fragment2.count(r1) or fragment2.count(r2)
		if vntr_count1 == vntr_count2 and fix == False:
			MLVA_code.append(vntr_count1)
			MLVA_code.append(0)
		elif vntr_count1 < vntr_count2:
			MLVA_code.append(vntr_count1)
			MLVA_code.append(vntr_count2)
		else:
			MLVA_code.append(vntr_count2)
			MLVA_code.append(vntr_count1)
		print >>seq, ">%s" % (name+"_VNTR3a"), "\n" + fragment2rev
		print >>q, "\n>>> VNTR3a"
		print >>q, "Fragment locus:", V3_start_1, "-", V3_end_1
		print >>q, "Fragment size:", fDistance
		print >>q, "Number of repeats:", vntr_count1
		print >>q, "\n>>> VNTR3b"
		print >>q, "Fragment locus:", V3_start_2, "-", V3_end_2
		print >>q, "Fragment size:", rDistance
		print >>q, "Number of repeats:", vntr_count2
	else:
		pass		

#------------------------------------------------------------------------------

VNTR1_search(head = "VNTR1", p1s = V1p1s,p1e = V1p1e, p2s = V1p2s, \
			p2e = V1p2e, r1 = VNTR1wt1,r2 = VNTR1wt2, r3 = VNTR1mut1_13, \
			r4 = VNTR1mut2_13, r5 = VNTR1mut1_1315,r6 = VNTR1mut2_1315)

VNTR3_search(p1s = V3p1s, p1e = V3p1e, p2s = V3p2s, p2e = V3p2e, \
			r1 = VNTR3_1, r2 = VNTR3_2)

VNTR_search(head = "VNTR4", p1s = V4p1s, p1e = V4p1e, p2s = V4p2s, \
			p2e = V4p2e, r1 = VNTR4_1, r2 = VNTR4_2)

VNTR_search(head = "VNTR5", p1s = V5p1s, p1e = V5p1e, p2s = V5p2s, \
			p2e = V5p2e, r1 = VNTR5_1, r2 = VNTR5_2)

VNTR_search(head = "VNTR6", p1s = V6p1s, p1e = V6p1e, p2s = V6p2s, \
			p2e = V6p2e, r1 = VNTR6_1, r2 = VNTR6_2)

#------------------------------------------------------------------------------

#Send MLVA number to output file
searchfile = open("%s/MLVAdatabase.txt" % mlvaDir, 'r')

for line in searchfile:
	s = open("%s/MLVA/summary.txt" % outDir, "a")
	if str(MLVA_code) in line:
		print >>s, name,':', line.rstrip('\n')
		print >>q, "\n>>> MLVA Type: ", line.rstrip('\n')
		test = True
	else:
		test = False

if test == False:
	print >>s, name,':', "MLVA Type not Found"
	print >>q, "\n>>> MLVA Type: MLVA Type not Found"

if keep == False:
	print ">>> Removing temperary files..."
	os.remove("%s/MLVA/MLVA_data/%s" % (outDir, outFile))
else:
	pass
print ">>> Search Complete."
