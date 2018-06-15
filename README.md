## whole genome sequence MLVA Program (wgsMLVA) version 1.0
Created by Dane Kania
06-14-2018

### ---What is this?---

The whole genome sequence Multiple Locus Variable Number Tandem Repeat Program (wgsMLVA) aims to recreate PCR-based MLVA typing using closed *Bordetella pertussis* genomes. It is currently limited to closed B. pertussis genomes with future plans to expand to draft genomes.


### ---How does it work?---

The program identifies six primer sites, labeled Variable Number Tandem Repeat (VNTR) sites, that each hold a set number of repeats. Each VNTR site has a unique, repeated DNA element that is counted by the program. This is done with all six VNTR sites, creating a list, similar to a barcode. wgsMLVA compares the specific list of VNTR counts to a database of possible MLVA types (included in this git page as MLVAdatabase.txt) and prints out the corrisponding MLVA type. The length of the fragment's flanking regions are also taken into account, to insure known consistanty across the VNTR's of each MLVA type.

The primers were developed by [Schouls et alia 2004](https://www.ncbi.nlm.nih.gov/pubmed/15292152).

The database was used in this program was downloaded and based on the [mlva.net](https://www.mlva.net/bpertussis/default.asp) database.

### ---How do I run wgsMLVA?---

Example: 
```
python wgsMLVA.py [other arguments] -in [input fasta file] -out [output directory]
```
Dependencies:

python 2.7 is needed to run.

wgsMLVA arguments:
```
-h or --help : list of arguments
		Prints out list of command line arguments.

-in or --input : fasta file
		Closed genome of a Bordetella pertussis sample. Required.

-out or --output : output directory.
		Optional, will use current directory if none given.

-n or --names : set output file names
		Optional, the program will use the name of the 
		original file if left blank.

-v or --version : Program Version		
		Prints out the program verion to terminal.	

-u or --unmask : Removes VNTR3b masking, shows actual repeat value.
		See README for more details.
```

Data output:

The program will create a new "MLVA" directory with two subdirectories ("MLVA_data" and "sequence_data") and onesummary file (summary.txt). MLVA output data is in "MLVA_data" and raw sequence information for each VNTRsite is stored in "sequence data".

examples MLVA type and VNTR count list output:

MLVA type = [VNTR1, VNTR3a, VNTR3b, VNTR4, VNTR5, VNTR6]

27 = [8, 7, 0, 7, 6, 7]
 
### ---Caveats and future work---

Traditional PCR MLVA is a less accurate than wgsMLVA. It estimates the number of repeats in each VNTR site by molecular weight, with larger PCR products having more repeats than smaller ones. This makes matching wgsMLVA and PCR MLVA difficult at times. Here are a list of various points of conflicts I have encountered developing wgsMLVA v1.0 and the remedies I have used. this list will update overtime as I find new issues and fix old ones.

#### I. VNTR3 is capable of duplicating in a B. pertussis genome and PCR MLVA has difficuties with this.

In the example lists above you may notice there are two VNTR3's. This is because VNTR3 has a chance to duplciate in the genome. These duplications can also have different numbers of repeat monomers too. However, because PCR MLVA relys on molecular weight and gel electrophoresus to esitmate the VNTR fragments there is no easy way to distinguish the VNTR3 duplications unless they are drastically different in size (eg. different number of repeats). wgsMLVA does not have this problem and can identify both VNTR3 fragments every time.

This produces an discrpency moving forward with wgsMLVA. wgsMLVA can literally produce more accurate data, but the data conflicts with standardized B. pertussis MLVA nomenclature. Currently, the program ignores this extra information. If the program finds a duplicate and they both have the same repeat count, the duplicate is changed to a zero (exactly how PCR MLVA works). Further work should be done to strengthen this argument too. 

An argument may be added in the program, -u or --unmask, which removes this, createing a more "real" MLVA type that accounts for 2 VNTR3's of the same length. Currently, this is not used because it would dratically change how the typing system works, so until these two schemes can be reconciled, the program will continue to use the standard MLVA typing system.

#### II. Variability in repeats.

A typical VNTR1 fragment for MLVA type 27 looks like this:

>VNTR1 
AAAATTGCGGCATGTGGGCTGACTCTGAAAGCGATGCTCACGAAAAGGGAACGCGGCGCCGTCGGGCGCCGCGCGCCGCTTAGGACTGCTGGCCTGCGGCCGGCGCCTGCTTGGCGGGTTCCTGCTTGGCGGGTTCCTGCTTGGCGGGCTCCTGCTTGGCGGGTTCCTGCTTGGCGGGTTCCTGCTTGGCGGGTTCCTGCTTGGCGGGCTGCTGCTTGGCGGGCTGCTGGGCCGGCGCCTGCTGGCCAGGAGCGGGCTGCTGGCCGGCAGGCGCCGCGCCCCCCTTGTTCCAGGGCGAAGCCTGCACCGGCGCCCCCGGACGGATCTTCTGGAAGCCTTCGACCACCACCACGTCTCCCGCCGCCAGG

Type 27's VNTR1 has 8 repeats. When looking for VNTR1's repeat monomer (CTGCTTGGCGGGTTC), there are only 5 exact matches. In order to get to 8, I hypothesize that there are repeat mutants that exist in the fragment. When I lower the stringency, I find the three other repeats, one with one SNP at base pair 13 and two with two SNPs at base pair 13 and 15. They also exist in between the wild type repeats (repeats 3 and 8), giving some crediblity to this hypothesis. 

>VNTR1 - Repeat variataion
CTGCTTGGCGGGTTC CTGCTTGGCGGGTTC
CTGCTTGGCGGG**C**TC CTGCTTGGCGGGTTC 
CTGCTTGGCGGG**C**T**G** CTGCTTGGCGGG**C**T**G**

While not confirmable, this has been the case in all VNTR1's tested with this program (over 500 closed and draft genomes). This may be one area where a length based method is better than an exact match, as the specificity of the repeats are not required. While the fact that the length of the fragment produces the same VNTR repeat count as my method of including these mutant repeats makes for a compelling arguement to include them, it may just be a coincidence. However, this program is confirmed on 464 closed genomes with corrisponding PCR MLVA data and has a 100% confirmation rate. This possiblility of repeat mutation should be considered going forward.
