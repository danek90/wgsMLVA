		in silico MLVA Program (MLVAisP) v1.1
			Created by Dane Kania
			     11/1/2017

	---What is this?---

	The MLVA in silco Program to type Bordetella pertussis isolates through the
	Multiple Locus Variable Number Tandem Repeat typing scheme. Using primers
	created from Schouls et al 2004 [https://www.ncbi.nlm.nih.gov/pubmed/15292152],
	this program counts number of Variable Number Tandem Repeats (VNTR) within six
	primer regions and compares the VNTR barcode (VNTR1, VNTR3, VNTR3b, VNTR4, VNTR5, and VNTR6)
	created to a MLVA database. From the known MLVA database, a MLVA type can be derived. 

	---How to run program---

	Example: 

	python isMLVA.py -i [input fasta file] -o [output directory] 
	
	Dependencies:

	1. python 2.7.X is needed to run. if running on biolinux add the module with this command:

		module add Python/2.7.13

	isMLVA arguments:
	
	-h or --help : list of arguments

		Prints out list of command line arguments.

	-i or --input : fasta file

		Closed genome of a Bordetella pertussis sample. 

	-o or --output : output directory
	
		Directory to create the output files. The program will create a new "MLVA"
		directory with two subdirectories "MLVA_data" and "sequence_data" and one
		summary file (summary.txt).

		MLVA output data is in "MLVA_data" and sequence information for each VNTR
		site is stored in "sequence_data"

		Diagram:

		/output/path/MLVA/
			--> /MLVA_data/
			--> /sequence_data/
			* summary.txt

	-v or --version : Program Version
	
		Prints out the program verion to terminal.		

 
