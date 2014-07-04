#!/usr/bin/python

#Author: Daniel Gomez-Sanchez <daniel.gomez.sanchez@hotmail.es>

#Create a tped and tfam with a vcf file (tab delimited)
#The missing data is 0 0
#The reference homozygote is 1 1
#The alternative homozygote is 2 2
#The heterozygote is 1 2

#Needs argparse-1.2.1.

#Returns the genotypes in the format "0 0", "1 1". "2 2", "1 2"
def geno(genotype, coding):
	'''genotype is a string'''
	if genotype == "0/0" or genotype == "0|0":
		code=[coding[0], coding[0]]
		return " ".join(code)
	elif genotype == "1/1" or genotype == "1|1":
		if len(coding[1]) > 1:
			return False
		else:
			code=[coding[1], coding[1]]
			return " ".join(code)
	elif genotype == "1/0" or genotype == "0/1" or genotype == "1|0" or genotype == "0|1":
		if len(coding[1]) > 1:
			return False
		else:
			code=[coding[0], coding[1]]
			return " ".join(code)
	elif genotype == "./." or genotype == ".|.":
		return "0 0"
	else:
		return False

def change_geno(genotypes, position, coding):
	'''genotype is a list, position a number
	coding is a list with the ref and alt alleles'''
	genotype=[]
	for n in genotypes:
		original = n.split(":")[position]
		genot = geno(original, coding)
		if genot:
			genotype.append(genot)
		else:
			return False
	return genotype
			

def print_tfam(header, handle, printed):
	#Store the sample names
	names = header.rstrip("\n").split("\t")[9:]
	#For every name, print in the filehandle the line with the tfam format
	if printed:
		for name in names:
			print >> handle, str(name)+"\t"+str(name)+"\t0\t0\t0\t-9"
	return len(names)

def print_tped(line, handle, coding, names, printed):
	'''coding is "12" or "ATCG"'''
	#Store the line as a list
	line = line.rstrip("\n").split("\t")
	#If there are an ID for the SNP
	if names:
		#Store it
		name = line[2]
	else:
		#If not, put a name as the chromosome, dot, position
		name = str(line[0])+"."+str(line[1])
	#Store the chromosome
	chrom = line[0].strip("chr")
	#Store the position
	pos = line[1]
	#Store the format
	FORMAT = line[8].split(":")
	#If the coding is "ATCG"
	if coding=="ATCG":
		#Store the ref and the alt
		refalt=line[3:5]
	else:
		refalt=["1", "2"]
	#Look in the format for the position for the GT
	GT=-1
	for n in xrange(len(FORMAT)):
		if FORMAT[n] == "GT":
			GT=n
			break
	#If there are GT
	if GT==-1:
		return name, 2
	else:
		#Store the genotypes in tped format
		genotypes = change_geno(line[9:], GT, refalt)
		#If there are correct genotypes for all
		if genotypes:
			if printed:
				print >> handle, str(chrom)+"\t"+str(name)+"\t0\t"+str(pos)+"\t"+"\t".join(genotypes)
			return "Correct", 0
		else:
			return name, 1
		

	

if __name__ == '__main__':
	#module to parse the commands
	import argparse, sys, time
	#create the parser with the arguments and the information
	parser = argparse.ArgumentParser(description="Transform a VCF with only SNPs to a TPED and TFAM PLINK format.\nThe file must be tab-delimited.\nNote that non-bialelic SNPS and/or whitout GT information will be removed.\nRecode option to print out the result file.")
	#args.input have the input name
	parser.add_argument("--vcf", dest="input", metavar="{fileroot}", help="Specify VCF file", required=True)
	#args.out have the output name
	parser.add_argument("--out", dest="out", metavar="{fileroot}", help="Specify output root filename")
	#recode options, if they are not set, no print the new files	
	parser.add_argument("--recode", action="store_true", help="Output new TPED and TMAP files")
	parser.add_argument("--recode12", action="store_true", help="As above, with 1/2 for reference and alternative alleles")
	#args.table if there are IDs set
	parser.add_argument("--ids", action="store_true", help="ID column is specified")
	#arg.Nchr contain a list of chromosomes to exclude
	parser.add_argument("--non-chr", default=[], metavar="chrN", dest="Nchr", help="Exclude chromosome (allowed multiple times)", action="append")
	


	#parsing the commands and print in the standar error the information
	args = parser.parse_args()

	if not args.out:
		args.out=args.input

	sys.stderr.write("@----------------------------------------------------------@\n")
	sys.stderr.write("|     VCF2tPLINK     |     v.1.1     |     10/Apr/2014     |\n")
	sys.stderr.write("|----------------------------------------------------------|\n")
	sys.stderr.write("|               Author: Daniel Gomez-Sanchez               |\n")
	sys.stderr.write("|            <daniel.gomez.sanchez@hotmail.es>             |\n")
	sys.stderr.write("@----------------------------------------------------------@\n")
	sys.stderr.write("\nAnalysis started: "+time.strftime("%c")+"\n")
	sys.stderr.write("\nOptions in effect:\n")
	sys.stderr.write("\t--vcf\t"+args.input+"\n")
	if args.recode:
		sys.stderr.write("\t--out\t"+args.out+"\n")
		sys.stderr.write("\t--recode\n")
		formato="ATCG"
		recode=True
	elif args.recode12:
		sys.stderr.write("\t--out\t"+args.out+"\n")
		sys.stderr.write("\t--recode12\n")
		formato="12"
		recode=True
	else:
		args.recode=False
		args.recode12=False
		formato="12"
		recode=False
	if args.Nchr:
		sys.stderr.write("\t--non-chr\t"+" ".join(args.Nchr)+"\n")
	if args.ids:
		sys.stderr.write("\t--ids\n")
	else:
		args.ids=False

	sys.stderr.write("\n** For gPLINK compatibility, do not use '.' in --out **\n")

	tfam=str(args.out)+".tfam"
	tped=str(args.out)+".tped"
	vcf=str(args.input)+".vcf"
	snperr=str(args.out)+".snps.failed"		
		
	excluded=0
	correct=0
	error1=0
	error2=0
	samples=0

	#Open the input and tped output file
	INPUT = open(vcf, "r")
	if recode:
		TPED=open(tped, "w")
		TFAM=open(tfam, "w")
	else:
		TPED=None
		TFAM=None

	#For every line
	for line in INPUT:
		#If the line is the header
		if "#CHROM" in line:
			#Write the tfam
			samples = print_tfam(line, TFAM, recode)
		elif "##" not in line:
			if line.rstrip("\n").split("\t")[0] in args.Nchr:
				log = "excluded"
			else:
				name, log = print_tped(line, TPED, formato, args.ids, recode)
			if log == 0:
				correct+=1
			elif log == 1:
				if error1 == 0 and error2 == 0:
					SNP=open(snperr, "w")
				error1+=1
				print >> SNP, name
			elif log == 2:
				if error1 == 0 and error2 == 0:
					SNP=open(snperr, "w")
				error2+=1
				print >> SNP, name
			elif log == "excluded":
				excluded+=1

	INPUT.close()
	if error1 > 0 or error2 > 0:
		SNP.close()

	total=correct+error1+error2+excluded
	sys.stderr.write(str(correct)+" (of "+str(total)+") markers to be included from [ "+str(vcf)+" ]\n")
	sys.stderr.write("Warning, found "+str(samples)+" individuals with ambiguous sex codes\n")
	sys.stderr.write(str(samples)+" individuals read from [ "+str(vcf)+" ]\n")
	sys.stderr.write("0 individuals with nonmissing phenotypes\n")
	sys.stderr.write("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n")
	sys.stderr.write("Missing phenotype value is also -9\n")
	sys.stderr.write("0 cases, 0 controls and "+str(samples)+" missing\n")
	sys.stderr.write("0 males, 0 females, and "+str(samples)+" of unspecified sex\n")
	if error2 > 0:
		sys.stderr.write("Warning, "+str(error2)+" SNPs without GT information\n")
	if error1 > 0:
		sys.stderr.write("Warning, "+str(error1)+" non-bialelic SNPs\n")
	if error2 > 0 or error1 >0:
		sys.stderr.write("Writing these SNPs to [ "+str(snperr)+" ]\n")
	if args.Nchr:
		sys.stderr.write("Excluded "+str(excluded)+" SNPs from "+", ".join(args.Nchr)+"\n")
	sys.stderr.write("After conversion, there are "+str(correct)+" SNPs\n")
	if recode:
		sys.stderr.write("Writing recoded tped file to [ "+str(tped)+" ]\n")
		TPED.close()
		sys.stderr.write("Writing recoded tfam file to [ "+str(tfam)+" ]\n")
		TFAM.close()
	sys.stderr.write("\nAnalysis finished: "+time.strftime("%c")+"\n")
