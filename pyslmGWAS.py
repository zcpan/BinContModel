import pdb,time,sys
from scipy.stats import mode
def printOutHead():
	out.write("\t".join(["SNP_ID","BETA_1","BETA_2","P_VALUE(beta2)"]) + "\n")
def outputResult(id,beta1,beta2,pvalue):
	out.write("\t".join([str(x) for x in [id,beta1,beta2,pvalue]]) + "\n")
from optparse import OptionParser,OptionGroup
usage = """
This program provides basic genome-wide association (GWAS) functionality.  
You provide a phenotype and genotype file and the program outputs a result file with information about each SNP, including the association p-value. 
The input file are all standard plink formatted with the first two columns specifiying the individual and family ID.  
For the phenotype file, we accept either NA or -9 to denote missing values.  

Usage:

      python pyslmGWAS.py --tfile plinkFile --f population_frequency --phenofile plinkFormattedPhenotypeFile resultFile

"""
parser = OptionParser(usage=usage)
basicGroup = OptionGroup(parser,"Basic Options")
advancedGroup = OptionGroup(parser,"Advanced Options")

basicGroup.add_option("--tfile",dest="tfile",help="The base for a PLINK tped file")
basicGroup.add_option("--phenofile", dest="phenoFile",default=None)
basicGroup.add_option("--f",dest="population_frequency",default=None)
advancedGroup.add_option("-v", "--verbose",action="store_true", dest="verbose", default=False,help="Print extra info")

parser.add_option_group(basicGroup)
parser.add_option_group(advancedGroup)

(options, args) = parser.parse_args()

import os
import numpy as np
from scipy import linalg
from slm.slm import SLM
from slm import input
import random

if len(args) != 1:  
   parser.print_help()
   sys.exit()

outFile = args[0]

if not options.tfile:
	parser.error("You must provide PLINK input file base (--tfile)")
if not options.population_frequency:
	parser.error("You must define the population frequency (--f)")
if options.verbose:
	sys.stderr.write("Reading SNP input...\n")
f = float(options.population_frequency)
IN = input.plink(options.tfile,phenoFile=options.phenoFile,f = f,normGenotype=True)

if not os.path.isfile(options.phenoFile or IN.fbase + '.phenos'):
   parser.error("No .pheno file exist for %s.  Please provide a phenotype file using the --phenofile" % (options.phenoFile or IN.fbase + '.phenos'))
if options.phenoFile and os.path.isfile(options.phenoFile):
	IN.getPhenos(options.phenoFile)
else:
	IN.getPhenos(IN.fbase + '.phenos')
count = 0
out = open(outFile,'w')
printOutHead()
n = len(IN.keep)

Y1 = IN.phenos[:,0].reshape((n,1))
Y2d = IN.phenos[:,1].reshape((n,1))
Y1[np.isnan(Y1[:,0]),0] = Y1[True - np.isnan(Y1[:,0]),0].mean()
Y2d[np.isnan(Y2d[:,0]),0] = int(mode(Y2d[True - np.isnan(Y2d[:,0]),0])[0][0])

for snp, id in IN:
	count += 1
	if options.verbose and count % 1000 == 0:
		sys.stderr.write("At SNP %d\n" % count)
	x = snp.reshape((n,1))
	Ls = SLM(Y1,Y2d,x,f,verbose=options.verbose)
	Ls.fit()
	outputResult(id,Ls.beta1,Ls.beta2_norm,Ls.pvalue)
out.close()
