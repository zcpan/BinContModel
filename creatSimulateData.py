import numpy as np
import random
from scipy.stats import norm

from optparse import OptionParser,OptionGroup
def get_power(N,power):
	none_centra_para = 0
	while abs(norm.cdf(np.sqrt(N)*none_centra_para + norm.ppf(0.0025)) + 1 - norm.cdf(none_centra_para * np.sqrt(N) - norm.ppf(0.0025)) - power) > 1e-3:
		none_centra_para += 0.000001
	return none_centra_para

usage = """

	  python creatSimulateData.py --power 0.5 --n 5000 --f 0.4 --maf 0.2 --base ./data/simulate --cor 0.2 --ez 0.3

"""

parser = OptionParser(usage=usage)
basicGroup = OptionGroup(parser,"Basic Options")
basicGroup.add_option('--power',dest="Power",help="The power you want to set based on univariate test")
basicGroup.add_option('--n',dest="sample_size",help="The total number of case and control studies")
basicGroup.add_option('--f',dest="population_frequency",help="The population frequency")
basicGroup.add_option('--maf',dest='maf',help="Minor Allele Frequency")
basicGroup.add_option('--base',dest='base',help="The file base for simulation files")
basicGroup.add_option('--cor',dest='correlation',help="The correlation of phenotypes data and disease liability")
basicGroup.add_option('--ez',dest='effect_size',help="The effect size of phenotypes on SNPs")
(options, args) = parser.parse_args()

if not options.Power:
	parser.error("You must specific the power")
power = float(options.Power)
if not options.sample_size:
	parser.error('Number of case and control studies should be provided')
N = int(options.sample_size)
if not options.population_frequency:
	parser.error('Population frequency has not been provided')
f = float(options.population_frequency)
if not options.maf:
	parser.error('Minor Allele Frequency is missing')
maf = float(options.maf)
if not options.base:
	parser.error('file base should be provided')
fbase = str(options.base)
if not options.correlation:
	parser.error('Correlation of phenotypes data and disease liability should be provided')
cor = float(options.correlation)
if not options.effect_size:
	parser.error('phenotypes effect size on SNPs should be provided')
beta1 = float(options.effect_size)
beta2 = get_power(N,power)
case, control = N/2, N/2
sampleSize = int(max(case/f,control/(1-f)))
np.random.seed(1)
s = np.random.binomial(2, maf, sampleSize)
norm_s = (s - 2*maf) / np.sqrt(2*maf*(1-maf))
mean = (0,0)
cov =  [[1, cor], [cor, 1]]
x = np.random.multivariate_normal(mean,cov,s.shape[0])
liability_y = norm_s * beta2 + x[:,0]
pheno_y = norm_s * beta1 + x[:,1]
threshold = norm.ppf(1-f)

where = np.where(liability_y > threshold)
neg = np.random.choice(np.where(liability_y <= threshold)[0],len(where[0]),replace=False)
choicelist = []
choicelist.extend(where[0])
choicelist.extend(neg)

liability_y = liability_y[choicelist]
s = s[choicelist]
pheno_y = pheno_y[choicelist]
norm_s = norm_s[choicelist]

File = fbase + '.map'
fileWrite = open(File,'w')
fileWrite.write('1\trs3115860\t0\t743268\n')
fileWrite.close()

File = fbase + '.tfam'
fileWrite = open(File,'w')
for i in range(0,s.shape[0]):
	fileWrite.write(str(1000 + i) + '\t' + str(1000 + i) + '\t' + '0\t0\t1\t1\n')
fileWrite.close()

File = fbase + '.tped'
fileWrite = open(File,'w')
fileWrite.write('1\trs3115860\t0\t743268\t')
for i in range(0,s.shape[0]):
	if s[i] == 0:
		fileWrite.write('1\t1\t')
	elif s[i] == 1:
		fileWrite.write('1\t2\t')
	elif s[i] == 2:
		fileWrite.write('2\t2\t')
	else:
		print 'error...'
fileWrite.close()

File = fbase + '.phenos'
fileWrite = open(File,'w')
fileWrite.write('FID\tIID\ttgres\tstatus\n')
for i in range(0,s.shape[0]):
	if liability_y[i] > threshold:
		fileWrite.write(str(1000 + i) + '\t' + str(1000 + i) + '\t' + str(pheno_y[i]) + '\t'+'1\n')
	else:
		fileWrite.write(str(1000 + i) + '\t' + str(1000 + i) + '\t' + str(pheno_y[i]) + '\t'+'0\n')
fileWrite.close()



