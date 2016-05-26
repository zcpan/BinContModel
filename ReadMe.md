## Simultaneous modeling of disease status and clinical phenotypes to increase power in GWAS

pyslmGWAS, a method method for simultaneous association testing of genetic variants with both case status and a clinical covariate that resolves both of these issues.

### Tutorial

  Usage: python pyslmGWAS.py --tfile plinkFile --f population_frequency --phenofile plinkFormattedPhenotypeFile resultFile (result file should be indicated in the command line)
  Get beta1, beta2 estimation from bincontmodel. 
  
  Options:
  -v, --verbose <string>    Print extra info
  
  --tfile       <string>    the base for a PLINK tped file
  
  --phenofile   <string>    phenofile which contains continual phenotypes and discrete disease status
                            if you do not set the parameters, it would automatically search the file 
                            whose name starts with the base for a PLINK tped file, and ends with ".phenos"
                            
  --f           <int>       the population frequency
  
  Here is the toy example to show how to run the program
  
  $ python pyslmGWAS.py -v --tfile ./data/simulate --f 0.5 ./data/test.output.txt
  
  
### Simulation

  I also provide the simulation script to further validation our methods. 
  
  Usage: python creatSimulateData.py --power 0.5 --n 5000 --f 0.4 --maf 0.2 --base ./data/simulate --cor 0.2 --ez 0.3
  --power  <float>    The power you want to set based on univariate test
  --n      <int>      The total number of case and control studies
  --f      <float>    The population frequency
  --maf    <float>    Minor Allele Frequency
  --base   <string>   The file base for simulation data
  --cor    <float>    The correlation of phenotypes data and disease liability
  --ez     <float>    The effect size of phenotypes on SNPs
  
  $ python creatSimulateData.py --power 0.5 --n 5000 --f 0.4 --maf 0.2 --base ./data/simulate --cor 0.2 --ez 0.3
  