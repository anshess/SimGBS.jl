using SimGBS

## parameters
### Step One: Generate GBS Fragments
genofile =  "ref.fa.gz"
re = [SimGBS.ApeKI]; # specify the restriction enzyme to be used for virtual digestion
useChr = [1]; # specify either the number of chromosome(s) or a set of chromosome(s) to be used for simulating GBS data
useChrLen = Array{Float64}(undef,0); # specify the length of chromsome(s) in cM to be simulated, otherwise using the entire chromosome
lower = 65; # lower bound of fragment size selection
upper = 195; # upper bound of fragment size selection
winSize = 1000000 # window size to be used to sample SNP density

### Step Two: Define Population Structure
numFounders = 100; # number of founders in the base population
endSize = 1000; # number of individuals to end up in the changingPopSize step
numGenCha = 20; # number of generations for changingPopSize function 
numGenCon = 50; # number of generations for constantPopSize function 
numGenFinal = 4; # number of final generations to be used to select individual 
numInd = 96; # number of individuals to be simulated
useWeights = Array{Float64}(undef,0) # weights of each of contributing genetarion in the fianal population composition 
usePedigree = false;  # false if you don't use pedigree, otherwise specify the pedigree file to be used
pedFile = "newPed.ped"; # file stores the pedigree (default = "sim.ped")
pedOutput = false; # true if print out pedigree (default = false)

### Step Three: Simulate GBS Process
totalQTL = 100 # total (across all chromosomes) number of QTLs to be simulated
totalSNP = 0 # total (across all chromosomes) number of QTLs to be simulated, set [totalSNP = 0] if sampling based on density
muDensity = -2.0 # location parameter of log-Laplace distribution (for sampling SNP density)
sigmasqDensity = 0.001 # scale parameter of log-Laplace distribution (for sampling SNP density)
winSize = 1000000 # window size to be used to sample SNP density
muAlleleFreq = 0.5 # mean of sampled allele frequency
sigmasqAlleleFreq = 0.0061 # variance of sampled allele frequency
meanDepth = 20.0 # expected sequencing depth
barcodeFile = "GBS_Barcodes.txt" # file contains GBS barcodes
useChr = [1]; # specify either the number of chromosome(s) or a set of chromosome(s) to be used for simulating GBS data 

### miscellaneous
plotOutput = false; # true if plots are required
writeOutput = true; # true if outputs in text file are required
onlyOutputGBS = true; # true if only GBS data is required
 
## Run SimGBS
### Step One: Generate GBS Fragments
@time digestGenome(genofile, re, useChr, useChrLen, lower ,upper, winSize, plotOutput, writeOutput);

### Step Two: Define Population Structure
@time definePopulation(numFounders, endSize, numGenCha, numGenCon, numGenFinal, numInd, useWeights, usePedigree, pedFile, pedOutput);

### Step Three: Simulate GBS Process
@time GBS(totalQTL, totalSNP, muDensity, sigmasqDensity, winSize, muAlleleFreq, sigmasqAlleleFreq, re, meanDepth, barcodeFile, useChr, plotOutput, writeOutput, onlyOutputGBS);