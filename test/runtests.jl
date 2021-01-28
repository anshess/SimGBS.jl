using simGBS
using Test

@testset "simGBS.jl" begin
## parameters
## step 1. genome
re = [ApeKI]; # specify the restriction enzyme to be used for virtual digestion
useChr = [1]; # specify either the number of chromosome(s) or a set of chromosome(s) to be simulated
useChrLen = Array{Float64}(undef,0); # specify the length of chromsome(s) in cM to be simulated, otherwise using the entire chromosome
lower = 65; # lower bound of fragment size selection
upper = 195; # upper bound of fragment size selection

## step 2. population structure
numFounders = 100; # number of founders in the base population
endSize = 500; # number of individuals to end up in the changingPopSize step
numInd = 96; # number of individuals to be simulated (default = 96)
numGenCha = 20; # number of generations for changingPopSize function (default = 20)
numGenCon = 100; # number of generations for constantPopSize function (default = 100)
numGenFinal = 4; # number of final generations to be used to select individual (default = 4)
useWeights = Array{Float64}(undef,0) # weights of each of contributing genetarion in the fianal population composition (default = false)
usePedigree = true;  # false if you don't use pedigree, otherwise specify the pedigree file to be used
pedFile = "sim.ped"; # file stores the pedigree (default = "sim.ped")
pedOutput = false; # true if print out pedigree (default = false)

## step 3. data
totalQTL = 1000 # total (across all chromosomes) number of QTLs to be simulated
totalSNP = 0 # total (across all chromosomes) number of QTLs to be simulated, set [totalSNP = 0] if sampling based on density
muDensity = 0.05 # expected gloabl SNP density
sigmasqDensity = 0.0001 # varaince of global SNP density
winSize = 1000000 # window size to be used to sample SNP positions
muAlleleFreq = 0.175 # expected allele frequency
sigmasqAlleleFreq = 0.0065 # varaince of allele frequency
seqDepth = 10 # expected sequencing depth
barcodeFile = "GBS_Barcodes.txt" # file contains GBS barcodes

## miscellaneous
plotOutput = true; # true if plots are required
writeOutput = true; # true if outputs in text file are required



## simulation
## step 1. geneate GBS fragments
# virtualDigest(genofile, re, whichChr, chrLen, lower ,upper, plotOutput, writeOutput);
digestGenome(genofile, re, useChr, useChrLen, lower ,upper, plotOutput, writeOutput);

## step 2. define population
definePopulation(numFounders, endSize, numGenCha, numGenCon, numGenFinal, numInd, useWeights, usePedigree, pedFile, pedOutput);

## step 3. simulate GBS data
simGBS(totalQTL, totalSNP, muDensity, sigmasqDensity, winSize, muAlleleFreq,sigmasqAlleleFreq, re, seqDepth, barcodeFile, plotOutput, writeOutput)

end



