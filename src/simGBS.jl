module simGBS

using Compat, Random, Distributions, DelimitedFiles, GZip, StatsBase, Distributed, SharedArrays, Plots

## 01.Genome
## This file contains functions used for virtual digestioon in simGBS
## define restriction enzyme [name, cut-sequence, cut-site, overhang]
mutable struct restrictionEnzyme
    name::String # name of the restriction enzyme
    cutSeq::Array{String} # recognition sequence
    cutPlace::Array{Int64} # where it gets cut
    overhang::Array{String} # overhang of the restriction enzyme if exists
    restrictionEnzyme(name, cutSeq, cutPlace, overhang) = new(name, cutSeq, cutPlace, overhang)
end;

## list of restraiction enzymes
ApeKI = restrictionEnzyme("ApeKI", ["GCAGC", "GCTGC"], [2, 2], ["CTG", "CAG"]); # Elshire et al. 2011
PstI = restrictionEnzyme("PstI", ["CTGCAG"], [2], ["CTGCA"]); # Poland et al. 2012
MspI = restrictionEnzyme("MspI", ["CCGG"],[2],["CGG"]); # Poland et al. 2012


## nucleotides and its reverse complements
nucl = Dict('R' => rand(['A', 'G']), 'Y' => rand(['C', 'T']), 'S' => rand(['G', 'C']), 'W' => rand(['A', 'T']),
    'K' => rand(['G', 'T']), 'M' => rand(['A', 'C']), 'B' => rand(['C', 'G', 'T']), 'D' => rand(['A', 'G', 'T']),
    'H' => rand(['A', 'C', 'T']), 'V' => rand(['A', 'C', 'G'])); # randomly assign IUPAC codes for uncertain ones in the assembly
comp = Dict('A' => 'T', 'C' => 'G','G' => 'C','T' => 'A'); # complements of DNA nucleotides

## split sequence
splitSequence = function(seq, cutSeq, startSeq, endSeq)
    cuts = split(uppercase(seq), cutSeq) # split DNA sequence
    if length(cuts) > 1
        cuts[2:length(cuts)] = startSeq .* cuts[2:length(cuts)] # re-attach the strat seq to the fragment after cut-place
        cuts[1:(length(cuts)-1)] = cuts[1:(length(cuts)-1)] .* endSeq # re-attach the end seq to the fragment before cut-place
    else
        cuts
    end
    cuts
end;

## digest sequence with restriction enzyme
digest = function(sequence, re::restrictionEnzyme)
    cutSeq = re.cutSeq
    cutPlace = re.cutPlace
    overhang = re.overhang
    startSeq = [SubString(cutSeq[i], cutPlace[i]) for i = 1:length(cutSeq)]
    endSeq = [SubString(cutSeq[i], 1, cutPlace[i] - 1) * overhang[i] for i = 1:length(cutSeq)]
    for k = 1:length(cutSeq)
        cuts = sequence
        sequence = []
        for s in cuts
            sequence = [sequence; splitSequence(s, cutSeq[k], startSeq[k], endSeq[k])]
        end
        sequence
    end
    if length(sequence) > 2
        sequence[2:length(sequence)]
    else
        sequence
    end
end;

## function for double-digestion
digestSecond = function(sequence, re::restrictionEnzyme)
    cutSeq = re.cutSeq
    cutPlace = re.cutPlace
    overhang = re.overhang
    startSeq = [SubString(cutSeq[i], cutPlace[i]) for i = 1:length(cutSeq)]
    endSeq = [SubString(cutSeq[i], 1, cutPlace[i] - 1) .* overhang[i] for i = 1:length(cutSeq)]
    goodSeq = []
    for k = 1:length(sequence)
        s = sequence[k]
        sp = splitSequence(uppercase(s), cutSeq[1], startSeq[1], endSeq[1])
        if length(sp) > 1
            if k == 1
                s2 = startSeq .* sp[length(sp)]
                goodSeq = [goodSeq; s2]
            elseif k == length(sequence)
                s1 = sp[1] .* endSeq
                goodSeq = [goodSeq; s1]
            else
                s1 = sp[1] .* endSeq
                s2 = startSeq .* sp[length(sp)]
                goodSeq = [goodSeq; s1; s2]
            end
        else
            goodSeq = [goodSeq; sp]
        end
    end
    goodSeq
end;

## virtual digestion
digestGenome = function(genofile::String, re::Array{restrictionEnzyme,1}, useChr::Array{Int64}, useChrLen::Array{Float64}, lowerThresh::Int64, upperThresh::Int64, plotOutput::Bool , writeOutput::Bool)
    ## module 1: read genome and cut it into fragments
    genome = readdlm(GZip.open(genofile), '\t') # load genome
    index = LinearIndices(genome) # index genome
    header = index[findall(x -> x[1] == '>', genome)] # identify headers within each FASTA file
    global chrLen = Array{Float64}(undef,0) # empty array to store estimated chromosome length
    seqFrag = Array{String}(undef,0) # empty array to store the DNA sequence
    chrFrag = Array{Int64}(undef,0) # empty array to store chromosome
    idFrag = Array{Int64}(undef,0) # empty array to store the ID
    posFrag = Array{Int64}(undef,0) # empty array to store the position
    lenFrag = Array{Int64}(undef,0) # empty array to store the position
    ## (optional) subsetting genome
    if length(useChr) == 1
        global numChr = useChr[1] # useChr is considered to be # of chromosomes if an integer is supplied
        useChr = collect(1:useChr[1]) # select first n chromosomes
        starts = (header.+1)[useChr] # subset genome
        ends = [header[2:end] .- 1; size(genome)[1]][useChr]
        println("INFO: Using chromsome $(useChr) in the simulation!")
    elseif length(useChr) == length(header)
        global numChr = length(useChr) # useChr is considered to be a set of chromosomes if a vector is supplied
        starts = (header.+1)[useChr] # select chromosomes (specified by [useChr])
        ends = [header[2:end] .- 1; size(genome)[1]][useChr] # subset genome
        println("INFO: Using chromsome $(useChr) in the simulation!")
    else
        error("Wrong chromosome number/set specified, please check your parameters!")
    end
    for c = 1:numChr
        seq = [join(genome[starts[c]:ends[c]])] # read genomic sequence
        ## virtual digestion
        if size(re, 1) == 1
            frag = digest(seq, re[1]) # single-digestion
            reName = re[1].name
        else
            raw = digest(seq, re[1])
            frag = digestSecond(raw, re[2]) # double-digestion if more than one enzymes is supplied
            reName = [re[1].name re[2].name]
        end
        numFrag = length(frag) # number of fragments
        println("CHROMOSOME $c: $numFrag GBS fragments are generated through virtual digestion using $(reName).")
        len = length.(frag) # length of fragments
        pos = [0; cumsum(len)[1:(end-1)]] .+ 1 # starting position of each fragment
        chrLenEst = pos[end] / 1e6 # estimate chromosome length
        println("CHROMOSOME $c: Estimated chromosome length equals to $chrLenEst Mb.")
        ## (optional) subsetting chromosome
        if useChrLen == [] || size(useChrLen, 1) < numChr
            keep = [1:length(frag)...] # keep all fragments if subsetting is not required
            chrLen = [chrLen; chrLenEst] # store chromosome lengths
            println("CHROMOSOME $c: Subsetting chromosome is not required. Keeping all $numFrag available GBS fragments.")
        elseif size(useChrLen, 1) == numChr || useChrLen[c] <= estLen
            keep = findall(x -> x < useChrLen[c] * 1e6, pos) # keep only fragments found within specified length (varied at each chromosome)
            chrLen = useChrLen # store user-defined chromosome lengths
            println("CHROMOSOME $c (up to $(useChrLen[c]) cM): $(length(keep)) out of $numFrag available GBS fragments are kept.")
        elseif size(useChrLen, 1) == 1 # keep only fragments found within specified length (uniform across chromosomes)
            keep = findall(x -> x < useChrLen[1] * 1e6, pos)
            chrLen = [chrLen; useChrLen[1]]
            println("CHROMOSOME $c (up to $(useChrLen[c]) cM): $(length(keep)) out of $numFrag available GBS fragments are kept.")
        else
            error("CHROMOSOME $c: Wrong chromosome length(s) specified, please check your parameters!")
        end
        idFrag = [idFrag; [1:length(keep)...]] # identifier of each fragment
        chrFrag = [chrFrag; fill(c, length(keep))] # origin of chromosome per fragment
        posFrag = [posFrag; pos[keep]] # starting position per fragment
        lenFrag = [lenFrag; len[keep]] # length of each fragment
        seqFrag = [seqFrag; frag[keep]] # genomic sequence of each fragment
    end
    ## module 2: fragment size-selection
    fragRaw = [idFrag chrFrag posFrag lenFrag seqFrag] # dataset stores all fragments
    numFragTotal = length(seqFrag) # total number of GBS fragments
    selected = findall(x -> lowerThresh <= x <= upperThresh, lenFrag) # fragment size-selection
    numFragSelected = length(selected) # number of fragments selected in size-selection
    println("INFO: $numFragSelected out of $numFragTotal GBS fragments are selected after size-selection with lower and upper thresholds equals to [$lowerThresh, $upperThresh].")
    fragGBS = fragRaw[selected, :] # dataset conatins GBS fragments (after size-selection)
    global GBSFrag = (frag = fragGBS[:, 5], len = fragGBS[:, 4], pos = fragGBS[:, 3], chr = fragGBS[:, 2], id = fragGBS[:, 1])
    ## module 3: estimate coverage of selected GBS fragments (per Mb)
    GBSCoverChr = Array{Int64}(undef,0)
    GBSFragWinStatrs = Array{Int64}(undef,0)
    GBSFragWinEnds= Array{Int64}(undef,0)
    GBSFragCover =  Array{Float64}(undef,0)
    for c in 1:numChr
        numWin = Int(round(chrLen[c]))
        winStarts = Int.([0:(numWin-1)...]  .* 1e6)
        winEnds = Int.([1:numWin...] .* 1e6)
        GBSFragChr = findall(x -> x == c, GBSFrag.chr)
        GBSFragPos = GBSFrag.pos[fragChr]
        GBSFragLen = GBSFrag.len[fragChr]
        GBSFragWin = [findall(x -> x <= winEnds[w] && x >= winStarts[w], GBSFragPos) for w in 1:numWin]
        GBSFragCoverWin = [sum(GBSFragLen[GBSFragWin[w]])/1e6 for w in 1:numWin]
        GBSFragCover = [GBSFragCover; GBSFragCoverWin]
        GBSFragWinEnds = [GBSFragWinEnds; winEnds]
        GBSFragWinStatrs = [GBSFragWinStatrs; winStarts]
        GBSCoverChr = [GBSCoverChr; fill(c,numWin)]
    end
    GBSCover = [GBSCoverChr GBSFragWinStatrs GBSFragWinEnds GBSFragCover]
    println("INFO: Expected sequencing coverage based on $numFragSelected selected GBS fragments is approximately $(round(mean(GBSFragCover)*100;digits=2))%.")
    ## module 4: writing out and plotting fragment-szie distribution and GBS coverage if required
    if plotOutput == true
        histogram(fragRaw[:, 4], normalize = :probability, title = "Size Distribution of Raw GBS Fragements", label = "", fmt = :png); savefig("RawFragSize")
        histogram(fragGBS[:,4], normalize = :probability, title = "Size Distribution of GBS Fragements (after Size-Selection)", label = "", fmt = :png); savefig("GBSFragSize")
        histogram(GBSFragCover, normalize = :probability, title = "Coverage of GBS Fragments per Mb", label = "", fmt = :png); savefig("GBSCoverage")
    end
    if writeOutput == true
        writedlm("RawFrag.txt", fragRaw)
        writedlm("GBSFrag.txt", fragGBS)
        writedlm("GBSCoverage.txt", GBSCover)
    end
end;




## 02.Population
## define chromosome [recombination site, ancestor, SNP position, QTL position]
mutable struct chromosome
    Position::Array{Float64}
    Origin::Array{Int64}
    SNPs::Array{Int64}
    QTL::Array{Int64}
    chromosome(Position, Origin) = new(Position, Origin, [], [])
end;

## define individual [ID, marternal chromosome, paternal chromosome]
mutable struct individual
    ID::Int64
    MatChrs::Array{chromosome}
    PatChrs::Array{chromosome}
    individual(ID, MatChrs, PatChrs) = new(ID, MatChrs, PatChrs)
end;

## generate founders
generateFounders = function(numFounders::Int64)
    nChr = 2 * numFounders # each founder has two chromosomes
    f = Array{individual}(undef, numFounders)
    ## define founder chromosomes
    for i = 1:size(f, 1)
        f[i] = individual(i,[chromosome([0.0], [2 * i - 1]) for c = 1:numChr],[chromosome([0.0], [2 * i]) for c = 1:numChr],)
    end
    f
end;

## sample chromosomesfor each progeny
sampleChromosome = function (ind::individual)
    newChrs = Array{chromosome}(undef, numChr)
    for c = 1:numChr
        binomN = Int64(round(chrLen[c] * 3, digits = 0))
        recombs = sort(Array{Int64}(sample(1:round(1e8 * chrLen[c] / 100),rand(Binomial(binomN, chrLen[c] / 100 / binomN)),replace = false)))
        startOrigin = sample(1:2, 1)[1]
        if size(recombs, 1) == 0
            if startOrigin == 1
                newChrs[c] = ind.MatChrs[c]
            else
                newChrs[c] = ind.PatChrs[c]
            end
        else
            pos = [ind.MatChrs[c].Position, ind.PatChrs[c].Position]
            ori = [ind.MatChrs[c].Origin, ind.PatChrs[c].Origin]
            startPos = 0.0
            chrPos = Array{Float64}(undef, 0)
            chrOrig = Array{Int64}(undef, 0)
            for p in [recombs; 1e8 * chrLen[c] / 100]
                endPos = p
                a = [maximum([1:size(pos[startOrigin],1)...][findall(pos[startOrigin] .<= startPos)]):maximum([1:size(pos[startOrigin], 1)...][findall(pos[startOrigin] .< endPos)])...]
                tmpPos = pos[startOrigin][a]
                tmpPos[1] = startPos
                chrPos = [chrPos; tmpPos]
                chrOrig = [chrOrig; ori[startOrigin][a]]
                startPos = endPos
                startOrigin = (startOrigin - 3) * -1
            end
            newChrs[c] = chromosome(chrPos, chrOrig)
        end
    end
    newChrs
end;

## generate offSpring by inheriting one chromosome from each parent
sampleOffspring = function (sire::individual, dam::individual, id = indCount[1])
    sireChrs = sampleChromosome(sire)
    damChrs = sampleChromosome(dam)
    indCount[1] = indCount[1] + 1
    off = individual(id, damChrs, sireChrs)
    off
end;

## change population to a specified size (= $endSize) through a certain number (= $numGen) of geneartions
changingPopSize = function (founders::Array{individual,1}, endSize::Int64, numGen::Int64)
    startSize = size(founders, 1)
    popSize = Int.(round.(range(startSize, stop = endSize, length = numGen + 1)))[2:end]
    parents = copy(founders)
    for gen = 1:numGen
        offSpring = [sampleOffspring(parents[sample(1:(size(parents, 1)), 1)[1]], parents[sample(1:(size(parents, 1)), 1)[1]]) for i = 1:popSize[gen]]
        parents = deepcopy(offSpring)
        if gen % 10 == 0
            println("CHANING POP SIZE GEN $gen: Done!")
        end
    end
    parents
end;

## generate popluation at a fixed size (= $numGenFinal) for a certain number (=$numGen) of geneartions
constantPopSize = function (founders::Array{individual,1}, numGen::Int64, numGenFinal ::Int64, numIndFinal::Int64, useWeights::Array{Float64})
    parents = copy(founders)
    final = Array{individual}(undef, numIndFinal)
    for gen = 1:numGen-numGenFinal
        offSpring = [sampleOffspring(parents[sample(1:(size(parents, 1)), 1)[1]], parents[sample(1:(size(parents, 1)), 1)[1]]) for i = 1:size(founders, 1)]
        parents = deepcopy(offSpring)
        if gen % 10 == 0
            println("CONSTANT POP SIZE GEN $gen: Done")
        end
    end
    if useWeights != [] && size(useWeights,1) == numGenFinal
        useWeights = useWeights
    else
        useWeights = repeat([1 / numGenFinal], numGenFinal)
    end
    index = sample([1:numGenFinal...],Weights(useWeights),numIndFinal)
    count = [sum(index .== c) for c in 1:numGenFinal]
    for ind = 1:numGenFinal
        offSpring = [sampleOffspring(parents[sample(1:(size(parents, 1)), 1)[1]], parents[sample(1:(size(parents, 1)), 1)[1]]) for i = 1:size(founders, 1)]
        parents = deepcopy(offSpring)
        final[findall(x -> x == ind, index)] = parents[randperm(size(parents, 1))[1:count[ind]]]
        println("INFO: Collecting $(count[ind]) individual at Gen $(numGen -numGenFinal + ind)")
    end
    final
end;

## geneate a breeding population follows a specified (user-defined) pedigree
samplePedigree = function(pedFile::String, pedFounders::Array{individual,1}, output::Bool)
    println("INFO: Simulating from the given pedigree file: $pedFile.")
    ped = readdlm(pedFile, '\t', Int64, header = false)
    numInd = size(ped, 1)
    println("INFO: Sampling $numInd Individuals from the specified pedigree.")
    ind = Array{individual}(undef, numInd)
    for i = 1:numInd
        p = ped[i, :]
        id = p[1]
        sire = p[2]
        dam = p[3]
        if sire == 0
            sireInd = pedFounders[sample(1:(size(pedFounders, 1)), 1)[1]]
        else
            sireInd = ind[sire]
        end
        if dam == 0
            damInd = pedFounders[sample(1:(size(pedFounders, 1)), 1)[1]]
        else
            damInd = ind[dam]
        end
        ind[i] = sampleOffspring(sireInd, damInd, id)
        if output
            println("Individual $id: S = $(sireInd.ID); D = $(damInd.ID)")
        end
    end
    ind
end;

## define the population sturcture
definePopulation = function(numFounders::Int64, endSize::Int64, numGenCha::Int64, numGenCon::Int64, numGenFinal::Int64, numInd::Int64, useWeights::Array{Float64}, usePedigree::Bool, pedFile::String, pedOutput::Bool)
    global founders = generateFounders(numFounders) # generate founders
    global indCount = [numFounders + 1]
    inc = changingPopSize(founders, endSize, numGenCha); # increase population size over geneartions
    off = constantPopSize(inc, numGenCon, numGenFinal, numInd, useWeights); # create population over fixed number of generation and select individauls in the final few generations
    if usePedigree != false # simulate population structure follows a given pedigree - optional
        global ind = samplePedigree(pedFile, off, pedOutput)
    else
        global ind = off
    end
end;


## 03.Data
## randomly sample QTL positions across the genome
sampleQTLPosition = function(totalQTL::Int64)
    println("INFO: A total of $totalQTL QTLs are sampled randomly across $numChr chromosome(s).")
    qtlChr = sample([1:numChr...], totalQTL);
    global numQTL = [sum(qtlChr .== c) for c = 1:numChr]
    qtlPos = [sort(Int64.(sample([1:(chrLen[c]*1e6)...], numQTL[c]))) for c = 1:numChr]
    qtlPos
end;

## sampel SNP positions through a two-stage appraoch
sampleSNPPosition = function(totalSNP::Int64, winSize::Int64, mu::Float64, sigmasq::Float64)
    if totalSNP != 0 ## set totalSNP = 0 to trigger option 1
        println("INFO: A total of $totalSNP SNPs are sampled randomly across $numChr chromosome(s).")
        snpChr = sample([1:numChr...], totalSNP); # same method as sampling QTL positions
        global numSNP = [sum(snpChr .== c) for c = 1:numChr]
        ## option 1: sample SNP positions randomly across the genome
        snpPos = [sort(Int64.(sample([1:(chrLen[c]*1e6)...], numSNP[c]))) for c = 1:numChr]
    else
        ## option 2: sample SNP positions using densities
        global numSNP = Array{Int64}(undef,0)
        snpPos = Array{Int64}(undef,0)
        for c = 1:numChr
            numWin = Int(round(chrLen[c] * 1e6 / winSize)) # number of windows to be split in chromosome c (fixed window size, variable window number)
            theta = sigmasq/mu # calcualte Gamma paramter: scale
            alpha = mu/theta # calcualte Gamma paramter: shape
            density = rand(Gamma(alpha, theta), numWin) # sample SNP density for each window
            win = fill(winSize, numWin) # generate a vector (with length = numWin) contating windown size for each window
            winPos = [0; cumsum(win)[1:end-1]] # starting postion of each window
            sam = [rand(Bernoulli(density[i]), winSize) for i = 1:numWin] # sampling the occurence of SNP at each bp using local SNP density (sampled from Gamma model for each win)
            sam2 = map(x -> findall(x .> 0), sam) # record SNP positions within each window
            sam3 = reduce(vcat, [sam2[i] .+ winPos[i] for i = 1:numWin]) # calcualte SNP position within chromosome
            numSampled = length(sam3) # number of sampled SNP postions
            numSNP = [numSNP;numSampled] # record the number of sampled SNP postions for each chromosome
            snpPos = [snpPos;[sort(Int64.(sam3[randperm(numSampled)[1:numSampled]]))]] # sort SNP positions
            println("CHROMOSOME $c: $numSampled SNPs sampled with average SNP density = $(round(mean(Gamma(alpha,theta)),digits=3)) (window size = $winSize bp).")
        end
    end
    snpPos
end;

## sample SNP/QTL allele frequency
sampleAlleleFrequency = function(numLoci::Array{Int64}, mu::Float64,sigmasq::Float64)
    beta = ((mu - 1) * (mu^2 - mu + sigmasq)) / (sigmasq)
    alpha = (-beta * mu) / (mu - 1)
    af = Array{Float64}(undef,0)
    for c = 1:size(numLoci, 1)
        af = [af; [rand(Beta(alpha, beta), numLoci[c])]] # sample allele frequency from a Beta distribution
    end
    af
end;

## generate foudner SNPs
makeFounderSNPs = function(founders::Array{individual,1}, snpAF::Array{Any,1})
    for c = 1:size(snpAF, 1)
        numHaps = 2 * size(founders, 1)
        haps = [rand(Bernoulli(snpAF[c][i]),1)[1] for i = 1:size(snpAF[c], 1), j = 1:numHaps]
        for i = 1:size(founders, 1)
            founders[i].MatChrs[c].SNPs = haps[:, (2*i)-1]
            founders[i].PatChrs[c].SNPs = haps[:, 2*i]
        end
    end
    founders
end;

## generate foudner QTLs
makeFounderQTL = function(founders::Array{individual,1}, qtlAF::Array{Any,1})
    for c = 1:size(qtlAF, 1)
        numHaps = 2 * size(founders, 1)
        haps = [rand(Bernoulli(qtlAF[c][i]), 1)[1] for i = 1:size(qtlAF[c], 1), j = 1:numHaps]
        for i = 1:size(founders, 1)
            founders[i].MatChrs[c].QTL = haps[:, (2*i)-1]
            founders[i].PatChrs[c].QTL = haps[:, 2*i]
        end
    end
    founders
end;

## fill haplotypes
fillHaplotypes = function (samples::Array{individual,1}, founders::Array{individual,1}, numChr::Int64, snpPos, qtlPos::Array{Array{Int64,1},1})
    if size(founders[1].MatChrs[1].QTL, 1) == 0
        error("Haplotypes not made - need to run makeFounderSNPs(founders,snpAF)")
    else
        for c = 1:numChr
            println("CHROMOSOME $c: Filling haplotypes!")
            founderHaps = Array{Int64}(undef, size(snpPos[c], 1), 2 * size(founders, 1))
            for i = 1:size(founders, 1)
                founderHaps[:, 2*i-1] = founders[i].MatChrs[c].SNPs
                founderHaps[:, 2*i] = founders[i].PatChrs[c].SNPs
            end
            for an = 1:size(samples, 1)
                matOrig = [samples[an].MatChrs[c].Origin[maximum(findall(samples[an].MatChrs[c].Position .< snpPos[c][i],))] for i = 1:size(snpPos[c], 1)]
                patOrig = [samples[an].PatChrs[c].Origin[maximum(findall(samples[an].PatChrs[c].Position .< snpPos[c][i],))] for i = 1:size(snpPos[c], 1)]
                samples[an].MatChrs[c].SNPs = [founderHaps[i, matOrig[i]] for i = 1:size(snpPos[c], 1)]
                samples[an].PatChrs[c].SNPs = [founderHaps[i, patOrig[i]] for i = 1:size(snpPos[c], 1)]
            end
            if size(qtlPos, 1) > 0
                founderHaps = Array{Int64}(undef, size(qtlPos[c], 1), 2 * size(founders, 1))
                for i = 1:size(founders, 1)
                    founderHaps[:, 2*i-1] = founders[i].MatChrs[c].QTL
                    founderHaps[:, 2*i] = founders[i].PatChrs[c].QTL
                end
                for an = 1:size(samples, 1)
                    matOrig = [samples[an].MatChrs[c].Origin[maximum(findall(samples[an].MatChrs[c].Position .< qtlPos[c][i],))] for i = 1:size(qtlPos[c], 1)]
                    patOrig = [samples[an].PatChrs[c].Origin[maximum(findall(samples[an].PatChrs[c].Position .< qtlPos[c][i],))] for i = 1:size(qtlPos[c], 1)]
                    samples[an].MatChrs[c].QTL = [founderHaps[i, matOrig[i]] for i = 1:size(qtlPos[c], 1)]
                    samples[an].PatChrs[c].QTL = [founderHaps[i, patOrig[i]] for i = 1:size(qtlPos[c], 1)]
                end
            end
            println("CHROMOSOME $c: DONE!")
        end
    end
end;

## fucntion to extract QTL haplotypes can be included
## extract (x2) haplotypes from each (diploid) individual
getHaplotypes = function(samples::Array{individual,1} = ind)
    numInd = size(samples, 1)
    numChr = size(samples[1].MatChrs, 1)
    numSNP = [size(samples[1].MatChrs[i].SNPs, 1) for i = 1:numChr]
    totalSNP = cumsum(numSNP)
    starts = [2; totalSNP[1:(end-1)] .+ 2]
    ends = totalSNP .+ 1
    haplotypes = Array{Int64}(undef, 2*size(samples, 1), sum(numSNP) + 1)
    for i = 1:numInd
        haplotypes[2*i-1:2*i, 1] .= samples[i].ID
        for c = 1:numChr
            mat = samples[i].MatChrs[c].SNPs
            pat = samples[i].PatChrs[c].SNPs
            haplotypes[collect(1:2:2*numInd)[i], starts[c]:ends[c]] = mat
            haplotypes[collect(2:2:2*numInd)[i], starts[c]:ends[c]] = pat
        end
    end
    haplotypes
end;

## geneate SNP genotypes
getSNPGenotypes = function(samples::Array{individual,1} = ind)
    numChr = size(samples[1].MatChrs, 1)
    numSNP = [size(samples[1].MatChrs[i].SNPs, 1) for i = 1:numChr]
    totalSNP = cumsum(numSNP) # count total number of SNPs across the genome
    starts = [2; totalSNP[1:(end-1)] .+ 2]
    ends = totalSNP .+ 1
    genotypes = Array{Int64}(undef, size(samples, 1), sum(numSNP) + 1) # define empty array to store SNP genotypes
    for i = 1:size(samples, 1)
        genotypes[i, 1] = samples[i].ID
        for c = 1:numChr
            genos = samples[i].MatChrs[c].SNPs + samples[i].PatChrs[c].SNPs
            genotypes[i, starts[c]:ends[c]] = genos
        end
    end
    genotypes
end;

## geneate QTL genotypes
getQTLGenotypes = function(samples::Array{individual,1} = ind)
    numChr = size(samples[1].MatChrs, 1)
    numQTL = [size(samples[1].MatChrs[c].QTL, 1) for c = 1:numChr]
    totalQTL = cumsum(numQTL)
    starts = ifelse(size(totalQTL,1) == 1, 2, [2; totalQTL[1:(end-1)] .+ 2])
    ends = totalQTL .+ 1
    genotypes = Array{Int64}(undef, size(samples, 1), sum(numQTL) + 1)
    for i = 1:size(samples, 1)
        genotypes[i, 1] = samples[i].ID
        for c = 1:numChr
            genos = samples[i].MatChrs[c].QTL + samples[i].PatChrs[c].QTL
            genotypes[i, starts[c]:ends[c]] = genos
        end
    end
    genotypes
end;


## function to replicate values over array (credit:Milan Bouchet-Valat (https://github.com/JuliaLang/julia/issues/16443))
rep = function (x, lengths)
    if length(x) != length(lengths)
        throw(DimensionMismatch("vector lengths must match"))
    end
    res = similar(x, sum(lengths))
    i = 1
    for idx = 1:length(x)
        tmp = x[idx]
        for kdx = 1:lengths[idx]
            res[i] = tmp
            i += 1
        end
    end
    return res
end;


## (optional) generate key file
getKeyFile = function(barcode_file,numInd, numRow = 24, numCol = 16)
    barcode = readdlm(barcode_file, '\t', header = false);
    barcodes = barcode[1:numInd];
    lane = fill(1, numRow*numCol);
    row = reduce(vcat, rep(string.(['A':'P'...]), fill(numRow, numCol)));
    col = reduce(vcat, fill([1:numRow...],numCol));
    key = [fill("ABC12AAXX", numRow*numCol)[1:numInd] lane[1:numInd] barcode[1:numInd] [
    "Ind_" * string(i) for i = 1:(numInd)]  fill("Plate1", 384)[1:numInd] row[1:numInd] col[1:numInd]]
    writedlm("keyFile.txt", key)
    barcodes
end;



## sample GBS depth based on sequencing depth
sampleGBSDepth = function (numInd::Int64 = 64, sequecningDepth::Int64 = 10)
    # numLoci = length()
    # numSNP = length()
    cols = rand(Gamma(meandepth * 0.5, 2), numLoci)
    rows = rand(Gamma(meandepth * 3, 1 / 3), numInd)
    anis = reshape(anis, numInd, 1)
    counts = (anis * locus') ./ meandepth


    readDepth = [rand(NegativeBinomial(counts[i, j], 0.5), 1)[1] for i = 1:numInd, j = 1:numLoci]
    snpDepth =[rep(readDepth[i, :], collect(Iterators.flatten(hapSizeGBS))) for i = 1:numInd]
    snpDepth = transpose(reshape(vcat(snpDepth...), numSNP, numInd))
    ref = [0 for i = 1:numSNP, j = 1:numInd]
    for j = 1:numInd
        genotype = snpGenoGBS[j, 2:numSNP+1]
        copies = snpDepth[j, :]
        success = Array{Int64}(round.(100 * meandepth .* (genotype / 2)))
        failure = Array{Int64}(round.(100 * meandepth .* (abs.(2 .- genotype) / 2)))
        dist = [Hypergeometric(success[i], failure[i], copies[i]) for i = 1:size(genotype, 1)]
        ref[:, j] = [rand(dist[i], 1)[1] for i = 1:size(genotype, 1)]
        println("SAMPLE $j of $numInd: Done!")
    end
    ref = transpose(ref)
    alt = snpDepth - ref
    writedlm("readDepth.txt", readDepth)
    writedlm("snpDepth.txt", snpDepth)
    writedlm("ref.txt", ref)
    writedlm("alt.txt", alt)
end;

## simulate GBS
simGBS = function(totalQTL::Int64, totalSNP::Int64, muSNPdensity::Float64, sigmasqSNPdensity::Float64, winSize::Int64, muAlleleFreq::Float64,sigmasqAlleleFreq::Float64, re::Array{restrictionEnzyme,1}, seqDepth::Int64, barcodeFile::String, plotOutput = true, writeOutput = true)
    ## module 1. sample GBS variants
    ## step 1.1: sample QTL positions and allele frequencies
    qtlPos = sampleQTLPosition(totalQTL)
    numQTL = length.(qtlPos)
    qtlAF = sampleAlleleFrequency(numQTL, muAlleleFreq, sigmasqAlleleFreq)
    ##  step 1.2: sample SNP positions and allele frequencies
    snpPos = sampleSNPPosition(totalSNP,winSize,muDensity,sigmasqDensity)
    numSNP = length.(snpPos)
    snpAF = sampleAlleleFrequency(numSNP, muAlleleFreq, sigmasqAlleleFreq)
    ## step 2.1: generate founder QTL and SNPs
    makeFounderQTL(founders, qtlAF)
    makeFounderSNPs(founders, snpAF)
    ## step 2.2: fill haplotypes
    fillHaplotypes(ind, founders, numChr, snpPos, qtlPos)
    ## step 3: extract haplotypes, SNP and QTL genotypes
    haps = getHaplotypes(ind)
    snpGeno = getSNPGenotypes(ind)
    qtlGeno = getQTLGenotypes(ind)
    ## step 4: generate GBS data
    numInd = size(ind)[1] # number of individuals to be simulated
    # fragSNP = Array{Int64}(undef,0) # empty array to store GBS SNP positions
    snpGBS = Array{Int64}(undef,0) # empty array to store only SNP that found on GBS fragments
    hapSize = Array{Int64}(undef,0) # empty array to store haplotype size
    fragSNP = Array{Int64}(undef,0) # empty array to store GBS fragments (conataining SNPs) index
    readDepth = Array{Int64}(undef,numInd,0) # empty array to store readDepth
    io = open("ABC12AAXX_1_fastq.txt", "w")
    totalReads =  Array{String}(undef,0)
    ## Let's rock!!!
    for c = 1:numChr
        chr = findall(x -> x == c, GBSFrag.chr) # index fragments from chromosome c
        frags = GBSFrag.frag[chr] # select only fragments from chromosome c
        starts = GBSFrag.pos[chr] # starting position of selected GBS fragments
        len = GBSFrag.len[chr] # length of selected GBS fragments
        ends = starts + len[chr] .- 1 # ending position of selected GBS fragments
        numFrag = length(chr) # total number of selected GBS fragments
        snpSelected = [findall(x -> x == 1, (snpPos[c] .< ends[t]) .& (snpPos[c] .> starts[t])) for t = 1:numFrag] # select only SNPs that captured by each GBS fragment
        fragSelected = findall(x -> x != [], [snpPos[c][snpSelected[t]] for t = 1:numFrag]) # selected only non-empty (containing SNP(s)) GBS fragments
        numSNPperFrag = length.(snpSelected)[fragSelected] # number of SNPs within each non-empty GBS fragment
        snpGBS = [snpGBS; [(reduce(vcat,snpSelected[fragSelected]) .+ [0;numFrag][c])]] # index of SNPs only found on GBS fragments
        numSNPChr = sum(numSNPperFrag)
        numHapLoci = length(fragSelected) # total number of GBS fragments contating SNPs (or equivelently, number of shortHaps loci)
        fragSNP = [fragSNP; fragSelected .+ ifelse(c == 1,0,numFrag)] # save only GBS fragments containing SNPs
        println("CHROMOSOME $c: Found $numSNPChr SNPs on $numHapLoci GBS fragments, with an average of $(round(mean(numSNPperFrag);digits=2)) SNPs per GBS fragment.")
        ## sample depth
        cols = rand(Gamma(seqDepth * 0.5, 2), numHapLoci) # sample mean SNP depth
        rows = reshape(rand(Gamma(seqDepth * 3, 1 / 3), numInd),numInd,1)  # sample mean sample depth
        depth = (rows * cols') ./ seqDepth # genearte raw depth amtrix
        readDepthChr = [rand(NegativeBinomial(depth[i, j], 0.5)) for i = 1:numInd, j = 1:numHapLoci] # genearte read depth amtrix
        readDepth = [readDepth readDepthChr] # store read depth
        matHapCopyChr = [rand(Binomial(readDepthChr[i,j])) for i = 1:numInd, j = 1:numHapLoci] # sample number of copies of marternal haplotypes
        patHapCopyChr = readDepthChr - matHapCopyChr # compute number of copiesof parternal haplotypes
        barcodes = getKeyFile(barcode_file,numInd)
        for k in fragSelected
            frag = frags[k] # extract kth GBS fragments in chromsome c
            hapSizeChr = length.(snpSelected)[k] # number of SNPs found in this GBS fragment
            hapSize = [hapSize; hapSizeChr] # store haplotype size information
            sites = snpPos[c][snpSelected[k]] .- (starts[k].-1) # SNP sites on this GBS fragment
            startPos = [1; sites .+ 1] # starting position of segments of the specified GBS fragment
            endPos = [sites .- 1; len[k]] # ending position of segments of the specified GBS fragment
            seg = [SubString(frag, startPos[j], endPos[j]) for j = 1:length(startPos)] # poistional info of small segments of the specified GBS fragment
            bases = [split(frag[sites],"") split(map(x -> comp[x],frag[sites]),"")] # define SNPs (ref. allele = ref. allele found on the genome; alt. allele = complement of the ref. allele. E.g., A/T,C/G only)
            hap =  haps[:,snpSelected[k]].+1 # extract shortHaps, i.e. short haplotypes defined by GBS fragmentation
            if length(frag) < 101
                stops = length(frag) # read the entire fragment if its less than 101 bp
            else
                stops = 101 # if the length of fragment is longer than 101, stop at 101
            end
            matHap =  [join([[seg[j] * bases[j,hap[i,j]] for j = 1:hapSizeChr]; seg[hapSizeChr+1]; re[1].overhang[1]])[1:stops] for i in 1:2:numInd*2] # genearte GBS reads from martenal haplotypes
            patHap =  [join([[seg[j] * bases[j,hap[i,j]] for j = 1:hapSizeChr]; seg[hapSizeChr+1]; re[1].overhang[1]])[1:stops] for i in 2:2:numInd*2] # generate GBS reads from partenal haplotypes
            matHapRevComp = [reverse(map(x -> comp[x], matHap[i]))[1:stops] for i in 1: numInd] # reverse complement GBS reads of maternal haplotypes
            patHapRevComp = [reverse(map(x -> comp[x], patHap[i]))[1:stops] for i in 1: numInd] # reverse complement GBS reads of paternal haplotypes
            matRead = barcodes .* matHap # attach barcodes to the GBS reads
            patRead = barcodes .* patHap # attach barcodes to the GBS reads
            matReadRevComp = barcodes .* matHapRevComp # attach barcodes to the GBS reads
            patReadRevComp = barcodes .* patHapRevComp # attach barcodes to the GBS reads
            matCopyTotal = matHapCopyChr[:,  findall(x -> x == k, fragSelected)]
            matCopy = [rand(Binomial(matCopyTotal[i])) for i in 1:numInd]
            matCopyRevComp = matCopyTotal - matCopy
            patCopyTotal = patHapCopyChr[:,  findall(x -> x == k, fragSelected)]
            patCopy = [rand(Binomial(patCopyTotal[i])) for i in 1:numInd]
            patCopyRevComp = patCopyTotal - patCopy
            readGBS = collect(Iterators.flatten([[repeat([matRead[j]], matCopy[j]); repeat([matReadRevComp[j]], matCopyRevComp[j]); repeat([patRead[j]], patCopy[j]); repeat([patReadRevComp[j]], patCopyRevComp[j])] for j = 1:numInd]))
            totalReads = [totalReads; readGBS]
        end
    end
    totalReads = shuffle(totalReads[:, 1]);
    writedlm(io,["@SIM001:001:ABC12AAXX:1:0000:0000:$r 1:N:0:0\n" * totalReads[r] * "\n+\n" * repeat("I", length(totalReads[r])) for r in 1:length(totalReads)],quotes = false)
    close(io)
    # define QTL and SNP datasets [ID chromosome position allele_frequency]
    qtlData = [["QTL_"*string(i) for i in 1:totalQTL] reduce(vcat,[repeat([i],inner=numQTL[i]) for i in 1:numChr]) reduce(vcat, qtlPos) reduce(vcat, qtlAF)]
    snpData = [["SNP_"*string(i) for i in 1:sum(numSNP)] reduce(vcat,[repeat([i],inner=numSNP[i]) for i in 1:numChr]) reduce(vcat, snpPos) reduce(vcat, snpAF)]
    # define GBS SNP dataset, genotypes and haplotypes
    snpDataGBS = snpData[collect(Iterators.flatten(snpGBS)), :] # info about SNPs that captured by GBS fragments
    snpGenoGBS = snpGeno[:, ([1; (collect(Iterators.flatten(snpGBS)) .+ 1)])] # genotypes of SNPs that captured by GBS fragments
    hapIndex = [0;rep(fragSNP,hapSize)] # indexing shortHaps
    hapGBS = vcat(hapIndex',haps[:,([1; (collect(Iterators.flatten(snpGBS)) .+ 1)])]) # haplotypes that captured by GBS fragments
    if plotOutput == true
        histogram(snpAF,normalize = :probability, title = "SNP Allele Frequency", label = "",fmt = :png); savefig("snpAF")
        histogram(qtlAF,normalize = :probability, title = "QTL Allele Frequency", label = "",fmt = :png); savefig("qtlAF")
        histogram(snpDataGBS[:,4], normalize = :probability, title = "GBS SNP Allele Frequency", label = "",fmt = :png); savefig("snpAFGBS")
        histogram(hapSize, normalize = :probability, title = "Number of SNP per GBS Fragment",label = "",fmt = :png); savefig("snpPerTag")
    end
    if writeOutput == true
        writedlm("snpGeno.txt", snpGeno)
        writedlm("qtlGeno.txt",qtlGeno)
        writedlm("snpInfo.txt",snpData)
        writedlm("qtlInfo.txt",qtlData)
        writedlm("snpGenoGBS.txt", snpGenoGBS)
        writedlm("snpInfoGBS.txt", snpDataGBS)
        writedlm("shortHap.txt", hapGBS)
    end
 end;



end
