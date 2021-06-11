"restriction enzyme used for digesting genome [name, cut-sequence, cut-site, overhang]"
mutable struct restrictionEnzyme
    name::String # name of the restriction enzyme
    cutSeq::Array{String} # recognition sequence
    cutPlace::Array{Int64} # where it gets cut
    overhang::Array{String} # overhang of the restriction enzyme if exists
    restrictionEnzyme(name, cutSeq, cutPlace, overhang) = new(name, cutSeq, cutPlace, overhang)
end

## list of restraiction enzymes
ApeKI = restrictionEnzyme("ApeKI", ["GCAGC", "GCTGC"], [2, 2], ["CTG", "CAG"]); # Elshire et al. 2011
PstI = restrictionEnzyme("PstI", ["CTGCAG"], [2], ["TGCA"]); # Poland et al. 2012
MspI = restrictionEnzyme("MspI", ["CCGG"], [2], ["CGG"]); # Poland et al. 2012
PstI_MspI = restrictionEnzyme("PstI-MspI", ["CTGCAG","CCGG"], [2,2], ["TGCA", "CGG"]); # Poland et al. 2012


## nucleotides and its complements
"nucleotides"
nucl = Dict('A' => 'A', 'C' => 'C','G' => 'G','T' => 'T','N' => 'N','R' => rand(['A', 'G']), 'Y' => rand(['C', 'T']), 'S' => rand(['G', 'C']), 'W' => rand(['A', 'T']),
    'K' => rand(['G', 'T']), 'M' => rand(['A', 'C']), 'B' => rand(['C', 'G', 'T']), 'D' => rand(['A', 'G', 'T']),
    'H' => rand(['A', 'C', 'T']), 'V' => rand(['A', 'C', 'G'])); # randomly assign IUPAC codes for uncertain ones in the assembly

"complements"
comp = Dict('A' => 'T', 'C' => 'G','G' => 'C','T' => 'A','N' => 'N'); # complements of DNA nucleotides



"split sequence by any specific sequence"
function splitSequence(seq, cutSeq, startSeq, endSeq)
    cuts = split(uppercase(seq), cutSeq) # split DNA sequence
    if length(cuts) > 1
        cuts[2:length(cuts)] = startSeq .* cuts[2:length(cuts)] # re-attach the strat seq to the fragment after cut-place
        cuts[1:(length(cuts)-1)] = cuts[1:(length(cuts)-1)] .* endSeq # re-attach the end seq to the fragment before cut-place
    else
        cuts
    end
    cuts
end

"digest sequence using specified restriction enzyme"
function digest(sequence, re::restrictionEnzyme)
    cutSeq = re.cutSeq
    cutPlace = re.cutPlace
    # overhang = re.overhang
    startSeq = [SubString(cutSeq[i], cutPlace[i]) for i = 1:length(cutSeq)]
    endSeq =  [SubString(cutSeq[i], 1, cutPlace[i] - 1) for i = 1:length(cutSeq)] # [SubString(cutSeq[i], 1, cutPlace[i] - 1) * overhang[i] for i = 1:length(cutSeq)]
    for k = 1:length(cutSeq)
        cuts = sequence
        sequence = []
        for s in cuts
            sequence = [sequence; splitSequence(s, cutSeq[k], startSeq[k], endSeq[k])]
        end
        sequence
    end
    # if length(sequence) > 2
    #    sequence[2:length(sequence)]
    # else
    # end
    sequence
end

"double-digestion (or enzyme has more than one recongnition site)"
function digestSecond(sequence, re::restrictionEnzyme)
    cutSeq = re.cutSeq
    cutPlace = re.cutPlace
    # overhang = re.overhang
    startSeq = [SubString(cutSeq[i], cutPlace[i]) for i = 1:length(cutSeq)]
    endSeq =  [SubString(cutSeq[i], 1, cutPlace[i] - 1) for i = 1:length(cutSeq)] #[SubString(cutSeq[i], 1, cutPlace[i] - 1) .* overhang[i] for i = 1:length(cutSeq)]
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
end

"""
	digestGenome(genofile, re, useChr, useChrLen, lower ,upper, plotOutput, writeOutput)

Perform virtual digestion using restriction enzyme and generate GBS fragments.

This function uses specified restriction enzyme(s) to digest genome and therefore generate GBS fragments. Fragment size-selection step is also included.

# Arguments
- `genofile`: file containing the reference genome
- `re`: restriction enzyme(s) to be used
- `useChr`: either the number of chromosome or a set of chromosome(s) to be simulated
- `useChrLen`: length of chromosome in cM to be used in simulation, otherwise using entire chromosome
- `lower`: lower threshold of fragment size-selection
- `upper`: upper threshold of fragment size-selection
- `winSize`: size of window used for calculating average genomic coverage
- `plotOutput`: set to true if graphical outputs are required
- `writeOutput`: set to true if text outputs are required
...

# Examples
```julia
julia> digestGenome("ref.fa.gz", [SimGBS.ApeKI], [1], Array{Float64}(undef,0), 65 ,195, 1000000, false, true)
```
"""
function digestGenome(genofile::String, re::Array{restrictionEnzyme,1}, useChr::Array{Int64}, useChrLen::Array{Float64}, lowerThresh::Int64, upperThresh::Int64, winSize::Int64, plotOutput::Bool , writeOutput::Bool)
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
    if length(useChr) <= length(header)
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
        println("CHROMOSOME $(useChr[c]): $numFrag GBS fragments are generated through virtual digestion using $(reName).")
        len = length.(frag) # length of fragments
        pos = [0; 1 .+ cumsum(len)[1:(end-1)]] # starting position of each fragment
        chrLenEst = pos[end] / 1e6 # estimate chromosome length
        println("CHROMOSOME $(useChr[c]): Estimated chromosome length equals to $chrLenEst Mb.")
        ## (optional) subsetting chromosome
        if useChrLen == [] || size(useChrLen, 1) < numChr
            keep = [1:length(frag)...] # keep all fragments if subsetting is not required
            chrLen = [chrLen; chrLenEst] # store chromosome lengths
            println("CHROMOSOME $(useChr[c]): Subsetting chromosome is not required. Keeping all $numFrag available GBS fragments.")
        elseif size(useChrLen, 1) == numChr || useChrLen[c] <= estLen
            keep = findall(x -> x < useChrLen[c] * 1e6, pos) # keep only fragments found within specified length (varied at each chromosome)
            chrLen = useChrLen # store user-defined chromosome lengths
            println("CHROMOSOME $(useChr[c]) (up to $(useChrLen[c]) cM): $(length(keep)) out of $numFrag available GBS fragments are kept.")
        elseif size(useChrLen, 1) == 1 # keep only fragments found within specified length (uniform across chromosomes)
            keep = findall(x -> x < useChrLen[1] * 1e6, pos)
            chrLen = [chrLen; useChrLen[1]]
            println("CHROMOSOME $(useChr[c]) (up to $(useChrLen[c]) cM): $(length(keep)) out of $numFrag available GBS fragments are kept.")
        else
            error("CHROMOSOME $(useChr[c]): Wrong chromosome length(s) specified, please check your parameters!")
        end
        idFrag = [idFrag; [1:length(keep)...]] # identifier of each fragment
        chrFrag = [chrFrag; fill(useChr[c], length(keep))] # origin of chromosome per fragment
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
        numWin = Int(round(chrLen[c]) .* 1e6 / winSize)
        winStarts = Int.([0:(numWin-1)...]  .* winSize)
        winEnds = Int.([1:numWin...] .* winSize)
        GBSFragChr = findall(x -> x == useChr[c], GBSFrag.chr)
        GBSFragPos = GBSFrag.pos[GBSFragChr]
        GBSFragLen = GBSFrag.len[GBSFragChr]
        GBSFragCoverWin = Array{Float64}(undef,numWin)
        for w in 1:numWin
            GBSFragWin = findall(x -> x <= winEnds[w] && x >= winStarts[w], GBSFragPos)
            if size(GBSFragWin,1) != 0
                GBSFragCoverWin[w] = sum(GBSFragLen[GBSFragWin])/winSize
            else
                GBSFragCoverWin[w] = 0
            end
        end
        GBSFragCover = [GBSFragCover; GBSFragCoverWin]
        GBSFragWinEnds = [GBSFragWinEnds; winEnds]
        GBSFragWinStatrs = [GBSFragWinStatrs; winStarts]
        GBSCoverChr = [GBSCoverChr; fill(c,numWin)]
    end
    GBSCover = [GBSCoverChr GBSFragWinStatrs GBSFragWinEnds GBSFragCover]
    println("INFO: Expected sequencing coverage based on $numFragSelected selected GBS fragments is approximately $(round(mean(GBSFragCover)*100;digits=2))%.")
    ## module 4: writing out and plotting fragment-szie distribution and GBS coverage if required
    if plotOutput == true
        histogram(fragRaw[:, 4], normalize = :probability, title = "Size Distribution of Raw GBS Fragments", label = "", fmt = :png); savefig("RawFragSize")
        histogram(fragGBS[:,4], normalize = :probability, title = "Size Distribution of GBS Fragments (after Size-Selection)", label = "", fmt = :png); savefig("GBSFragSize")
        histogram(GBSFragCover, normalize = :probability, title = "Coverage of GBS Fragments per Mb", label = "", fmt = :png); savefig("GBSCoverage")
    end
    if writeOutput == true
        writedlm("RawFrag.txt", fragRaw)
        writedlm("GBSFrag.txt", fragGBS)
        writedlm("GBSCoverage.txt", GBSCover)
    end
end
