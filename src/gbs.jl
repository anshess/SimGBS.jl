## This file contains functions and types used for GBS data geneartion
"function yp sample QTL positions (randomly across the genome)"
function sampleQTLPosition(totalQTL::Int64)
    println("INFO: A total of $totalQTL QTLs are sampled randomly across $numChr chromosome(s).")
    qtlChr = sample([1:numChr...], totalQTL);
    global numQTL = [sum(qtlChr .== c) for c = 1:numChr]
    qtlPos = [sort(Int64.(sample([1:(chrLen[c]*1e6)...], numQTL[c]))) for c = 1:numChr]
    qtlPos
end

"function to sample SNP positions (either randomly or via a two-stage appraoch)"
function sampleSNPPosition(totalSNP::Int64, winSize::Int64, mu::Float64, sigmasq::Float64)
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
            # theta = sigmasq/mu # calcualte Gamma paramter: scale
            # alpha = mu^2/sigmasq # calcualte Gamma paramter: shape
            density = exp.(rand(Laplace(mu, sigmasq), numWin)) # sample SNP density for each window
            win = fill(winSize, numWin) # generate a vector (with length = numWin) contating windown size for each window
            winPos = [0; cumsum(win)[1:end-1]] # starting postion of each window
            sam = [rand(Bernoulli(1.4*density[i]), winSize) for i = 1:numWin] # sampling the occurence of SNP at each bp using local SNP density (sampled from Gamma model for each win)
            sam2 = map(x -> findall(x .> 0), sam) # record SNP positions within each window
            sam3 = reduce(vcat, [sam2[i] .+ winPos[i] for i = 1:numWin]) # calcualte SNP position within chromosome
            numSampled = length(sam3) # number of sampled SNP postions
            numSNP = [numSNP;numSampled] # record the number of sampled SNP postions for each chromosome
            snpPos = [snpPos;[sort(Int64.(sam3[randperm(numSampled)[1:numSampled]]))]] # sort SNP positions
            println("CHROMOSOME $c: $numSampled SNPs sampled with average SNP density = $(round(mean(density),digits=3)) (window size = $winSize bp).")
        end
    end
    snpPos
end

"function to sample variants (incl. SNP and QTL) allele frequency"
function sampleAlleleFrequency(numLoci::Array{Int64}, mu::Float64,sigmasq::Float64)
    alpha = 0.953 # (mu * (1 - mu) / sigmasq -1) * mu # (-beta * mu) / (mu - 1)
    beta = 3.631 #(mu * (1 - mu) / sigmasq -1) * (1 - mu) # ((mu - 1) * (mu^2 - mu + sigmasq)) / (sigmasq)
    af = Array{Float64}(undef,0)
    for c = 1:size(numLoci, 1)
        af = [af; [rand(Beta(alpha, beta), numLoci[c])]] # sample allele frequency from a Beta distribution
    end
    af
end

"function to generate foudner SNPs"
function makeFounderSNPs(founders, snpAF::Array{Any,1})
    for c = 1:size(snpAF, 1)
        numHaps = 2 * size(founders, 1)
        haps = [rand(Bernoulli(snpAF[c][i]),1)[1] for i = 1:size(snpAF[c], 1), j = 1:numHaps]
        for i = 1:size(founders, 1)
            founders[i].MatChrs[c].SNPs = haps[:, (2*i)-1]
            founders[i].PatChrs[c].SNPs = haps[:, 2*i]
        end
    end
    founders
end

"function to generate foudner QTLs"
function makeFounderQTL(founders, qtlAF::Array{Any,1})
    for c = 1:size(qtlAF, 1)
        numHaps = 2 * size(founders, 1)
        haps = [rand(Bernoulli(qtlAF[c][i]), 1)[1] for i = 1:size(qtlAF[c], 1), j = 1:numHaps]
        for i = 1:size(founders, 1)
            founders[i].MatChrs[c].QTL = haps[:, (2*i)-1]
            founders[i].PatChrs[c].QTL = haps[:, 2*i]
        end
    end
    founders
end

"function to fill haplotypes at individual level"
function fillHaplotypes(samples, founders, numChr::Int64, useChr::Array{Int64}, snpPos, qtlPos::Array{Array{Int64,1},1})
    if size(founders[1].MatChrs[1].QTL, 1) == 0
        error("Haplotypes not made - need to run makeFounderSNPs(founders,snpAF)")
    else
        for c = 1:numChr
            println("CHROMOSOME $(useChr[c]): Filling haplotypes!")
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
            println("CHROMOSOME $(useChr[c]): DONE!")
        end
    end
end


"function to extract haplotypes from each (diploid x2) individual"
function getHaplotypes(samples = ind)
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
end

"function to geneate SNP genotypes"
function getSNPGenotypes(samples = ind)
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
end

"function to geneate QTL genotypes"
function getQTLGenotypes(samples = ind)
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
end


"function to replicate values over array - Milan Bouchet-Valat (https://github.com/JuliaLang/julia/issues/16443)"
function rep(x, lengths)
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
    res
end

"function to sample read depth"
function sampleReadDepth(numLoci::Int64, numInd::Int64, meanDepth::Float64)
# function sampleReadDepth(numLoci::Int64, numInd::Int64, meanDepth::Float64, sigmasqReadDepth::Float64, sigmasqSampleDepth::Float64, meanCallRate::Float64, sigmasqReadCallRate::Float64, sigmasqSampleCallRate::Float64)
    readDepth = rand(Gamma(meanDepth/11,11), numLoci)
    sampleDepth = reshape(rand(Gamma(meanDepth*3.1,1/3), numInd),numInd,1)
    depthRaw = (sampleDepth * readDepth') ./ meanDepth
    depthMat = [rand(NegativeBinomial(depthRaw[i, j], 0.5)) for i = 1:numInd, j = 1:numLoci]
    # Caution: changes are needed here!
    dp = log10.(Base.dropdims(mean(depthMat, dims = 1); dims = 1))
    dp[isinf.(dp)] .= -1
    # sampleCallRate = reshape(rand(Laplace(meanCallRate, sqrt(sigmasqSampleCallRate / 2)), numInd), numInd, 1)
    # readCallRate = rand(Beta((meanCallRate * (1- meanCallRate) / sigmasqReadCallRate - 1) * meanCallRate, (meanCallRate * (1- meanCallRate) / sigmasqReadCallRate - 1) * (1 - meanCallRate)), numLoci)
    cr = Base.dropdims(mean(depthMat .!= 0, dims = 1); dims = 1)
    @. model(x, p) = p[1] / (1 + exp(-(x - p[2]) * p[3]))
    p0 = [1, 0.5, 5]
    fit_cr = curve_fit(model, dp, cr, p0)
    xm = rand(Exponential(0.3), numLoci) # rand(Uniform(-1,1), numLoci)
    for i = 1:numLoci
        par1 = fit_cr.param[1]
        par2 = fit_cr.param[2]
        par3 = fit_cr.param[3]
        output = par1 / (1 + exp(-(dp[i] - par2 - xm[i]) * par3))
        if output > 1
            output = 1
        end
        if output > cr[i]
            output = cr[i]
        end
        zeros = findall(x -> x == 0, depthMat[:, i]) #depthMat[i, sample([1:numLoci...],Int(round(numLoci*rand(Laplace(1-callRate,0.sigmasqCallRate),1))), replace=false)] .= 0
        newzeros = Int64(floor((1 - output) * numInd))
        potentials = deleteat!([1:numInd...], zeros)
        diff = newzeros - length(zeros)
        if diff > 0
            setzero = sample(potentials, diff, replace = false)
            s = sum(depthMat[setzero, i])
            allzeros = sort([setzero; zeros...])
            k = deleteat!([1:numInd...], allzeros)
            l = sample(k, s, replace = true)
            depthMat[setzero,i] .= 0
            for j in l
                depthMat[j, i] = depthMat[j, i] + 1
            end
        end
    end
    depthMat
end


"function to generate key file for simulated GBS data"
function getKeyFile(barcodeFile, numInd, flowcell, lane, numRow = 24, numCol = 16)
    barcode = readdlm(barcodeFile, '\t', header = false)
    useBarcode = barcode[1:numInd]
    row = reduce(vcat, rep(string.(['A':'P'...]), fill(numRow, numCol)))
    col = reduce(vcat, fill([1:numRow...],numCol))
    key = [fill("$flowcell", numRow*numCol)[1:numInd] fill(lane, numRow*numCol)[1:numInd] useBarcode ["Ind_" * string(i) for i = 1:(numInd)]  fill("Plate1", numInd) row[1:numInd] col[1:numInd]]
    writedlm("keyFile_$(flowcell)_$(lane).txt", key)
    useBarcode
end

"""
	GBS(totalQTL, totalSNP, muDensity, sigmasqDensity, winSize, muAlleleFreq, sigmasqAlleleFreq, re, meanDepth, barcodeFile, useChr, plotOutput, writeOutput, onlyOutputGBS)

Simulate Genotyping-by-Sequencing (GBS) data.

This function generates GBS reads by inserting genomic variants into _in silico_ digested genomic fragments, ligates the polymorphic sequence with barcodes and replicates based on sequencing depth.

# Arguments
* `totalQTL`: total number of QTL to be simulated
* `totalSNP`: total number of SNPs to be simulated (set to "0" if sampling SNP positions based on density)
* `muDensity`: location parameter of log-Laplace distribution (for sampling SNP density)
* `sigmasqDensity`: scale parameter of log-Laplace distribution (for sampling SNP density)
* `winSize`: Size of window and bin for sampling SNP positions
* `muAlleleFreq`: mean of sampled allele frequency
* `sigmasqAlleleFreq`: variance of sampled allele frequency
* `re`: restriction enzyme(s) to be used
* `barcodeFile`: file containing GBS barcodes
* `useChr`: either the number of chromosomes or a set of chromosome(s) to be simulated
* `plotOutput`: set to true if graphical outputs are required
* `writeOutput`: set to true if text outputs are required
* `onlyOutputGBS`: set to true if only GBS data is kept


# Examples
```julia
julia> GBS(100, 0, -2.0, 0.001, 1000000, 0.5, 0.001, [SimGBS.ApeKI], 20.0, "GBS_Barcodes.txt", [1], false, true, true)
```
"""
function GBS(totalQTL::Int64, totalSNP::Int64, muSNPdensity::Float64, sigmasqSNPdensity::Float64, winSize::Int64, muAlleleFreq::Float64,sigmasqAlleleFreq::Float64, re, meanDepth::Float64, barcodeFile::String,useChr::Array{Int64}, plotOutput::Bool, writeOutput::Bool,onlyOutputGBS::Bool)
# function GBS(totalQTL::Int64, totalSNP::Int64, muSNPdensity::Float64, sigmasqSNPdensity::Float64, winSize::Int64, muAlleleFreq::Float64,sigmasqAlleleFreq::Float64, re, meanDepth::Float64, sigmasqReadDepth::Float64, sigmasqSampleDepth::Float64, meanCallRate::Float64, sigmasqReadCallRate::Float64, sigmasqSampleCallRate::Float64, barcodeFile::String, plotOutput::Bool, writeOutput::Bool,outputOnlyGBS::Bool)
     ## 1. sample variants
     ### 1.1 QTL positions, number and allele frequencies
     qtlPos = sampleQTLPosition(totalQTL)
     numQTL = length.(qtlPos)
     qtlAF = sampleAlleleFrequency(numQTL, muAlleleFreq, sigmasqAlleleFreq)
     ###  1.2 SNP positions
     snpPos = sampleSNPPosition(totalSNP,winSize,muSNPdensity,sigmasqSNPdensity)
     numSNP = length.(snpPos)
     numFrag = [length(findall(x -> x == useChr[c], GBSFrag.chr)) for c in 1:numChr] # number of GBS fragments found within each chromosome
     ### 1.3 SNP captured by GBS fragments
     snpGBS = Array{Int64}(undef,0) # empty array to store only SNP found on GBS fragments
     hapSize = Array{Int64}(undef,0) # empty array to store haplotype size (number of SNPs within each hap)
     fragSNP = Array{Int64}(undef,0) # empty array to store GBS fragments (conataining SNPs) index
     for c = 1:numChr
         #### 1.3.1 extract GBS fragments from chromosome c
         chr = findall(x -> x == useChr[c], GBSFrag.chr) # index fragments from chromosome c
         numFragChr = numFrag[c] # total number of selected GBS fragments
         starts = GBSFrag.pos[chr] # starting position of selected GBS fragments
         len = GBSFrag.len[chr] # length of selected GBS fragments
         ends = starts + len .- 1 # ending position of selected GBS fragments
         #### 1.3.2 select only SNPs captured by GBS fragments and non-empty GBS fragments (contain SNPs)
         snpSelected = [findall(x -> x == 1, (snpPos[c] .< ends[t]) .& (snpPos[c] .> starts[t])) for t = 1:numFragChr]
         fragSelected = findall(x -> x != [], snpSelected) # findall(x -> x != [], [snpPos[c][snpSelected[t]] for t = 1:numFragChr])
         hapSizeChr = length.(snpSelected)[fragSelected] # size of shortHaps (i.e., number of SNPs within each)
         numSNPChr = sum(hapSizeChr) # total number of SNPs found on GBS fragments
         numHapLoci = length(fragSelected) # total number of GBS fragments contating SNPs (or equivelently, number of shortHaps loci)
         #### 1.3.3 store SNP index, shortHap size and fragment index
         snpGBS = [snpGBS; [reduce(vcat,snpSelected[fragSelected])]] # store index of SNPs found on GBS fragments
         # [snpGBS; [(reduce(vcat,snpSelected[fragSelected]) .+ [0;length.(snpPos)][c])]] # store index of SNPs found on GBS fragments
         hapSize = [hapSize; hapSizeChr] # store shortHap size (number of SNPs with each shortHap)
         fragSNP = [fragSNP; fragSelected .+ [0;numFrag][c]] # store only fragments contain SNPs
         println("CHROMOSOME $(useChr[c]): Found $numSNPChr SNPs on $numHapLoci GBS fragments, with an average of $(round(mean(hapSizeChr),digits=2)) SNPs per GBS fragment.")
     end
     #### 1.3.4 (Optional) keep only SNPs captured by GBS fragments for the subsequent steps
     if onlyOutputGBS == true
         println("[INFO]: a total of $(sum(numSNP)) SNPs were sampled.")
         numSNP = length.(snpGBS) # update number of SNPs
         println("[INFO]: $(sum(numSNP)) SNPs captured by selected GBS SNPs.")
         snpPos = [snpPos[c][snpGBS[c]] for c in 1:numChr]
     end
     ### 1.4 SNP allele frequencies
     snpAF = sampleAlleleFrequency(numSNP, muAlleleFreq, sigmasqAlleleFreq)
     ### 1.5 founder QTL and SNPs
     makeFounderQTL(founders, qtlAF)
     makeFounderSNPs(founders, snpAF)
     ### 1.6 fill haplotypes
     @time fillHaplotypes(ind, founders, numChr, useChr, snpPos, qtlPos)
     ### 1.7 extract haplotypes
     haps = getHaplotypes(ind)
     ## 2. generate GBS reads
     ### 2.1 sample read depth
     totalInd = size(ind,1) # number of samples to be simulated
     totalFrag = length(fragSNP) # number of loci
     totalDepth = Array{Int64}(undef,totalInd,totalFrag)
     ### 2.2 inster variants
     #### 2.2.1 subset haplotype alleel by chromosome
     hapStart = [1; numSNP[1:end-1] .+ 1] .+ 1
     hapEnd = cumsum(numSNP) .+1
     #### 2.2.2 genearte FASTQ file(s)
     flowcell = "ABC12AAXX" # flowcell name
     multiplex = 384 # number of samples per lane
     numLane = Int(ceil(totalInd/multiplex)) # numebr of lanes in total
     for lane in 1:numLane
         if totalInd < multiplex*lane
             numInd = totalInd - multiplex*(lane-1) # number of samples to be 'sequenced'
         else
             numInd = multiplex # number of samples to be 'sequenced'
         end
         readDepth = sampleReadDepth(totalFrag, numInd, meanDepth) # depth
         useInd = [i + (lane - 1) * multiplex for i in 1:numInd] # which individual gets sequenced (ascending order)
         totalDepth[useInd,:] = readDepth
         println("[INFO] On Lane $lane of Flowcell $flowcell: Average GBS read depth equals to $(round(mean(readDepth),digits=2)).")
         matHapDepth = [rand(Binomial(readDepth[i,j])) for i in 1:numInd, j in 1:totalFrag]
         patHapDepth = readDepth - matHapDepth # depth of parternal haplotypes
         matDepthFwd = [rand(Binomial(matHapDepth[i,j])) for i in 1:numInd, j in 1:totalFrag]
         matDepthRevComp = matHapDepth - matDepthFwd
         patDepthFwd = [rand(Binomial(patHapDepth[i,j])) for i in 1:numInd,j in 1:totalFrag]
         patDepthRevComp = patHapDepth - patDepthFwd
         println("[INFO] On Lane $lane of Flowcell $flowcell: Generating GBS data for $numInd samples.")
         barcodes = getKeyFile(barcodeFile, numInd, flowcell, lane) # extract barcodes
         io = GZip.open("$(flowcell)_$(lane)_fastq.txt.gz", "w") # output file
         numReadsTotal = 0 # reads counter
         for c = 1:numChr
             ##### 2.2.2.1 extract GBS fragments and SNPs
             chr = findall(x -> x == useChr[c], GBSFrag.chr) # index fragments from chromosome c
             numFragChr = numFrag[c] # total number of selected GBS fragments
             starts = GBSFrag.pos[chr] # starting position of selected GBS fragments
             len = GBSFrag.len[chr] # length of selected GBS fragments
             ends = starts + len .- 1 # ending position of selected GBS fragments
             frags = GBSFrag.frag[chr] # select only fragments from chromosome c
             ##### 2.2.2.2 select only SNPs that captured by each GBS fragment and non-empty (containing SNP(s)) GBS fragments
             snpSelected = [findall(x -> x == 1, (snpPos[c] .< ends[t]) .& (snpPos[c] .> starts[t])) for t = 1:numFragChr]
             fragSelected = findall(x -> x != [], snpSelected) # findall(x -> x != [], [snpPos[c][snpSelected[t]] for t = 1:numFragChr])
             hapSizeChr = length.(snpSelected[fragSelected]) # number of SNPs within each non-empty GBS fragment
             numSNPChr = sum(hapSizeChr) # total number of SNPs found on GBS fragments
             numHapLoci = length(fragSelected) # total number of GBS fragments contating SNPs (or equivelently, number of shortHaps loci)
             ##### 2.2.2.3 genearte GBS reads with variants
             hapChr = haps[:,hapStart[c]:hapEnd[c]]
             for f in 1:length(fragSelected)
                 k = fragSelected[f]
                 frag = frags[k] # extract kth fragment in chromsome c
                 sites = snpPos[c][snpSelected[k]] .- starts[k] # (starts[k].-1) # SNP sites on this GBS fragment
                 startPos = [1; sites .+ 1] # starting position of segments of the specified GBS fragment
                 endPos = [sites .- 1; len[k]] # ending position of segments of the specified GBS fragment
                 seg = [SubString(frag, startPos[j], endPos[j]) for j = 1:length(startPos)] # poistional info of small segments of the specified GBS fragment
				 bases = [split(map(x -> nucl[x],frag[sites]),"") split(map(x -> comp[nucl[x]],frag[sites]),"")] # define SNPs (ref. allele = ref. allele found on the genome; alt. allele = complement of the ref. allele. E.g., A/T,C/G only)
                 hap = hapChr[:,snpSelected[k]] .+1 # extract shortHaps (i.e. short haplotypes defined by GBS fragmentation), adding 1 here to covert
                 matHap = hap[1:2:(totalInd*2),:]
                 patHap = hap[2:2:(totalInd*2),:]
                 if length(frag) < 101
                     stops = length(frag) # read the entire fragment if its less than 101 bp
                 else
                     stops = 101 # if the length of fragment is longer than 101, stop at 101
                 end
                 matHapChrTemp = [join([[seg[j] * bases[j, matHap[i,j]] for j = 1:length(sites)]; seg[length(sites)+1]; re[1].overhang[1]]) for i in useInd]
                 patHapChrTemp =  [join([[seg[j] * bases[j, patHap[i,j]] for j = 1:length(sites)]; seg[length(sites)+1]; re[1].overhang[1]]) for i in useInd]
                 matHapChr = [matHapChrTemp[i][1:stops] for i in 1: numInd]
                 patHapChr = [patHapChrTemp[i][1:stops] for i in 1: numInd]
                 matHapRevCompChr = [reverse(map(x -> comp[nucl[x]], matHapChrTemp[i]))[1:stops] for i in 1: numInd] # reverse complement GBS reads of maternal haplotypes
                 patHapRevCompChr = [reverse(map(x -> comp[nucl[x]], patHapChrTemp[i]))[1:stops] for i in 1: numInd] # reverse complement GBS reads of paternal haplotypes
                 readGBS = collect(Iterators.flatten([[rep([barcodes[i] .* matHapChrTemp[i]], matDepthFwd[i,f]); rep([barcodes[i] .* matHapRevCompChr[i]], matDepthRevComp[i,f]); rep([barcodes[i] .* patHapChr[i]], patDepthFwd[i,f]); rep([barcodes[i] .* patHapRevCompChr[i]], patDepthRevComp[i,f])] for i in 1:numInd]))
                 totalReads = shuffle(readGBS)
                 numReads = length(totalReads)
                 numReadsTotal = numReadsTotal + numReads
                 writedlm(io,["@SIM001:001:ABC12AAXX:$lane:0000:0000:$r 1:N:0:0\n" * totalReads[r] * "\n+\n" * repeat("I", length(totalReads[r])) for r in 1:length(totalReads)], quotes = false)
             end
         end
         println("INFO: A total of $numReadsTotal GBS reads generated.")
         close(io)
     end
     ## 3. GBS data
     #### 3.1.1 SNP and QTL genotypes
     snpGeno = getSNPGenotypes(ind)
     qtlGeno = getQTLGenotypes(ind)
     ### 3.2 SNP and QTL datasets [ID chromosome position allele_frequency]
     qtlData = [["QTL_"*string(i) for i in 1:totalQTL] reduce(vcat,[repeat([useChr[i]],inner=numQTL[i]) for i in 1:numChr]) reduce(vcat, qtlPos) reduce(vcat, qtlAF)]
     snpData = [["SNP_"*string(i) for i in 1:sum(numSNP)] reduce(vcat,[repeat([useChr[i]],inner=numSNP[i]) for i in 1:numChr]) reduce(vcat, snpPos) reduce(vcat, snpAF)]
     ### 3.3 GBS fragments conatin SNPs
     snpFragGBS = [GBSFrag.id[collect(Iterators.flatten(fragSNP))] GBSFrag.chr[collect(Iterators.flatten(fragSNP))] GBSFrag.pos[collect(Iterators.flatten(fragSNP))] GBSFrag.len[collect(Iterators.flatten(fragSNP))] GBSFrag.frag[collect(Iterators.flatten(fragSNP))]]
     ### 3.4 , shortHaps
     hapIndex = [0;rep(fragSNP,hapSize)] # indexing shortHaps
     hapGBS = vcat(hapIndex',haps)
     ### 3.5 SNP depth, call rates
     snpDepth = transpose(reshape(vcat([rep(totalDepth[i, :], collect(Iterators.flatten(hapSize))) for i = 1:totalInd]...), sum(hapSize), totalInd))
     if writeOutput == true
         writedlm("snpGeno.txt", snpGeno)
         writedlm("qtlGeno.txt",qtlGeno)
         writedlm("snpInfo.txt",snpData)
         writedlm("qtlInfo.txt",qtlData)
         writedlm("snpFragGBS.txt",snpFragGBS)
         writedlm("shortHap.txt", hapGBS)
         writedlm("readDepth.txt", totalDepth)
         writedlm("snpDepth.txt", snpDepth)
         # writedlm("allele.txt",allele[:,2:end])
     end
     if plotOutput == true
         histogram(snpAF,normalize = :probability, title = s"SNP Allele Frequency", bins =1000, label = "",fmt = :png); savefig("snpAF")
         histogram(qtlAF,normalize = :probability, title = "QTL Allele Frequency", bins =1000, label = "",fmt = :png); savefig("qtlAF")
         histogram(hapSize, normalize = :probability, title = "Number of SNP per GBS Fragment", bins =1000, label = "",fmt = :png); savefig("snpPerTag")
     end
  end
