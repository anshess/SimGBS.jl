## This file contains functions and types used for GBS data geneartion
## function yp sample QTL positions (randomly across the genome)
function sampleQTLPosition(totalQTL::Int64)
    println(
        "INFO: A total of $totalQTL QTLs are sampled randomly across $numChr chromosome(s).",
    )
    qtlChr = sample([1:numChr...], totalQTL)
    global numQTL = [sum(qtlChr .== c) for c = 1:numChr]
    qtlPos = [sort(Int64.(sample([1:(chrLen[c]*1e6)...], numQTL[c]))) for c = 1:numChr]
    qtlPos
end

## function to sample SNP positions (either randomly or via a two-stage appraoch)
function sampleSNPPosition(totalSNP::Int64, winSize::Int64, mu::Float64, sigmasq::Float64)
    if totalSNP != 0 ## set totalSNP = 0 to trigger option 1
        println(
            "INFO: A total of $totalSNP SNPs are sampled randomly across $numChr chromosome(s).",
        )
        snpChr = sample([1:numChr...], totalSNP) # same method as sampling QTL positions
        global numSNP = [sum(snpChr .== c) for c = 1:numChr]
        ## option 1: sample SNP positions randomly across the genome
        snpPos = [sort(Int64.(sample([1:(chrLen[c]*1e6)...], numSNP[c]))) for c = 1:numChr]
    else
        ## option 2: sample SNP positions using densities
        global numSNP = Array{Int64}(undef, 0)
        snpPos = Array{Int64}(undef, 0)
        for c = 1:numChr
            numWin = Int(round(chrLen[c] * 1e6 / winSize)) # number of windows to be split in chromosome c (fixed window size, variable window number)
            theta = sigmasq / mu # calcualte Gamma paramter: scale
            alpha = mu / theta # calcualte Gamma paramter: shape
            density = rand(Gamma(alpha, theta), numWin) # sample SNP density for each window
            win = fill(winSize, numWin) # generate a vector (with length = numWin) contating windown size for each window
            winPos = [0; cumsum(win)[1:end-1]] # starting postion of each window
            sam = [rand(Bernoulli(density[i]), winSize) for i = 1:numWin] # sampling the occurence of SNP at each bp using local SNP density (sampled from Gamma model for each win)
            sam2 = map(x -> findall(x .> 0), sam) # record SNP positions within each window
            sam3 = reduce(vcat, [sam2[i] .+ winPos[i] for i = 1:numWin]) # calcualte SNP position within chromosome
            numSampled = length(sam3) # number of sampled SNP postions
            numSNP = [numSNP; numSampled] # record the number of sampled SNP postions for each chromosome
            snpPos = [snpPos; [sort(Int64.(sam3[randperm(numSampled)[1:numSampled]]))]] # sort SNP positions
            println(
                "CHROMOSOME $c: $numSampled SNPs sampled with average SNP density = $(round(mean(Gamma(alpha,theta)),digits=3)) (window size = $winSize bp).",
            )
        end
    end
    snpPos
end

## function to sample variants (incl. SNP and QTL) allele frequency
function sampleAlleleFrequency(numLoci::Array{Int64}, mu::Float64, sigmasq::Float64)
    beta = ((mu - 1) * (mu^2 - mu + sigmasq)) / (sigmasq)
    alpha = (-beta * mu) / (mu - 1)
    af = Array{Float64}(undef, 0)
    for c = 1:size(numLoci, 1)
        af = [af; [rand(Beta(alpha, beta), numLoci[c])]] # sample allele frequency from a Beta distribution
    end
    af
end

## function to generate foudner SNPs
function makeFounderSNPs(founders, snpAF::Array{Any,1})
    for c = 1:size(snpAF, 1)
        numHaps = 2 * size(founders, 1)
        haps =
            [rand(Bernoulli(snpAF[c][i]), 1)[1] for i = 1:size(snpAF[c], 1), j = 1:numHaps]
        for i = 1:size(founders, 1)
            founders[i].MatChrs[c].SNPs = haps[:, (2*i)-1]
            founders[i].PatChrs[c].SNPs = haps[:, 2*i]
        end
    end
    founders
end

## function to generate foudner QTLs
function makeFounderQTL(founders, qtlAF::Array{Any,1})
    for c = 1:size(qtlAF, 1)
        numHaps = 2 * size(founders, 1)
        haps =
            [rand(Bernoulli(qtlAF[c][i]), 1)[1] for i = 1:size(qtlAF[c], 1), j = 1:numHaps]
        for i = 1:size(founders, 1)
            founders[i].MatChrs[c].QTL = haps[:, (2*i)-1]
            founders[i].PatChrs[c].QTL = haps[:, 2*i]
        end
    end
    founders
end

## function to fill haplotypes at individual level
function fillHaplotypes(
    samples,
    founders,
    numChr::Int64,
    snpPos,
    qtlPos::Array{Array{Int64,1},1},
)
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
                matOrig = [
                    samples[an].MatChrs[c].Origin[maximum(
                        findall(samples[an].MatChrs[c].Position .< snpPos[c][i]),
                    )] for i = 1:size(snpPos[c], 1)
                ]
                patOrig = [
                    samples[an].PatChrs[c].Origin[maximum(
                        findall(samples[an].PatChrs[c].Position .< snpPos[c][i]),
                    )] for i = 1:size(snpPos[c], 1)
                ]
                samples[an].MatChrs[c].SNPs =
                    [founderHaps[i, matOrig[i]] for i = 1:size(snpPos[c], 1)]
                samples[an].PatChrs[c].SNPs =
                    [founderHaps[i, patOrig[i]] for i = 1:size(snpPos[c], 1)]
            end
            if size(qtlPos, 1) > 0
                founderHaps = Array{Int64}(undef, size(qtlPos[c], 1), 2 * size(founders, 1))
                for i = 1:size(founders, 1)
                    founderHaps[:, 2*i-1] = founders[i].MatChrs[c].QTL
                    founderHaps[:, 2*i] = founders[i].PatChrs[c].QTL
                end
                for an = 1:size(samples, 1)
                    matOrig = [
                        samples[an].MatChrs[c].Origin[maximum(
                            findall(samples[an].MatChrs[c].Position .< qtlPos[c][i]),
                        )] for i = 1:size(qtlPos[c], 1)
                    ]
                    patOrig = [
                        samples[an].PatChrs[c].Origin[maximum(
                            findall(samples[an].PatChrs[c].Position .< qtlPos[c][i]),
                        )] for i = 1:size(qtlPos[c], 1)
                    ]
                    samples[an].MatChrs[c].QTL =
                        [founderHaps[i, matOrig[i]] for i = 1:size(qtlPos[c], 1)]
                    samples[an].PatChrs[c].QTL =
                        [founderHaps[i, patOrig[i]] for i = 1:size(qtlPos[c], 1)]
                end
            end
            println("CHROMOSOME $c: DONE!")
        end
    end
end


## function to extract haplotypes from each (diploid x2) individual
function getHaplotypes(samples = ind)
    numInd = size(samples, 1)
    numChr = size(samples[1].MatChrs, 1)
    numSNP = [size(samples[1].MatChrs[i].SNPs, 1) for i = 1:numChr]
    totalSNP = cumsum(numSNP)
    starts = [2; totalSNP[1:(end-1)] .+ 2]
    ends = totalSNP .+ 1
    haplotypes = Array{Int64}(undef, 2 * size(samples, 1), sum(numSNP) + 1)
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

## function to geneate SNP genotypes
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

## function to geneate QTL genotypes
function getQTLGenotypes(samples = ind)
    numChr = size(samples[1].MatChrs, 1)
    numQTL = [size(samples[1].MatChrs[c].QTL, 1) for c = 1:numChr]
    totalQTL = cumsum(numQTL)
    starts = ifelse(size(totalQTL, 1) == 1, 2, [2; totalQTL[1:(end-1)] .+ 2])
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
    return res
end

## function to generate key file for simulated GBS data
function getKeyFile(barcodeFile, numInd, numRow = 24, numCol = 16)
    barcode = readdlm(barcodeFile, '\t', header = false)
    barcodes = barcode[1:numInd]
    lane = fill(1, numRow * numCol)
    row = reduce(vcat, rep(string.(['A':'P'...]), fill(numRow, numCol)))
    col = reduce(vcat, fill([1:numRow...], numCol))
    key = [fill("ABC12AAXX", numRow * numCol)[1:numInd] lane[1:numInd] barcode[1:numInd] [
        "Ind_" * string(i) for i = 1:(numInd)
    ] fill("Plate1", 384)[1:numInd] row[1:numInd] col[1:numInd]]
    writedlm("keyFile.txt", key)
    barcodes
end


## function to simulate GBS
function GBS(
    totalQTL::Int64,
    totalSNP::Int64,
    muSNPdensity::Float64,
    sigmasqSNPdensity::Float64,
    winSize::Int64,
    muAlleleFreq::Float64,
    sigmasqAlleleFreq::Float64,
    re,
    seqDepth::Int64,
    barcodeFile::String,
    plotOutput::Bool,
    writeOutput::Bool,
)
    ## module 1. sample GBS variants
    ## step 1.1: sample QTL positions and allele frequencies
    qtlPos = sampleQTLPosition(totalQTL)
    numQTL = length.(qtlPos)
    qtlAF = sampleAlleleFrequency(numQTL, muAlleleFreq, sigmasqAlleleFreq)
    ##  step 1.2: sample SNP positions and allele frequencies
    snpPos = sampleSNPPosition(totalSNP, winSize, muDensity, sigmasqDensity)
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
    numInd = size(ind, 1) # number of individuals to be simulated
    barcodes = getKeyFile(barcodeFile, numInd) # extract barcodes for GBS reads generation
    snpGBS = Array{Int64}(undef, 0) # empty array to store only SNP that found on GBS fragments
    hapSize = Array{Int64}(undef, 0) # empty array to store haplotype size
    fragSNP = Array{Int64}(undef, 0) # empty array to store GBS fragments (conataining SNPs) index
    readDepth = Array{Int64}(undef, numInd, 0) # empty array to store read Depth
    snpDepth = Array{Int64}(undef, numInd, 0) # empty array to store snp Depth
    matDepth = Array{Int64}(undef, numInd, 0) # empty array to store copyies of maternal haplotypes
    patDepth = Array{Int64}(undef, numInd, 0) # empty array to store copyies of paternal haplotypes
    totalReads = Array{String}(undef, 0) # store simulated GBS reads
    io = GZip.open("ABC12AAXX_1_fastq.txt.gz", "w") # output file
    ## Let's party!!!
    for c = 1:numChr
        ## part a. extract fragments located in chromosome c
        chr = findall(x -> x == c, GBSFrag.chr) # index fragments from chromosome c
        numFrag = length(chr) # total number of selected GBS fragments
        starts = GBSFrag.pos[chr] # starting position of selected GBS fragments
        len = GBSFrag.len[chr] # length of selected GBS fragments
        ends = starts + len[chr] .- 1 # ending position of selected GBS fragments
        frags = GBSFrag.frag[chr] # select only fragments from chromosome c
        ## part b. select SNPs only SNPs that captured by each GBS fragment and non-empty (containing SNP(s)) GBS fragments
        snpSelected = [
            findall(x -> x == 1, (snpPos[c] .< ends[t]) .& (snpPos[c] .> starts[t])) for
            t = 1:numFrag
        ]
        fragSelected = findall(x -> x != [], [snpPos[c][snpSelected[t]] for t = 1:numFrag])
        hapSizeChr = length.(snpSelected)[fragSelected] # number of SNPs within each non-empty GBS fragment
        numSNPChr = sum(hapSizeChr) # total number of SNPs found on GBS fragments
        numHapLoci = length(fragSelected) # total number of GBS fragments contating SNPs (or equivelently, number of shortHaps loci)
        hapSize = [hapSize; hapSizeChr] # store haplotype size (number of SNPs per fragment)n
        snpGBS = [snpGBS; [(reduce(vcat, snpSelected[fragSelected]) .+ [0; numFrag][c])]] # store SNPs only found on GBS fragments
        fragSNP = [fragSNP; fragSelected .+ ifelse(c == 1, 0, numFrag)] # store only GBS fragments containing SNPs
        println(
            "CHROMOSOME $c: Found $numSNPChr SNPs on $numHapLoci GBS fragments, with an average of $(round(mean(hapSizeChr);digits=2)) SNPs per GBS fragment.",
        )
        ## part c. sample read depth based on user-defined sequencing depth
        cols = rand(Gamma(seqDepth * 0.5, 2), numHapLoci) # per locus depth
        rows = reshape(rand(Gamma(seqDepth * 3, 1 / 3), numInd), numInd, 1) # per sample depth
        depth = (rows * cols') ./ seqDepth # depth matrix
        readDepthChr =
            [rand(NegativeBinomial(depth[i, j], 0.5)) for i = 1:numInd, j = 1:numHapLoci] # read depth
        snpDepthChr = transpose(
            reshape(
                vcat(
                    [
                        rep(readDepthChr[i, :], collect(Iterators.flatten(hapSizeChr)))
                        for i = 1:numInd
                    ]...,
                ),
                sum(hapSizeChr),
                numInd,
            ),
        )# snp depth
        matHapCopyChr =
            [rand(Binomial(readDepthChr[i, j])) for i = 1:numInd, j = 1:numHapLoci] # sample number of copies of marternal haplotypes
        patHapCopyChr = readDepthChr - matHapCopyChr # compute number of copiesof parternal haplotypes
        readDepth = [readDepth readDepthChr] # store read depth
        snpDepth = [snpDepth snpDepthChr] # store snp depth
        matDepth = [matDepth matHapCopyChr] # store number of copies of each maternal haplotypes
        patDepth = [patDepth patHapCopyChr]  # store number of copies of each paternal haplotypes
        ## part d. genearte GBS reads
        for k in fragSelected
            frag = frags[k] # extract kth fragment in chromsome c
            sites = snpPos[c][snpSelected[k]] .- (starts[k] .- 1) # SNP sites on this GBS fragment
            startPos = [1; sites .+ 1] # starting position of segments of the specified GBS fragment
            endPos = [sites .- 1; len[k]] # ending position of segments of the specified GBS fragment
            seg = [SubString(frag, startPos[j], endPos[j]) for j = 1:length(startPos)] # poistional info of small segments of the specified GBS fragment
            bases = [split(frag[sites], "") split(map(x -> comp[x], frag[sites]), "")] # define SNPs (ref. allele = ref. allele found on the genome; alt. allele = complement of the ref. allele. E.g., A/T,C/G only)
            hap = haps[:, snpSelected[k]] .+ 1 # extract shortHaps, i.e. short haplotypes defined by GBS fragmentation
            if length(frag) < 101
                stops = length(frag) # read the entire fragment if its less than 101 bp
            else
                stops = 101 # if the length of fragment is longer than 101, stop at 101
            end
            matHap = [
                join(
                    [
                        [seg[j] * bases[j, hap[i, j]] for j = 1:length(sites)]
                        seg[length(sites)+1]
                        re[1].overhang[1]
                    ],
                )[1:stops] for i = 1:2:numInd*2
            ] # genearte GBS reads from martenal haplotypes
            patHap = [
                join(
                    [
                        [seg[j] * bases[j, hap[i, j]] for j = 1:length(sites)]
                        seg[length(sites)+1]
                        re[1].overhang[1]
                    ],
                )[1:stops] for i = 2:2:numInd*2
            ] # generate GBS reads from partenal haplotypes
            matHapRevComp =
                [reverse(map(x -> comp[x], matHap[i]))[1:stops] for i = 1:numInd] # reverse complement GBS reads of maternal haplotypes
            patHapRevComp =
                [reverse(map(x -> comp[x], patHap[i]))[1:stops] for i = 1:numInd] # reverse complement GBS reads of paternal haplotypes
            matRead = barcodes .* matHap # attach barcodes to the GBS reads
            patRead = barcodes .* patHap # attach barcodes to the GBS reads
            matReadRevComp = barcodes .* matHapRevComp # attach barcodes to the GBS reads
            patReadRevComp = barcodes .* patHapRevComp # attach barcodes to the GBS reads
            ## replicate GBS reads
            matCopyTotal = matHapCopyChr[:, findall(x -> x == k, fragSelected)]
            matCopy = [rand(Binomial(matCopyTotal[i])) for i = 1:numInd]
            matCopyRevComp = matCopyTotal - matCopy
            patCopyTotal = patHapCopyChr[:, findall(x -> x == k, fragSelected)]
            patCopy = [rand(Binomial(patCopyTotal[i])) for i = 1:numInd]
            patCopyRevComp = patCopyTotal - patCopy
            readGBS = collect(
                Iterators.flatten([
                    [
                        repeat([matRead[j]], matCopy[j])
                        repeat([matReadRevComp[j]], matCopyRevComp[j])
                        repeat([patRead[j]], patCopy[j])
                        repeat([patReadRevComp[j]], patCopyRevComp[j])
                    ] for j = 1:numInd
                ]),
            )
            totalReads = [totalReads; shuffle(readGBS)]
        end
    end
    #    totalReads = totalReads[:, 1])
    writedlm(
        io,
        [
            "@SIM001:001:ABC12AAXX:1:0000:0000:$r 1:N:0:0\n" *
            totalReads[r] *
            "\n+\n" *
            repeat("I", length(totalReads[r])) for r = 1:length(totalReads)
        ],
        quotes = false,
    )
    close(io)
    ## step4. returen GBS data info
    ## define QTL and SNP datasets [ID chromosome position allele_frequency]
    qtlData =
        [["QTL_" * string(i) for i = 1:totalQTL] reduce(vcat, [repeat([i], inner = numQTL[i]) for i = 1:numChr]) reduce(vcat, qtlPos) reduce(vcat, qtlAF)]
    snpData =
        [["SNP_" * string(i) for i = 1:sum(numSNP)] reduce(vcat, [repeat([i], inner = numSNP[i]) for i = 1:numChr]) reduce(vcat, snpPos) reduce(vcat, snpAF)]
    ## define GBS SNP dataset, genotypes and haplotypes
    snpDataGBS = snpData[collect(Iterators.flatten(snpGBS)), :] # info about SNPs that captured by GBS fragments
    snpGenoGBS = snpGeno[:, ([1; (collect(Iterators.flatten(snpGBS)) .+ 1)])] # genotypes of SNPs that captured by GBS fragments
    hapIndex = [0; rep(fragSNP, hapSize)] # indexing shortHaps
    hapGBS = vcat(hapIndex', haps[:, ([1; (collect(Iterators.flatten(snpGBS)) .+ 1)])]) # haplotypes that captured by GBS fragments
    ## determine the number of ref. and alt. alleles for GBS SNPs
    ref = [0 for i = 1:sum(hapSize), j = 1:numInd]
    for j = 1:numInd
        genotype = snpGenoGBS[j, 2:end]
        copies = snpDepth[j, :]
        success = Array{Int64}(round.(100 * seqDepth .* (genotype / 2)))
        failure = Array{Int64}(round.(100 * seqDepth .* (abs.(2 .- genotype) / 2)))
        dist =
            [Hypergeometric(success[i], failure[i], copies[i]) for i = 1:size(genotype, 1)]
        ref[:, j] = [rand(dist[i], 1)[1] for i = 1:size(genotype, 1)]
        # println("SAMPLE $j of $numInd: Done!")
    end
    ref = transpose(ref)
    alt = snpDepth - ref
    if plotOutput == true
        histogram(
            snpAF,
            normalize = :probability,
            title = "SNP Allele Frequency",
            label = "",
            fmt = :png,
        )
        savefig("snpAF")
        histogram(
            mean(snpGeno, dims = 1),
            normalize = :probability,
            title = "SNP Allele Frequency (Genotypes)",
            label = "",
            fmt = :png,
        )
        savefig("snpGenoAF")
        # plot(snpAF,mean(snpGeno,dims=1), title = "Sampled Allele Frequency vs. Simulated Allele Frequency (Genotypes)", label = "",fmt = :png); savefig("snpAFvsGenoAF")
        histogram(
            qtlAF,
            normalize = :probability,
            title = "QTL Allele Frequency",
            label = "",
            fmt = :png,
        )
        savefig("qtlAF")
        histogram(
            snpDataGBS[:, 4],
            normalize = :probability,
            title = "GBS SNP Allele Frequency",
            label = "",
            fmt = :png,
        )
        savefig("snpAFGBS")
        histogram(
            hapSize,
            normalize = :probability,
            title = "Number of SNP per GBS Fragment",
            label = "",
            fmt = :png,
        )
        savefig("snpPerTag")
        histogram(
            mean(readDepth, dims = 1),
            normalize = :probability,
            title = "Mean Read Depth",
            label = "",
            fmt = :png,
        )
        savefig("readDepthGBS")
        histogram(
            reduce(
                vcat,
                [count(i -> (i >= 0), snpDepth[:, j]) for j = 1:size(readDepth, 2)] / numInd,
            ),
            normalize = :probability,
            title = "Call Rates",
            label = "",
            fmt = :png,
        )
        savefig("callRates")
    end
    if writeOutput == true
        writedlm("snpGeno.txt", snpGeno)
        writedlm("qtlGeno.txt", qtlGeno)
        writedlm("snpInfo.txt", snpData)
        writedlm("qtlInfo.txt", qtlData)
        writedlm("snpGenoGBS.txt", snpGenoGBS)
        writedlm("snpInfoGBS.txt", snpDataGBS)
        writedlm("shortHap.txt", hapGBS)
        writedlm("readDepth.txt", readDepth)
        writedlm("snpDepth.txt", snpDepth)
        writedlm("ref.txt", ref)
        writedlm("alt.txt", alt)
    end
end
