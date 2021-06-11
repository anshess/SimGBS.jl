## This file contains functions and types used for defining population structure
## define chromosome [recombination site, ancestor, SNP position, QTL position]
mutable struct chromosome
    Position::Array{Float64}
    Origin::Array{Int64}
    SNPs::Array{Int64}
    QTL::Array{Int64}
    chromosome(Position, Origin) = new(Position, Origin, [], [])
end

## define individual [ID, marternal chromosome, paternal chromosome]
mutable struct individual
    ID::Int64
    MatChrs::Array{chromosome}
    PatChrs::Array{chromosome}
    individual(ID, MatChrs, PatChrs) = new(ID, MatChrs, PatChrs)
end

## function to generate the founding population
function generateFounders(numFounders::Int64)
    nChr = 2 * numFounders # each founder has two chromosomes
    f = Array{individual}(undef, numFounders)
    ## define founder chromosomes
    for i = 1:size(f, 1)
        f[i] = individual(
            i,
            [chromosome([0.0], [2 * i - 1]) for c = 1:numChr],
            [chromosome([0.0], [2 * i]) for c = 1:numChr],
        )
    end
    f
end

## function to sample chromosomes of each offSpring
function sampleChromosome(ind)
    newChrs = Array{chromosome}(undef, numChr)
    for c = 1:numChr
        binomN = Int64(round(chrLen[c] * 3, digits = 0))
        recombs = sort(
            Array{Int64}(
                sample(
                    1:round(1e8 * chrLen[c] / 100),
                    rand(Binomial(binomN, chrLen[c] / 100 / binomN)),
                    replace = false,
                ),
            ),
        )
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
                a = [
                    maximum(
                        [1:size(pos[startOrigin], 1)...][findall(
                            pos[startOrigin] .<= startPos,
                        )],
                    ):maximum(
                        [1:size(pos[startOrigin], 1)...][findall(
                            pos[startOrigin] .< endPos,
                        )],
                    )...,
                ]
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
end

## function to generate offSpring
function sampleOffspring(sire, dam, id = indCount[1])
    sireChrs = sampleChromosome(sire)
    damChrs = sampleChromosome(dam)
    indCount[1] = indCount[1] + 1
    off = individual(id, damChrs, sireChrs)
    off
end

## function to generate a new population at a SPECIFIED size (= $endSize) through a certain number (= $numGen) of geneartions
function changingPopSize(founders, endSize::Int64, numGen::Int64)
    startSize = size(founders, 1)
    popSize = Int.(round.(range(startSize, stop = endSize, length = numGen + 1)))[2:end]
    parents = copy(founders)
    for gen = 1:numGen
        offSpring = [
            sampleOffspring(
                parents[sample(1:(size(parents, 1)), 1)[1]],
                parents[sample(1:(size(parents, 1)), 1)[1]],
            ) for i = 1:popSize[gen]
        ]
        parents = deepcopy(offSpring)
        if gen % 10 == 0
            println("CHANING POP SIZE GEN $gen: Done!")
        end
    end
    parents
end

## function to generate popluation at a FIXED size (= $numGenFinal) for a certain number (=$numGen) of geneartions
function constantPopSize(
    founders,
    numGen::Int64,
    numGenFinal::Int64,
    numIndFinal::Int64,
    useWeights::Array{Float64},
)
    parents = copy(founders)
    final = Array{individual}(undef, numIndFinal)
    for gen = 1:numGen-numGenFinal
        offSpring = [
            sampleOffspring(
                parents[sample(1:(size(parents, 1)), 1)[1]],
                parents[sample(1:(size(parents, 1)), 1)[1]],
            ) for i = 1:size(founders, 1)
        ]
        parents = deepcopy(offSpring)
        if gen % 10 == 0
            println("CONSTANT POP SIZE GEN $gen: Done")
        end
    end
    if useWeights != [] && size(useWeights, 1) == numGenFinal
        useWeights = useWeights
    else
        useWeights = repeat([1 / numGenFinal], numGenFinal)
    end
    index = sample([1:numGenFinal...], Weights(useWeights), numIndFinal)
    count = [sum(index .== c) for c = 1:numGenFinal]
    for ind = 1:numGenFinal
        offSpring = [
            sampleOffspring(
                parents[sample(1:(size(parents, 1)), 1)[1]],
                parents[sample(1:(size(parents, 1)), 1)[1]],
            ) for i = 1:size(founders, 1)
        ]
        parents = deepcopy(offSpring)
        final[findall(x -> x == ind, index)] =
            parents[randperm(size(parents, 1))[1:count[ind]]]
        println(
            "INFO: Collecting $(count[ind]) individual at Gen $(numGen -numGenFinal + ind)",
        )
    end
    final
end

## function to geneate a breeding population follows a user-defined pedigree
function samplePedigree(pedFile::String, pedFounders, output::Bool)
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
end

## function to set up the population sturcture
function definePopulation(
    numFounders::Int64,
    endSize::Int64,
    numGenCha::Int64,
    numGenCon::Int64,
    numGenFinal::Int64,
    numInd::Int64,
    useWeights::Array{Float64},
    usePedigree::Bool,
    pedFile::String,
    pedOutput::Bool,
)
    global founders = generateFounders(numFounders) # generate founders
    global indCount = [numFounders + 1]
    inc = changingPopSize(founders, endSize, numGenCha) # increase population size over geneartions
    off = constantPopSize(inc, numGenCon, numGenFinal, numInd, useWeights) # create population over fixed number of generation and select individauls in the final few generations
    if usePedigree != false # simulate population structure follows a given pedigree - optional
        global ind = samplePedigree(pedFile, off, pedOutput)
    else
        global ind = off
    end
end
