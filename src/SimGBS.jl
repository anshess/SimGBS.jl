## Main code for running SimGBS
module SimGBS

export restrictionEnzyme,
    digestGenome,           # main function for in slico digestion using restriction enzyme
    definePopulation,       # main function for defining population structure
    GBS                    # main function for GBS simulation


using Compat,
    DelimitedFiles,         # for reading and writing text-delimited files
    GZip,                   # for reading and writing gzipped files
    Random,                 # for random permutation
    Distributions,          # for sampling distribution
    StatsBase,              # for weighted vectors
    Plots,                  # for plotting
    LsqFit  # fitting non-linear models


include("genome.jl")          # step 1. in slico digestion using restriction enzyme
include("pop.jl")             # step 2. set up population structure
include("gbs.jl")             # step 3. simultae GBS data

end
