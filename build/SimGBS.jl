module SimGBS

using Compat
using DelimitedFiles         # for reading and writing text-delimited files
using GZip                   # for reading and writing gzipped files
using Random                 # for random permutation
using Distributions          # for sampling distribution
using StatsBase              # for weighted vectors
using Plots                   # for plotting

export digestGenome  # main function for in slico digestion using restriction enzyme
export definePopulation       # main function for defining population structure
export GBS                    # main function for GBS simulation

include("genome.jl")          # step 1. in slico digestion using restriction enzyme
include("pop.jl")             # step 2. set up population structure
include("gbs.jl")             # step 3. simultae GBS data

end
