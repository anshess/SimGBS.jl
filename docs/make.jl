push!(LOAD_PATH,"../src/")

using Documenter, SimGBS

makedocs(
    sitename = "SimGBS.jl"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
