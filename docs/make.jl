using Documenter
using SimGBS

makedocs(
    sitename = "SimGBS",
    format = Documenter.HTML(),
    modules = [SimGBS]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
