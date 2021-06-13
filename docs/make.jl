using Documenter, SimGBS


makedocs(
 modules = [SimGBS], 
 sitename = "SimGBS.jl",
 # pages = ["Home" => "index.md"]
)


# deploy_config = "GitHubActions"

deploydocs(
 repo = "github.com/kanji709/SimGBS.jl.git"
)
