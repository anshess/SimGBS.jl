using Documenter, SimGBS


makedocs(modules = [SimGBS], sitename = "SimGBS.jl", pages = ["Home" => "index.md"])

deploydocs(repo = "github.com/kanji709/SimGBS.jl.git", devbranch = "master")
