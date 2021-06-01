using Documenter, SimGBS

makedocs(sitename = "SimGBS.jl", modules = [SimGBS], pages = ["Home" => "index.md"])

deploydocs(repo = "github.com/kanji709/SimGBS.jl.git")
