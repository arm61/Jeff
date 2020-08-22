import Pkg
Pkg.add("Documenter")

using Documenter, Jeff, Distributions

makedocs(sitename="Jeff.jl", 
         modules=[Jeff],
         pages = [
             "Home" => "index.md",
             "API" => [
                "reflect" => "reflect.md",
                "data" => "data.md",
                "model" => "model.md",
                "objective" => "objective.md",
             ]
         ])

deploydocs(repo = "github.com/arm61/Jeff.jl.git")