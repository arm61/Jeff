using Documenter, Jeff, Measurements

makedocs(sitename="Jeff", 
         modules=[Jeff],
         pages = [
             "API" => [
                "reflect" => "reflect.md",
                "data" => "data.md"
             ]
         ],
         repo = "github.com/arm61/Jeff.jl.git",)