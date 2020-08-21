tests = ["reflecttests.jl", "datatests.jl", "modeltests.jl", "objectivetests.jl"]

for test in tests
    include(test)
end