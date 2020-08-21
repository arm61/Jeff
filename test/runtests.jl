tests = ["reflecttests.jl", "datatests.jl", "modeltests.jl"]

for test in tests
    include(test)
end