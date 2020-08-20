tests = ["reflecttests.jl", "datatests.jl"]

for test in tests
    include(test)
end