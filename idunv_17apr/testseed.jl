
using Random

function testseed()

    Random.seed!(1234)


    println(randn())
end

