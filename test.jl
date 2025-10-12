println("Hello world")

function f(x, y)
    return x-y*y
end
x = 1
y = 2
println("x is $x and y is $y. f is $(f(x,y))")
println("f(x,y) is ",f(x,y))

A = [1 2 3; 4 3 5]
B = [2 1; 2 1; 2 3]
C = A*A'
println(C)
try
    println(inv(C))
    println(inv([0 1; 0 1]))
catch e
    println("Error: ", e)
end
