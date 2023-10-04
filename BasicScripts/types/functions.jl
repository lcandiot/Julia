function myFunc(a, b)
    println("This is a function")
    a + b
end

result = myFunc(5,3)
println(result)

function foo(a, b, z = 10)
    return (a+b)*z
end
println(foo(2,3,4))
println(foo(2,3))

function summit(args...)
    sum = 0
    for a in args
        sum += a 
    end

    return sum
end
println(summit(1,2))
