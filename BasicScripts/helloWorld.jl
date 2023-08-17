# My first Julia script on the MacBook Pro 14''
function helloWorld()
    println("Hello World!")
end

# Call the function
helloWorld()

# Writing another hello world function which will call the function as many
# times as the scalar input
function helloWorld(count::Int)
    for i in 1:count
        println("Hello World!", i)
    end
end

helloWorld(3)