using Distributed

Threads.@threads for i in 1:Threads.nthreads()
    println("i: ", i, "\t Thread ID: ", Threads.threadid())
end


