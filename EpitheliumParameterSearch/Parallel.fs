module Parallel

let threads = System.Environment.ProcessorCount *2

//Roughly 8% speedup in real time using the below approach relative to Array.Parallel

//BH02 (core i7 4650u, 8GB memory) benchmarks on kasumiCG

//Array.parallel
//let ans = CalculateProbability.parameterSearch Test.kasumiCG CalculateProbability.Simulation;;
//Real: 00:38:41.800, CPU: 01:51:21.171, GC gen0: 901896, gen1: 803919, gen2: 264

//Parallel.arrayMap
//let ans = CalculateProbability.parameterSearch Test.kasumiCG CalculateProbability.Simulation;;
//Real: 00:35:10.547, CPU: 01:45:06.156, GC gen0: 901920, gen1: 866762, gen2: 187

//Potential additional efficiency (makes most sense for big multicore machines:
//tasksPerThread = if tasks%threads <> 0 then tasks/threads+1 else tasks/threads
//(Avoids too many leftover tasks going on one thread)

let arrayMap f c =
    //Array.Parallel.map f c
    
    
    let tasks = Array.length c
    if tasks < threads then Array.Parallel.map f c else
    let tasksPerThread = tasks/threads
    let additionalTasks = tasks%threads
    let pstruct =   Array.init threads (fun i ->    let threadTasks = if i=(threads-1) then tasksPerThread + additionalTasks else tasksPerThread
                                                    Array.init threadTasks (fun j -> c.[i*tasksPerThread+j] )
                                                    )
                    |> Array.Parallel.map (fun taskBundle -> Array.map f taskBundle )
    //pstruct
    //Unbundle and return array
    Array.init tasks (fun taskID -> pstruct.[(min (taskID/tasksPerThread) (threads-1) )].[if (taskID/tasksPerThread) <= (threads-1) then (taskID%tasksPerThread) else taskID - tasksPerThread*(tasks/tasksPerThread-1)] )
    