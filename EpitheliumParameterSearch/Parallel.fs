module Parallel

let threads = System.Environment.ProcessorCount *2

let arrayMap f c =
    //Array.Parallel.map f c
    let tasks = Array.length c
    let tasksPerThread = tasks/threads
    let additionalTasks = tasks%threads
    let pstruct =   Array.init threads (fun i ->    let threadTasks = if i=(threads-1) then tasksPerThread + additionalTasks else tasksPerThread
                                                    Array.init threadTasks (fun j -> c.[i*tasksPerThread+j] )
                                                    )
                    |> Array.Parallel.map (fun taskBundle -> Array.map f taskBundle )
    //pstruct
    //Unbundle and return array
    Array.init tasks (fun taskID -> pstruct.[(min (taskID/tasksPerThread) (threads-1) )].[if (taskID/tasksPerThread) <= (threads-1) then (taskID%tasksPerThread) else taskID - tasksPerThread*(tasks/tasksPerThread-1)] )
    