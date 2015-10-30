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
    //Unbundle and return array
    Array.init tasks (fun taskID -> pstruct.[(min (tasks/taskID) (threads-1) )].[tasks%taskID] )
    