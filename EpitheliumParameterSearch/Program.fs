// Learn more about F# at http://fsharp.net
// See the 'F# Tutorial' project for more help.

type Job = Simulation

type programParam   =   {   clone   :   SimulationCloneSizeDistribution.clone
                            task    :   Job
                            output  :   string
                        }

let defaultParam    =   {clone=SimulationCloneSizeDistribution.specificClone;task=Simulation;output="data.dat"}

let parseJob t =
    match t with
    | "Simulation" -> Simulation
    | _ -> failwith "Task not recognised"

//Will use the defaults without warning if not explicitly handed to the program. 
//This should either give warnings or errors unless explicitly supressed
let rec parseArguments a p = 
    match a with
    | [] -> p
    | "-rho" :: rho :: remainder -> parseArguments remainder {p with clone={p.clone with rho=float(rho)}}
    | "-lambda" :: lambda :: remainder -> parseArguments remainder {p with clone={p.clone with lambda=(float(lambda)*1.<Types.cell/Types.week>)}}
    | "-r" :: r :: remainder -> parseArguments remainder {p with clone={p.clone with r=float(r)}}
    | "-task" :: t :: remainder -> parseArguments remainder {p with task=(parseJob t)}

let runSimulation c oup =
    //Runs a simulation of a single clone, gives a text file for an output
    SimulationCloneSizeDistribution.simulate c 
    |> List.map (fun (i:SimulationCloneSizeDistribution.populationState)->sprintf "%f %f %f" (float(i.time)) (float(i.population.basal)) (float(i.population.suprabasal))) //process
    |> fun a -> IO.quickOutput oup [a]

[<EntryPoint>]
let main argv = 
    printfn "%A" argv
    let p = parseArguments (List.ofArray argv) defaultParam
    match p.task with
    | Simulation -> runSimulation p.clone p.output
    0 // return an integer exit code
