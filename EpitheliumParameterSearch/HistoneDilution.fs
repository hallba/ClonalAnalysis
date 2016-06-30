module HistoneDilution

open SimulationCloneSizeDistribution

let rec addCells (n,a) pop =
    if n = 0<Types.cell> then pop else addCells ((n-1<Types.cell>),a) (a::pop) 

let dilutionTwoDivision numberOfClone1 clone1 numberOfClone2 clone2  = 
    let rec dilutionSim nClone1 clone1 nClone2 clone2 count acc =
        if count<>(nClone1+nClone2) then
            let clone = if count < nClone1 then clone1 else clone2
            let result = simulate {clone with rng=System.Random(count)}
            //Randomise the initial label concentration
            let rec initialLabel () =   let result = RandomNumbers.gaussianMP (System.Random(count)) 1. 0.2 
                                        if result <= 0. then initialLabel () else result
            let cellDilution = List.map (fun (i:populationState) -> 
                match i.population.dilution with
                | Divisions(n) -> (i.population.basal,(initialLabel ())/float(n+1)) 
                | _ -> failwith "Called division function with wrong type of clone" ) result
            dilutionSim nClone1 clone1 nClone2 clone2 (count+1) (cellDilution::acc)
        else
            acc 
            |> List.fold (fun population individual ->   match population with
                                                                                            | None ->  Some(List.map (fun i -> [i]) individual)
                                                                                            | Some(partial) -> Some(List.map2 (fun a i -> i::a ) partial individual)
                                                                                            ) None 
            |> (fun a -> match a with
                         | None     -> [] 
                         | Some(l)  -> List.map (fun t ->   List.filter (fun (p,i) -> p>0<Types.cell> ) t ) l
                                       |> List.map (fun t -> List.fold (fun p c -> addCells c p ) [] t  )  
                         )
    dilutionSim numberOfClone1 clone1 numberOfClone2 clone2 0 []

let dilutionSingleDivision numberOfClones clone = 
    dilutionTwoDivision numberOfClones clone 0 clone

let intensityToText (I:float list list) name = 
    use textFile = new System.IO.StreamWriter(name, true)
    List.iter (fun It -> List.iter (fun Iti -> textFile.WriteLine(sprintf "%f" Iti)) It; textFile.WriteLine("")) I
    textFile.Close()
    
let dilutionSinglePopulation numberOfClones (clone:clone) =
    match clone.state.population.dilution with
    | DontMeasure       -> failwith "Can't estimate dilution without a method"
    | Divisions(_)      -> dilutionSingleDivision numberOfClones clone
    | Population(_)     -> dilutionSinglePopulation numberOfClones clone