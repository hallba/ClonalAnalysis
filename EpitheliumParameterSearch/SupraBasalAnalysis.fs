module SupraBasalAnalysis

open SimulationCloneSizeDistribution

let countCellRatio total popState = 
    if popState.population.basal = 0<Types.cell> then total else
        let count = fst(total) + 1
        let ratio = snd(total) + float(popState.population.basal)/float(popState.population.suprabasal)
        (count,ratio)  

let averageSBBRatio clone =
    let sims = Array.init 100000 (fun i -> simulate {clone with rng=System.Random(i)})
    let neutral = List.map (fun i -> (0,0.)) sims.[0]
    sims |>  Array.fold ( fun acc observation -> List.map2 (fun tc tn -> countCellRatio tc tn) acc observation ) neutral
    |> List.map (fun (c,r) -> r/float(c) )
    
    
    
//    Array.fold ( fun (popState:populationState) acc ->         
//                                                                        if popState.population.basal = 0 then acc else
//                                                                            let count = snd(acc)
//                                                                            let ratio = fst(acc) + float(popState.population.basal)/float(popState.population.suprabasal)
//                                                                            (count,ratio)          
//                                                                        ) neutral