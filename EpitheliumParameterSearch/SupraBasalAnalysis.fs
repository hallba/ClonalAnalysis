module SupraBasalAnalysis

open SimulationCloneSizeDistribution

let countCellRatio total popState = 
    if popState.population.basal = 0<Types.cell> then total else
        let count = fst(total) + 1
        let ratio = snd(total) + float(popState.population.suprabasal)/float(popState.population.basal)
        (count,ratio)  

let countFloatingClones total popState = 
    if (popState.population.basal+popState.population.suprabasal) = 0<Types.cell> then total
    else if popState.population.basal > 0<Types.cell> then ((fst(total)+1),snd(total))
    else ((fst(total)+1),snd(total)+1)

let averageSBBRatio numberSimulations clone =
    let sims = Array.init numberSimulations (fun i -> simulate {clone with rng=System.Random(i)})
    let neutral = List.map (fun i -> (0,0.)) sims.[0]
    sims |>  Array.fold ( fun acc observation -> List.map2 (fun tc tn -> countCellRatio tc tn) acc observation ) neutral
    |> List.map (fun (c,r) -> r/float(c) )

type finishingCondition = Tolerance of float | Count of int
    
let floatingCloneProbability finishingCondition clone = 
    let sims = Seq.init numberSimulations (fun i -> simulate {clone with rng=System.Random(i)})
    let neutral = List.map (fun i -> (0,0)) (   match clone.reporting with 
                                                | Specified(a) -> ((0.<Types.week>)::a)
                                                | Regular(r) -> ((0.<Types.week>)::(List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)))
                                                )
    sims |>  Seq.fold ( fun acc observation -> List.map2 (fun tc tn -> countFloatingClones tc tn) acc observation ) neutral
    |> List.map (fun (c,fc) -> float(fc)/float(c) )
//    Array.fold ( fun (popState:populationState) acc ->         
//                                                                        if popState.population.basal = 0 then acc else
//                                                                            let count = snd(acc)
//                                                                            let ratio = fst(acc) + float(popState.population.basal)/float(popState.population.suprabasal)
//                                                                            (count,ratio)          
//                                                                        ) neutral