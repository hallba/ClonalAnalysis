﻿module SupraBasalAnalysis

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

let noFiltration result = Some(result)

let filterByFirstTimepointDivision minDivisions (result:populationState list) = 
    let t0 = result.[1]
    match t0.population.dilution with
    | Divisions(n) -> if n < minDivisions then None else Some(result)
    | _ -> failwith "Need simulations to count divisions in order to filter by divisions" 

let floatingCloneProbability' finishingCondition clone filtration = 
    if clone.SBRatio = None then failwith "Cannot calculate floating clone probabilities without an exit rate"
    let calcP (c,fc) = float(fc)/float(c)
    let rec fcsim count clone finishingCondition acc =
        let result = simulate {clone with rng=System.Random(count)} |> filtration
        match result with
        | None -> fcsim (count+1) clone finishingCondition acc
        | Some(result) ->
            let acc' = List.map2 (fun tc tn -> countFloatingClones tc tn) acc result
            match finishingCondition with
            | Count(i) ->   //printf "Count %A Acc %A\n" count acc'
                            if count >= (i-1) then acc' else fcsim (count+1) clone finishingCondition acc'
            | Tolerance(f) ->   let maxDiff = List.map2 (fun aE aE' -> abs((calcP aE')-(calcP aE))) acc acc' |> List.max
    //                            let t0::trest = acc'
    //                            let mino = List.minBy (fun (c,fc) -> fc) trest
                                //ignore (if count%10000=0 then printf "Count %A maxDiff %A\n" count maxDiff else ())
                                if (maxDiff > f|| System.Double.IsNaN(maxDiff) || maxDiff=0.) then fcsim (count+1) clone finishingCondition acc' else acc'
    
    let neutral = List.map (fun i -> (0,0)) (   match clone.reporting with 
                                                | Specified(a) -> ((0.<Types.week>)::a)
                                                | Regular(r) -> ((0.<Types.week>)::(List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)))
                                                )
    fcsim 0 clone finishingCondition neutral
    |> List.map calcP

let floatingCloneProbability numberSimulations clone = 
    let sims = Seq.init numberSimulations (fun i -> simulate {clone with rng=System.Random(i)})
    let neutral = List.map (fun i -> (0,0)) (   match clone.reporting with 
                                                | Specified(a) -> ((0.<Types.week>)::a)
                                                | Regular(r) -> ((0.<Types.week>)::(List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)))
                                                )
    let r = sims |>  Seq.fold ( fun acc observation -> List.map2 (fun tc tn -> countFloatingClones tc tn) acc observation ) neutral
    printf "R %A" r
    r |> List.map (fun (c,fc) -> float(fc)/float(c) )
//    Array.fold ( fun (popState:populationState) acc ->         
//                                                                        if popState.population.basal = 0 then acc else
//                                                                            let count = snd(acc)
//                                                                            let ratio = fst(acc) + float(popState.population.basal)/float(popState.population.suprabasal)
//                                                                            (count,ratio)          
//                                                                        ) neutral

let logLikelihood p o = 
    let fc = fst(o)
    let nfc = snd(o)
    if (fc+nfc=0) then 0. else (log(p)*float(fc) + log(1.-p)*float(nfc))

let countSBIndividual popState total basalSizeIndividual =
    if basalSizeIndividual <> popState.population.basal then total
    else if popState.population.suprabasal = 0<Types.cell> then ((fst(total)+1),(snd(total)))
    else ((1+fst(total)),(1+snd(total)))

let suprabasalObservationProbability finishingCondition basalSizes clone =
    let calcP (c,fc) = float(fc)/float(c)
    let rec sbsim clone finishingCondition count acc =
        let sim =   simulate {clone with rng=System.Random(count)}
        let acc' = List.map2 (fun popState currentState -> List.map2 (countSBIndividual popState) currentState basalSizes) sim acc

        match finishingCondition with
        | Count(n) -> if count >= (n-1) then acc' else sbsim clone finishingCondition (count+1) acc'
        | Tolerance(a) ->   let result = List.map (fun timePoint -> List.map calcP timePoint ) acc
                            let result' = List.map (fun timePoint -> List.map calcP timePoint ) acc'
                            let maxDiff = List.map2 (fun timePoint' timePoint -> List.map2 (fun b' b -> let r = abs(b' - b)
                                                                                                        if System.Double.IsNaN(r) then 0. else r ) timePoint' timePoint ) result' result
                                          |> List.map (fun e -> List.max e)
                                          |> List.max
                            if (maxDiff > a|| System.Double.IsNaN(maxDiff) || maxDiff=0.) then sbsim clone finishingCondition (count+1) acc' else acc'
    
    let timePoints =    match clone.reporting with 
                        | Specified(a) -> ((0.<Types.week>)::a)
                        | Regular(r) -> ((0.<Types.week>)::(List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)))
                                                                                                
    let neutral = List.map (fun j -> (List.map (fun i -> (0,0)) basalSizes))  timePoints
    
    sbsim clone finishingCondition 0 neutral
    |> List.map (fun timePoint -> List.map calcP timePoint )

let supraBasalRatio finishingCondition clone =
    let rec restructureResult result (inp:float option list list) =
        match inp with
        | [] -> result
        | sim::rest ->  let result' = List.map2 (fun timepoint indiSim ->   match indiSim with 
                                                                            | None -> timepoint
                                                                            | Some(i) -> i::timepoint) result sim
                        restructureResult result' rest
    let calculateRatio (i:populationState) = 
        if i.population.basal > 0<Types.cell> && i.population.suprabasal >0<Types.cell> then
            Some(float(i.population.suprabasal)/float(i.population.basal))
        else
            None
    let rec sbsim clone finishingCondition acc seed =
        let sim =   simulate {clone with rng=System.Random(seed)} |> List.map calculateRatio
        let acc' = sim::acc
        match finishingCondition with
        | Tolerance(a) -> failwith "Not implemented yet"
        | Count(n) -> if n>0 then sbsim clone (Count(n-1)) acc' (seed+1) else acc' |> restructureResult (List.map (fun i -> []) acc'.Head)
    sbsim clone finishingCondition [] 0

let fcPerSingle finishingCondition clone =
    let updateAcc acc (i:populationState) = 
        if i.population.basal > 1<Types.cell> then
            acc
        else if i.population.basal = 1<Types.cell> then
            (fst(acc),(snd(acc)+1))
        else if i.population.suprabasal > 0<Types.cell> then
            ((fst(acc)+1),snd(acc))
        else acc
    let rec sbsim clone finishingCondition acc seed =
        let acc' =   simulate {clone with rng=System.Random(seed)} |> List.fold updateAcc acc
        match finishingCondition with
        | Count(n) -> if n>0 then sbsim clone (Count(n-1)) acc' (seed+1) else acc'
        | _ -> failwith "Not implemented yet"
    sbsim clone finishingCondition (0,0) 0