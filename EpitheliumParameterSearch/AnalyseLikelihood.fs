module AnalyseLikelihood

type individualLikelihood = {
                            value       :   float
                            rIndex      :   int option
                            rhoIndex    :   int option
                            lambdaIndex :   int option
                            deltaIIndex :   int option
                            }

let maxIndex a =
    Array.mapi (fun i v -> (i,v) ) a
    |> Array.fold (fun acc pair ->  let acc0 =  match acc with
                                                | acc0::rest -> acc0
                                                | [] -> failwith "Empty accumulator"
                                    
                                    if snd(acc0) < snd(pair) then [pair] 
                                    else if snd(acc0) = snd(pair) && fst(acc0) <> fst(pair) then pair::acc
                                    else acc
                                     ) [(0,a.[0])]

let maxByIndex f a =
    Array.mapi (fun i v -> (i,v) ) a
    |> Array.fold (fun acc pair ->  let acc0 =  match acc with
                                                | acc0::rest -> acc0
                                                | [] -> failwith "Empty accumulator"
                                    
                                    if f (snd(acc0)) < f (snd(pair)) then [pair] 
                                    else if f (snd(acc0)) = f (snd(pair)) && fst(acc0) <> fst(pair) then pair::acc
                                    else acc
                                     ) [(0,a.[0])]

//Functions required

//Scale and sum 2D plane of L
let scaleAndSum lPlane =
    let max =   Array.map (fun lRow -> Array.max lRow) lPlane
                |> Array.max
    max

//Project into different planes

//Select a plane of a specific lambda

//Calculate avg and stdev of datasets

let compareListsOfIndexValuePairs iv =
    match iv with
    | first::rest -> snd(first)
    | [] -> failwith "EmptyResult"

let globalMax likelihoodMatrix =
    Array.map    ( fun delta -> Array.map ( fun lambda -> Array.map ( fun rho -> maxIndex rho ) lambda ) delta ) likelihoodMatrix 
    |> Array.map ( fun delta -> Array.map ( fun lambda -> ( fun rho -> maxByIndex ( fun iv ->  compareListsOfIndexValuePairs iv ) rho) lambda  ) delta )
    //|> Array.map ( fun delta -> Array.map ( fun lambda -> maxByIndex (fun a -> compareListsOfIndexValuePairs a) lambda) )
    