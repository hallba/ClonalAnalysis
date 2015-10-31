module AnalyseLikelihood

type individualLikelihood = {
                            value       :   float
                            rIndex      :   int
                            rhoIndex    :   int
                            lambdaIndex :   int
                            deltaIndex :   int
                            } with
                            member this.toParameterSet (ans:Types.parameterSearch) = {Types.testSystem with r=ans.rRange.[this.rIndex]
                                                                                                            rho=ans.rhoRange.[this.rhoIndex]
                                                                                                            lambda=ans.lambdaRange.[this.lambdaIndex]
                                                                                                            delta=  match ans.deltaRange with
                                                                                                                    | Types.Zero -> 0.<Types.probability>
                                                                                                                    | Types.Range(ran) -> ran.[this.deltaIndex]
                                                                                                            maxN=ans.maxN

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

let globalMax L =
    Types.resultsMapi (fun a b c d i -> {value=i;deltaIndex=a;lambdaIndex=b;rhoIndex=c;rIndex=d} ) L
    |> Array.map (fun d -> Array.map (fun l -> Array.map (fun rho -> Array.maxBy (fun i -> i.value) rho) l)  d)
    |> Array.map (fun d -> Array.map (fun l -> Array.maxBy (fun i -> i.value) l)  d)
    |> Array.map (fun d -> Array.maxBy (fun i -> i.value) d)
    |> Array.maxBy (fun i -> i.value)