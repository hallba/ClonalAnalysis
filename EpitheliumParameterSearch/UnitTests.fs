﻿module UnitTests

open SimulationCloneSizeDistribution

type Problem = RhoMismatch of populationState | MMismatch of populationState

type ElementResult = Fine of populationState | Bad of Problem

type TestResult = Pass | Fail of ElementResult list

let floatRatio a b n = 
    //Convienience function for calculating deviation of a ratios of ints (a,b) compared to another float (n)
    (float(a)/float(b))/n
    
//Homeostasis
let testHomeoState (c:clone) (p:populationState) =
    let sbH =   match c.SBRatio with 
                | None -> None
                | Some(n) ->    let actualRatio = floatRatio p.population.suprabasal p.population.basal n
                                if actualRatio <1.1 && actualRatio > 0.9 then None else Some(n,p.population.suprabasal,p.population.basal,(floatRatio p.population.suprabasal p.population.basal n))
    let actualRho = floatRatio p.population.A p.population.B (c.rho/(1.-c.rho))  
    let rhoH =  if actualRho <1.1 && actualRho > 0.9 then None else Some((c.rho/(1.-c.rho)),p.population.A,p.population.B,(floatRatio p.population.A p.population.B (c.rho/(1.-c.rho))) )
    //(sbH,rhoH)
    match sbH,rhoH with
    | None,None -> Fine(p)
    | Some(_),_ -> Bad(MMismatch(p))
    | _ -> Bad(RhoMismatch(p))

let testHomeostasis (result:populationState list) (c:clone) =
    List.map (testHomeoState c) result 

let parseHStasisResult i =  match i with
                            | Bad(n) -> true
                            | Fine(n) -> false

let homeostasis_rho50_r25 =
    let tissueSize = {initClone.state.population with A=1000<Types.cell>;B=1000<Types.cell>;C=2000<Types.cell>}
    let tissue = {initClone with    state={initClone.state with population=tissueSize};
                                    rho=0.5;
                                    SBRatio=Some(1.);
                                    r=0.25;
                                    lambda=2.<Types.cell/Types.week>
                                    }

    let result = simulate tissue
    //Need a function here to test that result- should have A~B C~A+B
    let test = testHomeostasis result tissue 
    test
    |> List.filter parseHStasisResult
    |> (fun i -> if List.length i > 0 then Fail(test) else Pass)

let homeostasis_rho65_r25 = 
    let tissueSize = {initClone.state.population with A=1300<Types.cell>;B=700<Types.cell>;C=2000<Types.cell>}
    let tissue = {initClone with    state={initClone.state with population=tissueSize};
                                    rho=0.65;
                                    SBRatio=Some(1.);
                                    r=0.25;
                                    lambda=2.<Types.cell/Types.week>
                                    }
    let result = simulate tissue
    let test = testHomeostasis result tissue
    test
    |> List.filter parseHStasisResult
    |> (fun i -> if List.length i > 0 then Fail(test) else Pass)