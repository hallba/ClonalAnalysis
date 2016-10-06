module UnitTests

open SimulationCloneSizeDistribution

type TestResult = Pass of populationState list| Fail of populationState list

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
    (sbH,rhoH)

let testHomeostasis (result:populationState list) (c:clone) =
    List.map (testHomeoState c) result 

let parseHStasisResult i =  fst(i) <> None && snd(i) <> None

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
    testHomeostasis result tissue 
    |> List.filter parseHStasisResult
    |> (fun i -> if List.length i > 0 then Fail(result) else Pass(result))

let homeostasis_rho65_r25 = 
    let tissueSize = {initClone.state.population with A=1300<Types.cell>;B=700<Types.cell>;C=2000<Types.cell>}
    let tissue = {initClone with    state={initClone.state with population=tissueSize};
                                    rho=0.65;
                                    SBRatio=Some(1.);
                                    r=0.25;
                                    lambda=2.<Types.cell/Types.week>
                                    }
    let result = simulate tissue
    testHomeostasis result tissue
    |> List.filter parseHStasisResult
    |> (fun i -> if List.length i > 0 then Fail(result) else Pass(result))