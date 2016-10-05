module UnitTests

open SimulationCloneSizeDistribution

//Homeostasis
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
    result

let homeostasis_rho65_r25 = 
    let tissueSize = {initClone.state.population with A=1400<Types.cell>;B=600<Types.cell>;C=2000<Types.cell>}
    let tissue = {initClone with    state={initClone.state with population=tissueSize};
                                    rho=0.65;
                                    SBRatio=Some(1.);
                                    r=0.25;
                                    lambda=2.<Types.cell/Types.week>
                                    }
    let result = simulate tissue
    result