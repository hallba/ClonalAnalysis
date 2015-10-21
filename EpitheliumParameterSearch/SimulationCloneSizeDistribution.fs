module SimulationCloneSizeDistribution

type cellPopulation = 
                {   A : int<Types.cell>;
                    B : int<Types.cell>;
                    C : int<Types.cell>;  }

type populationState =
                {   population  : cellPopulation;
                    time        : float<Types.week>;
                }

type reportStyle = Frequency of float<Types.week> | Specified of float<Types.week> list

type clone = {  state   : populationState;
                lambda  : float<Types.cell/Types.week>;
                rho     : float;
                r       : float;
                delta   : float;
                rng     : System.Random;
                reportFrequency  : float<Types.week>;
                lastReportTime   : float<Types.week>;
                report  : populationState option;
                finalState : populationState option
                }
                with
                member this.gamma = this.lambda * this.rho / (1. - this.rho)
                member this.R = this.lambda * float(this.state.population.A) * 1.<Types.cell> + this.gamma * float(this.state.population.B) * 1.<Types.cell>
                member this.pAA =   this.r*(1.+this.delta)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.R
                member this.pAB =   (1.-2.*this.r)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.R
                member this.pBB =   this.r*(1.-this.delta)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.R
                member this.pB2C =  float(this.state.population.B)*1.<Types.cell>*this.gamma/this.R
                member this.selectEvent =   let random = this.rng.NextDouble()
                                            //AA
                                            if random < this.pAA then {this.state.population with A=this.state.population.A+1<Types.cell>}
                                            //AB
                                            else if random < this.pAA + this.pAB then {this.state.population with B=this.state.population.B+1<Types.cell>}
                                            //BB
                                            else if random < this.pAA + this.pAB + this.pBB then {this.state.population with A=this.state.population.A-1<Types.cell>;B=this.state.population.B+2<Types.cell>}
                                            //Migration
                                            else {this.state.population with B=this.state.population.B-1<Types.cell>;C=this.state.population.C+1<Types.cell>}
                member this.update =    //If the system has run out of stem cells, just update the final state and return the clone
                                        match this.finalState with
                                        | None ->
                                                    //Update time, report value
                                                    let dt = - 1.<Types.week> * log (this.rng.NextDouble()/ (this.R*1.<Types.week Types.cell^-2> ) )
                                                    let time'  = this.state.time + dt;
                                                    let (lastReportTime',report') = if (time' - this.lastReportTime) > this.reportFrequency then 
                                                                                        let reportTime = (time'-(time'%this.reportFrequency))
                                                                                        let state = { population = this.state.population ; time=reportTime }
                                                                                        (reportTime,Some(state)) 
                                                                                    else (this.lastReportTime,None)
                                                    //Selection an action
                                                    let population'= this.selectEvent
                                                    if population'.A + population'.B > 0<Types.cell> then {this with state = {population = population'; time = time'} ; lastReportTime = lastReportTime' ; report = report' }
                                                    else 
                                                        let finalReportTime = (time'-(time'%this.reportFrequency)) + this.reportFrequency
                                                        {this with state = {population = population'; time = time'} ; lastReportTime = lastReportTime' ; report = report' ; finalState = Some({time=finalReportTime; population=population'}) }
                                        | Some(finalState) -> 
                                                    let dt = this.reportFrequency
                                                    let time' = this.state.time + dt
                                                    let (lastReportTime',report') = let reportTime = (time'-(time'%this.reportFrequency))
                                                                                    let state = {finalState with time=reportTime}
                                                                                    (reportTime,Some(state))
                                                    let population' = finalState.population
                                                    {this with state = {population = population'; time = time'} ; lastReportTime = lastReportTime' ; report = report' }

let initClone = {   state = {   population = {  A = 1<Types.cell>
                                                B = 0<Types.cell>
                                                C = 0<Types.cell>; }
                                time =0.<Types.week> }
                    lambda  = 2.<Types.cell/Types.week>
                    rho = 0.85
                    r = 0.15
                    delta = 0.
                    rng = System.Random()
                    reportFrequency = 4.<Types.week>
                    lastReportTime = 0.<Types.week>
                    report = None 
                    finalState = None}

let simulate clone timeLimit = 
    let rec core clone timeLimit trace =
        match clone.state.time > timeLimit with
        | true              ->  List.rev trace
        | false             ->      let clone' = clone.update 
                                    match clone'.report with
                                    | None -> core clone' timeLimit trace
                                    | Some(state) -> core clone' timeLimit (state::trace)
    core clone timeLimit [clone.state]

let cloneProbability clone number timeLimit=
    let simulations = Array.Parallel.init number (fun i -> simulate {clone with rng=System.Random(i)} timeLimit)
    ()
