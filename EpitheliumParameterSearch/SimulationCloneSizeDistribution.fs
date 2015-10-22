﻿module SimulationCloneSizeDistribution

type cellPopulation = 
                {   A : int<Types.cell>;
                    B : int<Types.cell>;
                    C : int<Types.cell>;  }

type populationState =
                {   population  : cellPopulation;
                    time        : float<Types.week>;
                }

type regularReporting = 
                {   frequency   :   float<Types.week>
                    timeLimit   :   float<Types.week>
                    lastReport  :   float<Types.week>
                }

type reportStyle = Regular of regularReporting | Specified of float<Types.week> list

type clone = {  state   : populationState;
                lambda  : float<Types.cell/Types.week>;
                rho     : float;
                r       : float;
                delta   : float;
                rng     : System.Random;
                reporting   : reportStyle
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
                                        match (this.finalState,this.reporting) with
                                        | (None,Regular(r)) ->
                                                                //Update time, report value
                                                                let dt = - 1.<Types.week> * log (this.rng.NextDouble()/ (this.R*1.<Types.week Types.cell^-2> ) )
                                                                let time'  = this.state.time + dt;
                                                                let (lastReportTime',report') = if (time' - r.lastReport) > r.frequency then 
                                                                                                    let reportTime = (time'-(time'%r.frequency))
                                                                                                    let state = { population = this.state.population ; time=reportTime }
                                                                                                    (reportTime,Some(state)) 
                                                                                                else (r.lastReport,None)
                                                                //Selection an action
                                                                let population'= this.selectEvent
                                                                if population'.A + population'.B > 0<Types.cell> then {this with state = {population = population'; time = time'} ; reporting=Regular({r with lastReport = lastReportTime'}) ; report = report' }
                                                                else 
                                                                    let finalReportTime = (time'-(time'%r.frequency)) + r.frequency
                                                                    {this with state = {population = population'; time = time'} ; reporting=Regular({r with lastReport = lastReportTime'}) ; report = report' ; finalState = Some({time=finalReportTime; population=population'}) }
                                        | (Some(finalState),Regular(r)) -> 
                                                                let dt = r.frequency
                                                                let time' = this.state.time + dt
                                                                let (lastReportTime',report') = let reportTime = (time'-(time'%r.frequency))
                                                                                                let state = {finalState with time=reportTime}
                                                                                                (reportTime,Some(state))
                                                                let population' = finalState.population
                                                                {this with state = {population = population'; time = time'} ; reporting=Regular({r with lastReport = lastReportTime'}) ; report = report' }
                                        | (None,Specified(l)) -> match l with
                                                                | [] -> failwith "Cannot update a completed simulation trace"
                                                                | nextReport::later -> 
                                                                                        let dt = - 1.<Types.week> * log (this.rng.NextDouble()/ (this.R*1.<Types.week Types.cell^-2> ) )
                                                                                        let time'  = this.state.time + dt
                                                                                        let (report',l') =  if time' > nextReport then 
                                                                                                                                let state = { population = this.state.population ; time=nextReport }
                                                                                                                                (Some(state),later) 
                                                                                                                            else (None,l)
                                                                                        let population' = this.selectEvent
                                                                                        if population'.A + population'.B > 0<Types.cell> then {this with state = {population = population'; time = time'} ; reporting=Specified(l') ; report = report' }
                                                                                        else
                                                                                            {this with state = {population = population'; time = time'} ; reporting=Specified(l') ; report = report' ; finalState = Some({time=time'; population=population'}) }
                                        | (Some(state),Specified(l)) -> match l with
                                                                        | [] -> failwith "Cannot extend a completed simulation trace"
                                                                        | nextReport::later ->
                                                                                        let time' = nextReport
                                                                                        let report' = {state with time=nextReport}
                                                                                        {this with state=report'; finalState = Some(report'); reporting=Specified(later); report=Some(report')}

                                                                
let initClone = {   state = {   population = {  A = 1<Types.cell>
                                                B = 0<Types.cell>
                                                C = 0<Types.cell>; }
                                time =0.<Types.week> }
                    lambda  = 2.<Types.cell/Types.week>
                    rho = 0.85
                    r = 0.15
                    delta = 0.
                    rng = System.Random()
                    reporting = Regular({timeLimit=200.<Types.week>;frequency=4.<Types.week>;lastReport=0.<Types.week>})
                    report = None 
                    finalState = None}

let hasSimulationFinished clone =
    match clone.reporting with 
    | Regular(a) -> clone.state.time > a.timeLimit
    | Specified(a) -> if a <> [] then false else true

let simulate clone = 
    let rec core clone trace =
        match hasSimulationFinished clone with
        | true              ->  List.rev trace
        | false             ->      let clone' = clone.update 
                                    match clone'.report with
                                    | None -> core clone' trace
                                    | Some(state) -> core clone' (state::trace)
    core clone [clone.state]

let cloneProbability clone number timeLimit=
    let simulations = Array.Parallel.init number (fun i -> simulate {clone with rng=System.Random(i)})
    ()
