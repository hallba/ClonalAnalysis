module SimulationCloneSizeDistribution

type cellPopulation = 
                {   A : int<Types.cell>;
                    B : int<Types.cell>;
                    C : int<Types.cell>;  } with
                    member this.basal = this.A + this.B
                    member this.suprabasal = this.C

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
                SBRatio : float;
                rng     : System.Random
                maxN    : int; //Maximum *requested* number of cells
                reporting   : reportStyle
                report  : populationState list option;
                finalState : populationState option
                }
                with
                member this.gamma = this.lambda * this.rho / (1. - this.rho) //Rate of B->C
                member this.mu = this.gamma/this.SBRatio //Rate of C loss
                member this.eventRate = this.lambda * float(this.state.population.A) * 1.<Types.cell> + this.gamma * float(this.state.population.B) * 1.<Types.cell> + this.mu * float(this.state.population.C) * 1.<Types.cell>
                member this.pAA =   this.r*(1.+this.delta)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.eventRate
                member this.pAB =   (1.-2.*this.r)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.eventRate
                member this.pBB =   this.r*(1.-this.delta)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.eventRate
                member this.pB2C =  float(this.state.population.B)*1.<Types.cell>*this.gamma/this.eventRate
                member this.pCExit = float(this.state.population.C)*1.<Types.cell>*this.mu/this.eventRate
                member private this.recordPreviousStates reportTimes currentTime =
                    let rec core reportTimes currentTime acc = 
                        match reportTimes with
                        | earliest::rest when earliest < currentTime -> core rest currentTime ({ population = this.state.population ; time=earliest }::acc)
                        | _ -> ((List.rev acc),reportTimes)
                    core reportTimes currentTime []
                member this.selectEvent =   
                    let random = this.rng.NextDouble()
                    //AA
                    if random < this.pAA then {this.state.population with A=this.state.population.A+1<Types.cell>}
                    //AB
                    else if random < this.pAA + this.pAB then {this.state.population with B=this.state.population.B+1<Types.cell>}
                    //BB
                    else if random < this.pAA + this.pAB + this.pBB then {this.state.population with A=this.state.population.A-1<Types.cell>;B=this.state.population.B+2<Types.cell>}
                    //Migration
                    else if random < this.pAA + this.pAB + this.pBB + this.pB2C then {this.state.population with B=this.state.population.B-1<Types.cell>;C=this.state.population.C+1<Types.cell>}
                    //Shedding
                    else {this.state.population with C=this.state.population.C-1<Types.cell>}
                member this.update =    //If the system has run out of cells, just update the final state and return the clone
                    match (this.finalState,this.reporting) with
                    | (None,Regular(r)) ->
                        //Skip an update step but convert the clone into a specified times
                        let timings = List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)
                        {this with reporting=Specified(timings)}
                    | (Some(finalState),Regular(r)) -> 
                        failwith "Regular reporting should never report a final state"
                    | (None,Specified(l)) -> 
                        match l with
                        | [] -> failwith "Cannot update a completed simulation trace"
                        | nextReport::later -> 
                            //let dt = - 1.<Types.week> * log (this.rng.NextDouble()/ (this.eventRate*1.<Types.week Types.cell^-2> ) )
                            let dt = -1. * log(this.rng.NextDouble()) / this.eventRate * 1.<Types.cell^2>
                            let time'  = this.state.time + dt
                            let (report',l') =  this.recordPreviousStates l time'
                            let report' =   match report' with
                                            | [] -> None
                                            | _ -> Some(report')
                            let population' = this.selectEvent
                            if (population'.basal+population'.suprabasal) > 0<Types.cell> then {this with state = {population = population'; time = time'} ; reporting=Specified(l') ; report = report' }
                            else
                                {this with state = {population = population'; time = time'} ; reporting=Specified(l') ; report = report' ; finalState = Some({time=time'; population=population'}) }
                    | (Some(state),Specified(l)) -> 
                        match l with
                        | [] -> failwith "Cannot extend a completed simulation trace"
                        | nextReport::later ->
                                        let time' = nextReport
                                        let report' = {state with time=nextReport}
                                        {this with state=report'; finalState = Some(report'); reporting=Specified(later); report=Some([report'])}

                                                                
let initClone = {   state = {   population = {  A = 1<Types.cell>
                                                B = 0<Types.cell>
                                                C = 0<Types.cell>; }
                                time =0.<Types.week> }
                    lambda  = 2.<Types.cell/Types.week>
                    rho = 0.85
                    r = 0.15
                    delta = 0.
                    maxN = 10
                    rng = System.Random(1982) //Specified a seed to make testing possible
                    SBRatio = 0.27 //100 B : 19 SB1 : 8 SB2
                    reporting = Regular({timeLimit=200.<Types.week>;frequency=4.<Types.week>;lastReport=0.<Types.week>})
                    report = None 
                    finalState = None}

let parameterSetToClone timePoints (inp : Types.parameterSet) = 
    //This needs to include all of the times to be tested!
    { initClone with lambda = inp.lambda; rho = inp.rho; r = 1.<Types.probability^-1>*inp.r; delta = 1.<Types.probability^-1>*inp.delta; reporting=timePoints ; maxN = inp.maxN}

let specificClone = {initClone with reporting = Specified([1.<Types.week>;2.<Types.week>;4.<Types.week>;8.<Types.week>;12.<Types.week>;26.<Types.week>;52.<Types.week>;78.<Types.week>])}

let hasSimulationFinished clone =
    match clone.reporting with 
    | Regular(a) -> clone.state.time > a.timeLimit
    | Specified(a) -> if a <> [] then false else true

let rec concatAllItems orderedList acc =
    match orderedList with 
    | [] -> acc
    | head::rest -> concatAllItems rest (head::acc)

let simulate clone = 
    let rec core clone trace =
        match hasSimulationFinished clone with
        | true              ->  List.rev trace
        | false             ->      let clone' = clone.update 
                                    match clone'.report with
                                    | None -> core clone' trace
                                    | Some(states) -> core clone' (concatAllItems states trace)
    core clone [clone.state]

let extendArrayForBigObservation arr bigObservation =
    Array.init (bigObservation+1) (fun i -> if i+1 < Array.length arr then arr.[i] else if i = bigObservation then 1 else 0)

let addObservation arr cloneSize = 
    if (Array.length arr) >= (cloneSize*(1<Types.cell^-1>)+1) 
    then    ignore (arr.[cloneSize*(1<Types.cell^-1>)] <- (arr.[cloneSize*(1<Types.cell^-1>)] + 1) )
            arr
    else    extendArrayForBigObservation arr (cloneSize*(1<Types.cell^-1>))

type cloneSizeProbability = 
    {   basalFraction           :   float array
        suprabasalFraction      :   float array
        time                    :   float<Types.week>
        }

let normaliseObservation arr =
    let sumArray = Array.sum arr |> float
    Array.map (fun obs -> float(obs)/sumArray )  arr

type cloneSizeDistribution =
    {   basalObservation        :   int array
        suprabasalObservation   :   int array
        time                    :   float<Types.week>
    } with 
    member this.add (observation:cellPopulation) =  let basal = observation.basal 
                                                    let suprabasal = observation.suprabasal
                                                    let basalObservation' = addObservation this.basalObservation basal
                                                    let suprabasalObservation' = addObservation this.suprabasalObservation suprabasal
                                                    {this with basalObservation=basalObservation';suprabasalObservation=suprabasalObservation'}
    member this.normalise = {basalFraction=normaliseObservation(this.basalObservation);suprabasalFraction=normaliseObservation(this.suprabasalObservation);time=this.time}
        

let noObservations = {basalObservation=[||];suprabasalObservation=[||];time=0.<Types.week>}

let forceArrayLength n instance =
    let n = n + 1 //This is intended to correct for the fact that the first value is the probability of no clones so we need n+1 obs to include prob for n clones
    let basalN = Array.length instance.basalFraction
    let suprabasalN = Array.length instance.suprabasalFraction
    let basal' =    if basalN < n || basalN > n then Array.init n (fun i -> if i < basalN then instance.basalFraction.[i] else 0. ) 
                    else instance.basalFraction
    let suprabasal' = if suprabasalN < n || suprabasalN > n then Array.init n (fun i -> if i < suprabasalN then instance.suprabasalFraction.[i] else 0. ) else instance.suprabasalFraction
    {instance with basalFraction=basal'; suprabasalFraction=suprabasal'} 

let enforceMax n probabilities =
    List.map (fun instance -> forceArrayLength n instance ) probabilities

let cloneProbability number (clone:clone)=
    let observations =  match clone.reporting with
                        | Specified(l)  -> List.map (fun time -> {noObservations with time=time}) ((0.<Types.week>)::l)
                        | Regular(r)    -> List.init (int(r.timeLimit/r.frequency)+1) (fun i -> {noObservations with time=float(i)*r.frequency} )
    let sims = Array.init number (fun i -> simulate {clone with rng=System.Random(i)})
    sims 
    |> Array.fold (fun observations simulation -> List.map2 (fun (o: cloneSizeDistribution) s -> o.add s.population) observations simulation ) observations
    |> List.map (fun timePoint -> timePoint.normalise)
    //The array must be as long or longer than MaxN- extend the clone size probabilities to reflect this
    |> enforceMax clone.maxN   //Add zeros to make arrays as long as maxN requires
                               //Ignore probabilities above maxN as they should never be observed
    |> (fun output ->
        match output with
        | [] -> failwith "No datapoints created!"
        | initPoint::rest -> rest //Discard t0 
        )

