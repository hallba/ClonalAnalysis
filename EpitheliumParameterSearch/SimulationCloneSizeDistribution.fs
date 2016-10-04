module SimulationCloneSizeDistribution

//open FSharp.Collections.ParallelSeq
open MathNet.Numerics

type dilutionMeasure = Divisions of int | DontMeasure | Population of float []

type basalEvent = Stratification | Division

type finishingCondition = Tolerance of float | Count of int

let dilutionUpdate dMeasure event (rng: System.Random) = 
    match dMeasure with
    |DontMeasure                ->  dMeasure
    |Divisions(n)               ->  if event = Division then Divisions(n+1) else dMeasure
    |Population(p)              ->  let l = Array.length p
                                    //Randomly pick one element
                                    let r = rng.Next(l)
                                    if event = Stratification then 
                                        Population(Array.append p.[..(r-1)] p.[(r+1)..]) 
                                        else    let rec mixture () = let i = RandomNumbers.gaussianMP rng 0.5 0.2
                                                                     if i >= 1. || i <= 0. then mixture () else i
                                                let balance = mixture ()
                                                let d' = p.[r]*balance
                                                let p' = Array.init (Array.length p) (fun i -> if i <> r then p.[i] else d')
                                                Population(Array.append p [|p.[r]*(1.-balance)|])


type cellPopulation = 
                {   A : int<Types.cell>;
                    B : int<Types.cell>;
                    C : int<Types.cell>;  
                    dilution: dilutionMeasure
                    } with
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


type randomBasal =
            {       time: float<Types.week>;
                    basalSum: int<Types.cell> array;
                    animalID:int;
                    parameterSet:string;
            }

type reportStyle = Regular of regularReporting | Specified of float<Types.week> list

type clone = {  state   : populationState;
                lambda  : float<Types.cell/Types.week>;
                rho     : float;
                r       : float;
                delta   : float;
                SBRatio : float option;
                rng     : System.Random
                maxN    : int; //Maximum *requested* number of cells
                reporting   : reportStyle
                report  : populationState list option;
                finalState : populationState option
                }
                with
                member this.gamma = this.lambda * this.rho / (1. - this.rho) //Rate of B->C
                member this.mu = match this.SBRatio with
                                 //BH this rate needs to be scaled with 1-rho to ensure complete replacement of SB population
                                 | Some(r) -> (1.-this.rho)*this.gamma/r //Rate of C loss
                                 | None -> 0.<Types.cell/Types.week> //If we aren't interested in SB cells we can pass Nne and ignore mu's contribution
                member this.eventRate = this.lambda * float(this.state.population.A) * 1.<Types.cell> + this.gamma * float(this.state.population.B) * 1.<Types.cell> + this.mu * float(this.state.population.C) * 1.<Types.cell>
                member this.pAA =   this.r*(1.+this.delta)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.eventRate
                member this.pAB =   (1.-2.*this.r)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.eventRate
                member this.pBB =   this.r*(1.-this.delta)*this.lambda*float(this.state.population.A)*1.<Types.cell>/this.eventRate
                member this.pB2C =  float(this.state.population.B)*1.<Types.cell>*this.gamma/this.eventRate
                member this.pCExit = float(this.state.population.C)*1.<Types.cell>*this.mu/this.eventRate

                member this.finishedCondition (p:cellPopulation) = 
                    match this.SBRatio with 
                    | Some(_) -> (p.basal+p.suprabasal) > 0<Types.cell>
                    | None -> p.basal > 0<Types.cell>

                member private this.recordPreviousStates reportTimes currentTime =
                    let rec core reportTimes currentTime acc = 
                        match reportTimes with
                        | earliest::rest when earliest < currentTime -> core rest currentTime ({ population = this.state.population ; time=earliest }::acc)
                        | _ -> ((List.rev acc),reportTimes)
                    core reportTimes currentTime []
                member this.selectEvent =   
                    let random = this.rng.NextDouble()
                    //AA
                    if random < this.pAA then 
                        {this.state.population with A=this.state.population.A+1<Types.cell>;dilution=(dilutionUpdate this.state.population.dilution Division this.rng) }
                    //AB
                    else if random < this.pAA + this.pAB then 
                        {this.state.population with B=this.state.population.B+1<Types.cell>;dilution=(dilutionUpdate this.state.population.dilution Division this.rng)}
                    //BB
                    else if random < this.pAA + this.pAB + this.pBB then 
                        {this.state.population with A=this.state.population.A-1<Types.cell>;B=this.state.population.B+2<Types.cell>;dilution=(dilutionUpdate this.state.population.dilution Division this.rng)}
                    //Migration
                    else if random < this.pAA + this.pAB + this.pBB + this.pB2C then 
                        {this.state.population with B=this.state.population.B-1<Types.cell>;C=this.state.population.C+1<Types.cell>;dilution=(dilutionUpdate this.state.population.dilution Stratification this.rng)}
                    //Shedding
                    else {this.state.population with C=this.state.population.C-1<Types.cell>}
                member this.update =    //If the system has run out of cells, just update the final state and return the clone
                    match (this.finalState,this.reporting) with
                    | (None,Regular(r)) ->
                        //Skip an update step but convert the clone into a specified times
                        let timings = List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)
                        {this with reporting=Specified(timings)}
                    | (Some(finalState),Regular(r)) -> 
                        //Regular reporting events are converted to Specified so this should never arise
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
                            //Simulation finished condition- if there are no more cells, report the final state else continue
                            if this.finishedCondition population' then {this with state = {population = population'; time = time'} ; reporting=Specified(l') ; report = report' }
                            else
                                {this with state = {population = population'; time = time'} ; reporting=Specified(l') ; report = report' ; finalState = Some({time=time'; population=population'}) }
                    | (Some(state),Specified(l)) -> 
                        //This should only be called if a simulation has completed but further reports are still required
                        match l with
                        | [] -> failwith "Cannot extend a completed simulation trace"
                        | nextReport::later ->
                                        let time' = nextReport
                                        let report' = {state with time=nextReport}
                                        {this with state=report'; finalState = Some(report'); reporting=Specified(later); report=Some([report'])}

                                                                
let initClone = {   state = {   population = {  A = 1<Types.cell>
                                                B = 0<Types.cell>
                                                C = 0<Types.cell>; 
                                                dilution = DontMeasure}
                                time =0.<Types.week> }
                    lambda  = 2.<Types.cell/Types.week>
                    rho = 0.85
                    r = 0.15
                    delta = 0.
                    maxN = 10
                    rng = System.Random(1982) //Specified a seed to make testing possible
                    SBRatio = None //Some(0.27) //100 B : 19 SB1 : 8 SB2
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

let summarizeObservations (allBasalSums: randomBasal [] )  =

    let filename = allBasalSums.[0].parameterSet
    let outputFile = @"//datacentre/Shares/Users/vk325/Desktop/test/"+filename

    let summary = 
      [| for i in allBasalSums ->
            let daytime = int (i.time * 7.0)
            let timeStr = "tExp(k) = "+daytime.ToString()+"/7;"
            let animalIDStr = "% "+daytime.ToString()+" days, mouse"+ i.animalID.ToString()
            let max = i.basalSum |> Seq.max   //maximum cell number observed
            let res = i.basalSum //|> Seq.filter (fun x-> int(x) <> 0)    // discard zero observations
                                 |> Seq.countBy (fun x-> x)
                                 |> Seq.toList
                                 |> List.sort

            let occs = Array.init (int(max)) (fun x-> x+1)
            //let occs = Array.init (int(max)+1) (fun x-> x)

            let obsSummary = 
                occs 
                |> Array.map (fun oc -> 
                    res
                    |>List.filter (fun (f,s)-> oc = (int)f))
                    |>Array.mapi (fun i el ->
                        match el with
                        | [] -> sprintf "%A %A" occs.[i] 0
                        | [(f,s)] -> sprintf "%A %A" f s
                        | _ -> failwith "No observation created!"
                        )
          
            Array.concat [| [|animalIDStr|]; [|"k=k+1;"|]; [|timeStr|]; [|"data{k} = ..."|]; [|"["|]; obsSummary; [|"];"|]; [|"\n"|] |]
             |]
        |> Array.concat
        
    summary|> (fun s -> System.IO.File.AppendAllLines(outputFile, s))
                                          
         

//calculate correlation coefficient (r^2) and linear equation
let calcLinear (x:float[], y:float[]) = 
    let linearParameters = Fit.Line(x,y)   //tuple:(fst:intercept,snd:slope)
    let rsq = GoodnessOfFit.RSquared(x,y)
    (linearParameters,rsq)
 
//calculate the slope from t=0, n=1 (x1,y1) and first non t=0 timepoint (x2,y2)    
let calcSlope (x2:float, y2:float) = 
    let x1 = 0.0
    let y1 = 1.0
    (y2 - y1) / (x2 - x1)

let shortTimeSanityCheck (timePoints: float [] list, avgBasalCells, clone) =
    //BH- hard coded variables should be passed to functions
    let outputFile = @"//datacentre/Shares/Users/vk325/Desktop/sanity_check/short_timescales.txt"

    let gamma = (clone.lambda*clone.rho) / (1.- clone.rho)
    let (intercept,slope),rsq = calcLinear(timePoints.[0].[1..], avgBasalCells)
    let earlySlope = calcSlope(timePoints.[0].[1], avgBasalCells.[0])
    let inputParameters = clone.r.ToString()+"\t"+clone.rho.ToString()+"\t"+clone.lambda.ToString()+"\t"

    let avg = [|for a in avgBasalCells -> a.ToString() |] |> String.concat "\t"

    System.IO.File.AppendAllLines(outputFile,[|inputParameters+(1./float(gamma)).ToString()+"\t"+System.Math.Round(rsq,4).ToString()+"\t"+slope.ToString()+"\t"+intercept.ToString()+"\t"+earlySlope.ToString()+"\t"+avg|])

    let  shortTimeRule =
        timePoints
        |> List.map (fun times ->
            times.[1..]
            |> Array.map (fun t ->
                (1. + float(clone.lambda) * t).ToString() )
            |> String.concat "\t")
  
    System.IO.File.AppendAllLines(outputFile,shortTimeRule)

let longTimeSanityCheck (timePoints: float [] list, avgBasalCells, clone ) =
    //BH- hard coded variables should be passed to functions
    let outputFile = @"//datacentre/Shares/Users/vk325/Desktop/sanity_check/long_timescales.txt"
    
    let inputParameters = clone.r.ToString()+"\t"+clone.rho.ToString()+"\t"+clone.lambda.ToString()+"\t"
    let ratio = (clone.r *clone.lambda) / clone.rho
    let (intercept,slope),rsq = calcLinear(timePoints.[0].[1..], avgBasalCells)

    let avg = [|for a in avgBasalCells -> a.ToString() |] |> String.concat "\t"
    System.IO.File.AppendAllLines(outputFile,[|inputParameters+ratio.ToString()+"\t"+(1./float(clone.rho)).ToString()+"\t"+System.Math.Round(rsq,4).ToString()+"\t"+System.Math.Round(slope,4).ToString()+"\t"+System.Math.Round(intercept,4).ToString()+"\t"+avg|])


type BasalAvgType = All | Surviving //select type of average number of Basal Cells: All (survivng + exctinct) or surviving only (i.e zero clones excluded)

let massiveSimulation numberSims (average:BasalAvgType) clone =
    let sims = Seq.init numberSims (fun i -> simulate {clone with rng=System.Random(i)})

    let totalTimePoints = match clone.reporting with
                            | Specified(n) -> List.length n
                            | Regular(r) -> (int(r.timeLimit/r.frequency))

    let acc = Array.init totalTimePoints (fun i -> int64(0<Types.cell>))
   
    let avgBasal = match average with
                        | All -> 
                                let addObservations runningTotal (result:populationState list) =
                                    let basalSummary = List.map (fun individualTimePoint -> individualTimePoint.population.basal) result
                                                        |> Array.ofList

                                    Array.map2 (fun r T -> int64(r)+int64(T)) basalSummary.[1..] runningTotal 

                                let totalBasalCells = Seq.fold addObservations acc sims
                                totalBasalCells |> Array.map(fun elem -> float(elem) / float(numberSims))
                        | Surviving ->
                                let addObservations runningTotal (result:populationState list) =
                                    let basalSummary = List.map (fun individualTimePoint -> individualTimePoint.population.basal) result
                                                        |> Array.ofList

                                    Array.map2 (fun r T -> int64(r)+int64(T)) basalSummary.[1..] runningTotal 

                                let totalBasalCells = Seq.fold addObservations acc sims

                                let countNonZeros runningTotal (result:populationState list) =
                                    let basalSummary = List.map (fun individualTimePoint -> individualTimePoint.population.basal) result
                                                        |>Array.ofList
                                    Array.map2 (fun r T -> 
                                                    match r with
                                                    |0<Types.cell> -> int64(T)
                                                    |_ -> int64(T)+int64(1)) basalSummary.[1..] runningTotal

                                let nonZeros = Seq.fold countNonZeros acc sims
                                Array.map2(fun tb nz -> float(tb) / float(nz)) totalBasalCells nonZeros

    let getTimePoints (result:populationState list) = 
        result|>List.map(fun t -> float(t.time))|>Array.ofList
        
    let timepoints = sims|> Seq.map (fun s -> getTimePoints s)|> Seq.take 1|> Seq.toList
    
    //BH- Commented out code that generates output to files with hardcoded filenames. Pass an option for output?
    //shortTimeSanityCheck(timepoints, avgBasal, clone)
    //longTimeSanityCheck(timepoints, avgBasal, clone)
    //BH- return timepoints and averages as arrays of equal length
    let processedTimepoints = timepoints.[0] |> List.ofArray |> fun i ->    match i with 
                                                                            | t0::rest -> Array.ofList rest
                                                                            | [] -> failwith "At least one timepoint should have been calculated"
    (processedTimepoints, avgBasal)
 

let getBasalSum timepoint id sims cloneNumberPerAnimal (rnd:System.Random) pfilename = 
    let timeIndex = timepoint + 1 //discard 0 time point, start from second time point (timpoint index: 1)
    let animalID = id + 1
    

    let allClonesPerTimepoint = sims|>Array.map (fun (el:populationState list)->el.[timeIndex]) // get all clone objects that correspond to each timepoint

    let timepoint = allClonesPerTimepoint.[0].time

    //filter empty clones out of sims

    let filteredClonesPerTimepoint = Array.filter (fun cl -> cl.population.basal > 0<Types.cell>) allClonesPerTimepoint

    // Get unique random simulation indices 
    let randomIndices = Seq.initInfinite (fun _ -> rnd.Next(filteredClonesPerTimepoint.Length))
    let uniqueRandomIndices = randomIndices
                            |> Seq.distinct
                            |> Seq.take(cloneNumberPerAnimal)
                            |> Seq.toArray

    let bs = Array.map (fun x -> filteredClonesPerTimepoint.[x].population.basal) uniqueRandomIndices   //get the number of basal cells at each random position

    let rb = {time=timepoint; basalSum=bs; animalID=animalID; parameterSet=pfilename}  
   
    rb
    

let getRandomBasal numberOfSims cloneNumberPerAnimal (clone:clone)=

    let sims = Array.init numberOfSims (fun i -> simulate {clone with rng=System.Random(i)})

    let pfilename = "synthetic"+clone.r.ToString()+"_"+clone.rho.ToString()+"_"+clone.lambda.ToString()+".m" //matlab output file name
    let numberOfTimepoints = sims.[0].Length;

    let animals = 4 //introducing 4 different animals to mimic the experimental data

    let outputFile = @"//datacentre/Shares/Users/vk325/Desktop/test/"+pfilename
    
    let parameters = "%PARAMETERS: "+"r: "+clone.r.ToString()+"; rho: "+clone.rho.ToString()+"; lamda: "+clone.lambda.ToString()+"\nk=0;"
    ignore (System.IO.File.AppendAllLines(outputFile, [|parameters|]))

    let allBasalSums = 
        Array.init (numberOfTimepoints-1) ( fun t ->  //the numberOfTimepoints includes the zero time point which has to be discarded 
          Array.init animals (fun id -> getBasalSum t id sims cloneNumberPerAnimal (System.Random(id)) pfilename) // group to 4 animal ids of equal number of elements
                |>summarizeObservations) 


    allBasalSums   


let createSyntheticDataset parameterSets numberOfSims cloneNumberPerAnimal =

    //let outputFile = @"//datacentre/Shares/Users/vk325/Desktop/test.txt"
    //System.IO.File.AppendAllText(outputFile, "k=0;\n")

    let completeSet = Types.createParameterSet {parameterSets with timePoints=[|0.<Types.week>|]}  
    let allRandomBasalSums = 
        Array.map (fun input ->     //parallel causes issues in printing
            input 
            |> parameterSetToClone (Specified(List.ofArray parameterSets.timePoints))
            |> getRandomBasal numberOfSims cloneNumberPerAnimal
        ) completeSet

    allRandomBasalSums

let validateSyntheticDataset parameterSets numberOfSims avgType =

    let completeSet = Types.createParameterSet {parameterSets with timePoints=[|0.<Types.week>|]}  
    let allRandomBasalSums = 
        Array.map (fun input -> 
            input 
            |> parameterSetToClone (Regular({timeLimit=100.0<Types.week>;frequency=4.0<Types.week>;lastReport=0.<Types.week>}))
            |> massiveSimulation numberOfSims avgType
        ) completeSet

    allRandomBasalSums
  
let rhoByCloneSize finishingCondition clone maxSize = 
    let addSim (simState:populationState) acc = 
        if simState.population.basal*(1<Types.cell^-1>) > maxSize then acc else
            //let old = acc.[(simState.population.basal*(1<Types.cell^-1>))]
            failwith "Incomplete"

    let updateAcc (sim:populationState list) acc =
        Array.map2 addSim (Array.ofList sim) acc
    let rec simLoop clone finishingCondition acc seed =
        let sim = simulate {clone with rng=System.Random(seed)}
        let acc' = updateAcc sim acc
        match finishingCondition with
        | Count(n) -> if n>0 then simLoop clone (Count(n-1)) acc' (seed+1) else acc'
        | _ -> failwith "Not implemented yet"

    let timePoints =    match clone.reporting with 
                        | Specified(a) -> ((0.<Types.week>)::a)
                        | Regular(r) -> ((0.<Types.week>)::(List.init (int(r.timeLimit/r.frequency)) (fun i -> float(i+1)*r.frequency)))
    let neutral = (Array.map (fun j -> Array.init maxSize (fun i -> (0,0)) )  (Array.ofList timePoints))
    simLoop clone finishingCondition neutral 0