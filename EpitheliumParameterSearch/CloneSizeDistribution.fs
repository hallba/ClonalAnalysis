module CloneSizeDistribution

open MathNet.Numerics

[<Measure>]
type cell

[<Measure>]
type week

[<Measure>]
type probability

let F x y t r gamma = 
    let v = sqrt(1.-4.*r)
    let w = (gamma*(1.-2.*r) - 2.*r) / (2.*gamma*sqrt(1.-4.*r))
    let gg = complex (sqrt(1.-4.*r)/gamma) 0.

    let u = ((complex 1. 0.)-y)*exp(complex (-gamma*t) 0.) 
    let u0 = ((complex 1. 0.)-y )
    
    //printf "v %A w %A gg %A u %A u0 %A\n" v w gg u u0

    let Q = (complex (1. + 2.*w) 0.) - gg*u0 + ( (complex (2.*r) 0.)*(x-y) + y - (complex 1. 0.) )/(complex gamma 0.)
    
    //Interesting bug. Forget to put the brackets around the  gg*u0- type checker does not find this, gives different result
    let C = ((-Q*(Whittaker.M w 0. (gg*u0)) ) + (complex (1.+2.*w) 0.)*(Whittaker.M (1.+w) 0. (gg*u0))) / (Q*(Whittaker.W w 0. (gg*u0) )+(complex 2. 0.)*(Whittaker.W (1.+w) 0. (gg*u0) ) )

    (complex 1. 0.) - u + 
    (u*(complex (1.+v) 0.)- (complex (gamma*(1.+2.*w)) 0.) ) / complex (2.*r) 0. +
    (complex (gamma/(2.*r)) 0.) * ( (complex (1.+2.*w) 0.)* (Whittaker.M (1.+w) 0. (u*gg) ) - (complex 2. 0.)*C*(Whittaker.W (1.+w) 0. (u*gg) ) ) /
    ( (Whittaker.M w 0. (u*gg)) + C* (Whittaker.W w 0. (u*gg)) )    

type parameterSet = {   time: float<week>;
                        rho: float; 
                        r: float<probability>;
                        lambda: float<cell/week>
                        } with 
                        member this.migration = this.rho/(1.-this.rho) * this.lambda
                        member this.correctPathologicalPoint = { this with rho = this.rho + 0.0001 ; r = this.r+0.00001<probability> }
                        override this.ToString() = sprintf "Rho: %A r:%A Time:%A Lambda:%A" this.rho this.r this.time this.lambda

let testSystem = {time=1.<week>; rho=0.85; r=0.15<probability>; lambda=2.<cell/week>}

let probabilityClonalSurvival inputParameters =
    //Note that the below line is incorrect, and the effects corrected later
    let gamma = 1./inputParameters.rho - 1.
    let T = inputParameters.lambda * inputParameters.time
    let p = (complex 1. 0.) - F (complex 0. 0.) (complex 0. 0.) (T*1.<cell^-1>) (inputParameters.r*1.<probability^-1>) gamma
    //printf "P: %A\n" p
    //The following line checks for "pathological points" as described in the original matlab implementation
    if p.r < 0. || p.r > 1. || System.Double.IsNaN(p.r) || System.Double.IsInfinity(p.r) then ((complex 1. 0.) - F (complex 0. 0.) (complex 0. 0.) (T*1.<cell^-1>) ((inputParameters.r*1.<probability^-1>)+0.00001) (gamma+0.00001)) else p

let complexSum c =
    let rec core c acc = 
        match c with
        | top::rest -> core rest (acc+top)
        | [] -> acc
    core c (complex 0. 0.)
    
let rec probabilityCloneSizes inputParameters nRange maxN =
    let N = if List.max nRange >= maxN then List.max nRange + 1 else maxN
    let gamma = 1./inputParameters.rho - 1.
    let T = inputParameters.lambda * inputParameters.time
    let k = List.init N (fun item -> exp(complex 0. (float(item)/float(N)*System.Math.PI*2.)) )
    let gVals = List.map (fun zMember -> if zMember = complex 1. 0. then complex 1. 0. else F zMember zMember (T*1.<cell^-1>) (inputParameters.r*1.<probability^-1>) gamma) k 
    let p = List.map (fun n -> List.mapi (fun index i -> i * exp(complex 0. (-2.*System.Math.PI*float(n)*float(index)/float(N)) ) ) gVals ) nRange 
            |>List.map (fun item -> (complexSum item).r/float(N) )
            |>List.map (fun item -> if item < 0. then 0. else item) //In the matlab code we had an additional test which could allow values less
                                                                    //than -1e-5 to be included in the output. Perhaps include in future?
    let totalP = List.sum p
    //Test for pathological points
    if totalP > 1. || totalP < 0. || System.Double.IsNaN(totalP) then 
        //NB: this is from the matlab code- note that it is a specified small perturbation, not 0.1%
        ignore (printf "Pathological point: %A. Making 0.1%% perturbation" inputParameters); probabilityCloneSizes inputParameters.correctPathologicalPoint nRange maxN else
        p

