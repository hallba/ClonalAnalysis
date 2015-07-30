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
    
    //printf "x %A t %A r %A gamma %A\n" x t r gamma
    //printf "v %A w %A gg %A u %A u0 %A\n" v w gg u u0

    let Q = (complex (1. + 2.*w) 0.) - gg*u0 + ( (complex (2.*r) 0.)*(x-y) + y - (complex 1. 0.) )/(complex gamma 0.)

    let C = ((-Q*(Whittaker.M w 0. gg*u0) ) + (complex (1.+2.*w) 0.)*(Whittaker.M (1.+w) 0. (gg*u0))) / (Q*(Whittaker.W w 0. gg*u0 )+(complex 2. 0.)*(Whittaker.W (1.+w) 0. (gg*u0) ) )
    
    (complex 1. 0.) - u + 
    (u*(complex (1.+v) 0.)- (complex (gamma*(1.+2.*w)) 0.) ) / complex (2.*r) 0. +
    (complex (gamma/(2.*r)) 0.) * ( (complex (1.+2.*w) 0.)* (Whittaker.M (1.+w) 0. (u*gg) ) - (complex 2. 0.)*C*(Whittaker.W (1.+w) 0. (u*gg) ) ) /
    ( (Whittaker.M w 0. (u*gg)) + C* (Whittaker.W w 0. (u*gg)) )    

let G z t r gamma =
    let rec core z t r gamma acc =
        match z with
        | [] -> acc
        | topZ::rest -> if topZ = complex 1. 0. then core z t r gamma ((complex 1. 0.)::acc) else core z t r gamma ((F topZ topZ t r gamma)::acc)
    core z t r gamma []

//function g = G(z,t,r,gamma)
//
//for k=1:length(z)
//    if(z(k)==1)
//        g(k) = 1;
//    else
//    
//        g(k) = F(z(k),z(k),t,r,gamma);
//        
//    end
//
//end
//end

//function z = u(y,t)
//
//z = (1-y)*exp(-t);
//
//end

type parameterSet = {   time: float<week>;
                        rho: float; 
                        r: float<probability>;
                        lambda: float<cell/week>
                        } with 
                        member this.migration = this.rho/(1.-this.rho) * this.lambda

let testSystem = {time=1.<week>; rho=0.85; r=0.15<probability>; lambda=2.<cell/week>}

let correctPathologicalPoints inputParameters =
    { inputParameters with rho = inputParameters.rho + 0.00001 }

let probabilityClonalSurvival inputParameters =
    //Note that the below line is incorrect, and the effects corrected later
    let gamma = 1./inputParameters.rho - 1.
    let T = inputParameters.lambda * inputParameters.time
    let p = (complex 1. 0.) - F (complex 0. 0.) (complex 0. 0.) (T*1.<cell^-1>) (inputParameters.r*1.<probability^-1>) gamma
    //printf "P: %A\n" p
    //The following line checks for "pathological points" as described in the original matlab implementation
    if p.r < 0. || p.r > 1. || System.Double.IsNaN(p.r) || System.Double.IsInfinity(p.r) then ((complex 1. 0.) - F (complex 0. 0.) (complex 0. 0.) (T*1.<cell^-1>) ((inputParameters.r*1.<probability^-1>)+0.00001) (gamma+0.00001)) else p
    
//let probabilityCloneSizes inputParameters nRange nIter maxN =
//    let N = if List.max nRange > nIter then List.max nRange + 1 else nIter
//    let gamma = 1./inputParameters.rho - 1.
//    let T = inputParameters.lambda * inputParameters.time
//    let k = List.init N (fun i -> System.Numerics.Complex 0. (float(i)/float(N)*System.Math.PI*2.) )
//
//    let gVals = 
    //Gvals = G(exp(2*pi*1i*k/N),T,r,gamma);

//
//    k=0:(N-1);


