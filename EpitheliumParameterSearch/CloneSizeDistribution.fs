module CloneSizeDistribution

[<Measure>]
type cell

[<Measure>]
type week

[<Measure>]
type probability

let F x y t r gamma = 
    let v = sqrt(1.-4.*r)
    let w = (gamma*(1.-2.*r) - 2.*r) / (2.*gamma*sqrt(1.-4.*r))
    let gg = sqrt(1.-4.*r)/gamma

    let u = (1.-y)*exp(-gamma*t)
    let u0 = (1.-y)

    let Q = 1. + 2.*w - gg*u0 + (2.*r*(x-y)+y-1.)/gamma

    let C = ((-Q*(Whittaker.M w 0. gg*u0) ) + (1.+2.*w)*(Whittaker.M (1.+w) 0. (gg*u0))) / (Q*(Whittaker.W w 0. gg*u0 )+2.*(Whittaker.W (1.+w) 0. (gg*u0) ) )
    1. - u + 
    (u*(1.+v)-gamma*(1.+2.*w))/(2.*r) +
    (gamma/(2.*r))*( (1.+2.*w)*(Whittaker.M (1.+w) 0. (u*gg)) - 2.*C*(Whittaker.W (1.+w) 0. (u*gg) ) ) /
    ( (Whittaker.M w 0. (u*gg)) + C* (Whittaker.W w 0. (u*gg)) )    

let G z t r gamma =
    let rec core z t r gamma acc =
        match z with
        | [] -> acc
        | topZ::rest -> if topZ = 1. then core z t r gamma (1.::acc) else core z t r gamma ((F topZ topZ t r gamma)::acc)
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

let correctPathologicalPoints inputParameters =
    { inputParameters with rho = inputParameters.rho + 0.00001 }

let probabilityClonalSurvival inputParameters =
    //Note that the below line is incorrect, and the effects corrected later
    let gamma = 1./inputParameters.rho - 1.
    let T = inputParameters.lambda * inputParameters.time
    let p = 1. - F 0. 0. (T*1.<cell^-1>) (inputParameters.r*1.<probability^-1>) gamma
    //printf "P: %A\n" p
    //The following line checks for "pathological points" as described in the original matlab implementation
    if p < 0. || p > 1. || System.Double.IsNaN(p) || System.Double.IsInfinity(p) then (1. - F 0. 0. (T*1.<cell^-1>) ((inputParameters.r*1.<probability^-1>)+0.00001) (gamma+0.00001)) else p
    
let probabilityCloneSizes inputParameters nRange nIter maxN =
    let N = if List.max nRange > nIter then List.max nRange + 1 else nIter
    let gamma = 1./inputParameters.rho - 1.
    let T = inputParameters.lambda * inputParameters.time
    let k = List.init N (fun i -> System.Numerics.Complex 0. (float(i)/float(N)*System.Math.PI*2.) )

    let gVals = 
    //Gvals = G(exp(2*pi*1i*k/N),T,r,gamma);

//
//    k=0:(N-1);


