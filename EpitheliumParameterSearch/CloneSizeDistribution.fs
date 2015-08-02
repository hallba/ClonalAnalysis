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
    printf "v %A w %A gg %A u %A u0 %A\n" v w gg u u0

    let Q = (complex (1. + 2.*w) 0.) - gg*u0 + ( (complex (2.*r) 0.)*(x-y) + y - (complex 1. 0.) )/(complex gamma 0.)
    
    //REALLY interesting bug. Forget to put the brackets around the  gg*u0- type checker does not find this, gives different result
    let C = ((-Q*(Whittaker.M w 0. (gg*u0)) ) + (complex (1.+2.*w) 0.)*(Whittaker.M (1.+w) 0. (gg*u0))) / (Q*(Whittaker.W w 0. (gg*u0) )+(complex 2. 0.)*(Whittaker.W (1.+w) 0. (gg*u0) ) )
    //
//Matlab
//    C = (-Q*whittakerM(w,0,gg*u0) + (1+2*w)*whittakerM(1+w,0,gg*u0)) / ...
//    (Q * whittakerW(w,0,gg*u0) + 2*whittakerW(1+w,0,gg*u0))
    printf "Q %A C %A" Q C
    //Examine C in depth
    (complex 1. 0.) - u + 
    (u*(complex (1.+v) 0.)- (complex (gamma*(1.+2.*w)) 0.) ) / complex (2.*r) 0. +
    (complex (gamma/(2.*r)) 0.) * ( (complex (1.+2.*w) 0.)* (Whittaker.M (1.+w) 0. (u*gg) ) - (complex 2. 0.)*C*(Whittaker.W (1.+w) 0. (u*gg) ) ) /
    ( (Whittaker.M w 0. (u*gg)) + C* (Whittaker.W w 0. (u*gg)) )    

let G z t r gamma =
    List.map (fun zMember -> if zMember = complex 1. 0. then complex 1. 0. else F zMember zMember t r gamma) z
//    let rec core z t r gamma acc =
//        match z with
//        | [] -> List.rev acc
//        | topZ::rest -> if topZ = complex 1. 0. then core z t r gamma ((complex 1. 0.)::acc) else core z t r gamma ((F topZ topZ t r gamma)::acc)
//    core z t r gamma []

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

let complexSum c =
    let rec core c acc = 
        match c with
        | top::rest -> core rest (acc+top)
        | [] -> acc
    core c (complex 0. 0.)
    
let probabilityCloneSizes inputParameters nRange maxN =
    let N = if List.max nRange >= maxN then List.max nRange + 1 else maxN
    let gamma = 1./inputParameters.rho - 1.
    let T = inputParameters.lambda * inputParameters.time
    let k = List.init N (fun item -> exp(complex 0. (float(item)/float(N)*System.Math.PI*2.)) )
    //CORRECT up to here- tested against matlab
    //let gVals = G k (T*1.<cell^-1>) (inputParameters.r*1.<probability^-1>) gamma
    let gVals = List.map (fun zMember -> if zMember = complex 1. 0. then complex 1. 0. else F zMember zMember (T*1.<cell^-1>) (inputParameters.r*1.<probability^-1>) gamma) k
    //printf "GVals %A\n" gVals
    let gValsSum = complexSum gVals
    gVals
    //let p = List.mapi (fun index g -> gValsSum * exp(g*(complex 0. (-2.*System.Math.PI*float(index)*float(nRange.[index]))) ) ) gVals 
    //Gvals = G(exp(2*pi*1i*k/N),T,r,gamma);


    //k=0:(N-1);
    //p

