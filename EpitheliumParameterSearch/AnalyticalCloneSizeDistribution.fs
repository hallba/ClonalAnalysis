module AnalyticalCloneSizeDistribution

open MathNet.Numerics

let F x y t r gamma = 
    //v needs to be complex for values of r over 0.25 
    //-> if v is complex, then w is complex
    //-> if w is complex then the "a" term of the Kummer function is complex
    //-> therefore we need to implement a complex gamma function in Kummer.fs
    let v = sqrt(complex (1.-4.*r) 0.)
    let w = (complex (gamma*(1.-2.*r) - 2.*r) 0.) / ( (complex (2.*gamma) 0.)*sqrt(complex (1.-4.*r) 0.))
    let gg = sqrt(complex (1.-4.*r) 0.)/(complex gamma 0.)

    let u = ((complex 1. 0.)-y)*exp(complex (-gamma*t) 0.) 
    let u0 = ((complex 1. 0.)-y )
    
    //Complex numbers shortcuts
    let c0 = complex 0. 0.
    let c1 = complex 1. 0.
    let c2 = complex 2. 0.
    let cGamma = complex gamma 0.
    let cR = complex r 0.

    //printf "v %A w %A gg %A u %A u0 %A\n" v w gg u u0

    //let Q = (complex (1. + 2.*w) 0.) - gg*u0 + ( (complex (2.*r) 0.)*(x-y) + y - (complex 1. 0.) )/(complex gamma 0.)
    let Q = c1 + c2*w - gg*u0 + ( (complex (2.*r) 0.)*(x-y) + y - c1 )/(complex gamma 0.)
    
    //Interesting bug. Forget to put the brackets around the  gg*u0- type checker does not find this, gives different result
    let C = ((-Q*(Whittaker.M w c0 (gg*u0)) ) + (c1+ (c2 * w) )*(Whittaker.M (c1+w) c0 (gg*u0))) / (Q*(Whittaker.W w c0 (gg*u0) )+c2*(Whittaker.W (c1+w) c0 (gg*u0) ) )
    //printf "Q %A C %A\n" Q C
    c1 - u + 
    (u*(v+c1)- cGamma *(c1+c2*w) ) / (c2*cR) +
    (cGamma/(c2*cR)) * ( ( c1 + c2*w ) * Whittaker.M (c1+w) c0 (u*gg)  - c2*C*(Whittaker.W (c1+w) c0 (u*gg) ) ) /
    ( (Whittaker.M w c0 (u*gg)) + C* (Whittaker.W w c0 (u*gg)) )  
    
//    1 - u + ...
//    (u*(1+v)-gamma*(1+2*w))/(2*r) + ...
//    (gamma/(2*r))*( (1+2*w)*whittakerM(1+w,0,u*gg) - 2*C*whittakerW(1+w,0,u*gg)) / ...
//    (whittakerM(w,0,u*gg) + C*whittakerW(w,0,u*gg));  

let probabilityCloneSurvival (inputParameters: Types.parameterSet) =
    //Note that the below line is incorrect, and the effects corrected later
    //let gamma = 1./inputParameters.rho - 1.
    //Corrected version. Note that the function will need to be inverted for comparison with equivalent matlab code
    let gamma = inputParameters.rho/(1.-inputParameters.rho)
    let T = inputParameters.lambda * inputParameters.time
    let p = (complex 1. 0.) - F (complex 0. 0.) (complex 0. 0.) (T*1.<Types.cell^-1>) (inputParameters.r*1.<Types.probability^-1>) gamma
    //printf "P: %A\n" p
    //The following line checks for "pathological points" as described in the original matlab implementation. Note that, following the matlab code, this is a different transformation than before
    if p.r < 0. || p.r > 1. || System.Double.IsNaN(p.r) || System.Double.IsInfinity(p.r) then ((complex 1. 0.) - F (complex 0. 0.) (complex 0. 0.) (T*1.<Types.cell^-1>) ((inputParameters.r*1.<Types.probability^-1>)+0.00001) (gamma+0.00001)) else p

let complexSum c =
    let rec core c acc = 
        match c with
        | top::rest -> core rest (acc+top)
        | [] -> acc
    core c (complex 0. 0.)

let probabilityCloneSizes (inputParameters : Types.parameterSet) nRange maxN =
    let rec core (inputParameters :Types.parameterSet) nRange maxN attempt maxAttempts =
        //No infinite loops
        if attempt = maxAttempts then failwith("Too many attempts to modify a pathological point")
        
        let N = if List.max nRange >= maxN then List.max nRange + 1 else maxN
        //Note that the below line is incorrect, and the effects corrected later
        //let gamma = 1./inputParameters.rho - 1.
        //Corrected version. Note that the function will need to be inverted for comparison with equivalent matlab code
        let gamma = inputParameters.rho/(1.-inputParameters.rho)
        let T = inputParameters.lambda * inputParameters.time
        //Some bug in k calculation
        let k = List.init N (fun item -> exp(complex 0. (float(item)/float(N)*System.Math.PI*2.)) )
        let gVals = List.map (fun zMember -> if zMember = complex 1. 0. then complex 1. 0. else F zMember zMember (T*1.<Types.cell^-1>) (inputParameters.r*1.<Types.probability^-1>) gamma) k 
        let p = List.map (fun n -> List.mapi (fun index i -> i * exp(complex 0. (-2.*System.Math.PI*float(n)*float(index)/float(N)) ) ) gVals ) nRange 
                |>List.map (fun item -> (complexSum item).r/float(N) )
                |>List.map (fun item -> if item < 0. then 0. else item) //In the matlab code we had an additional test which could allow values less
                                                                        //than -1e-5 to be included in the output. Perhaps include in future?
        let totalP = List.sum p
        //Test for pathological points
        if totalP > 1. || totalP < 0. || System.Double.IsNaN(totalP) then 
            //NB: this is from the matlab code- note that it is a specified small perturbation, not 0.1%
            ignore (printf "Pathological point: %A.\nsum%A=%A\n Making 0.1%% perturbation\n" inputParameters p totalP); core inputParameters.correctPathologicalPoint nRange maxN (attempt+1) maxAttempts else
            p |> Array.ofList
    core inputParameters nRange maxN 0 1
