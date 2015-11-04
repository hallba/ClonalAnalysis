module Kummer

open MathNet.Numerics

let debug = Gamma.debug

let factorial n =
    let rec core n acc =
        if n=0 then acc else core (n-1) (float(n)*acc)
    core n 1.

let cGamma = Gamma.complexGamma Gamma.lanczosGodfrey

let hyperGeometric0F1 a z =
    //The confluent hypergeometric limit function
    //known weakness- large z
    let accuracy = pown 10. -15
    let rec core a z n delta result =
        let delta' = if n=0 then complex 1. 0. else delta * z / (complex ((float(n)+a-1.)*(float(n))) 0. )
        let result' = if n=0 then complex 1. 0. else result + delta'
        match( ((delta'/result').Magnitude < accuracy), n<500) with
        | (true,_) -> result'
        | (false,true) -> core a z (n+1) delta' result'
        | (false,false) -> failwith("The confluent hypergeometric limit function failed to converge")
    core a z 0 (complex 0. 0.) (complex 0. 0.)

let besselJ v x =
    (x/(complex 2. 0.))**(complex v 0.) / (complex (Gamma.realGamma(v+1.)) 0.) * (hyperGeometric0F1 (v+1.) (-x*x/(complex 4. 0.)) )

let M (a:complex) (b:complex) (z:complex) = 
    //Test the input for anything untoward- we cannot cope with NaN
    if System.Double.IsNaN(a.Magnitude) || System.Double.IsNaN(b.Magnitude) || System.Double.IsNaN(z.Magnitude) then failwith "Kummer M NaN input: a %A b %A z %A\n" a b z
    
    //Different ways to calculate the hypergeometric function, based on their strengths. 
    //Current known weaknesses
    //-z > 100
    //-a or b > 50
    //-sign(Re(a)) = negative sign(Re(z))
    //
    //Things to implement in future
    //-Recurrance relation ((z> 100 or (a>50 or b>50))
    //-Asymptotic (z>100 and a<50 and b<50)
    //-Bucholz polynomial (sign(a) = negative sign(z), a<50, b<50, z<100
    //
    //See Masters Thesis of John Pearson for more details https://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf
    let accuracy = pown 10. -15
    let rec taylorExpansion accuracy a b z n numerator denominator result =
        //Taylor expansion of M. This has been reported as an efficent way of calculating M, but with known weaknesses
        //let numerator' = if n = 0 then complex 1. 0. else numerator * (complex (a + float(n-1)) 0.) * z
        //let denominator' = if n = 0 then complex 1. 0. else denominator * complex (float(n) * (b + float(n-1))) 0.
        let numerator' = if n = 0 then complex 1. 0. else numerator * (a + (complex (float(n-1)) 0.) ) * z
        let denominator' = if n = 0 then complex 1. 0. else denominator * (complex (float(n)) 0.)*(b + (complex (float(n-1)) 0.) )
        let delta = numerator' / denominator'
        let result' = result + delta
        if (delta/result).Magnitude < accuracy then ignore (if debug then printf "Taylor expansion in %A steps\na=%A b=%A c=%A result=%A\nhypergeom(%A+%Ai,%A+%Ai,%A+%Ai)\n######\n" n a b z result' a.r a.i b.r b.i z.r z.i); result' else 
            if n+1 < 500 then taylorExpansion accuracy a b z (n+1) numerator' denominator' result' else printf "a %A b %A c %A r %A\n" a b z result' ; failwith("Kummer M function failed to converge (Taylor)")
    let rec singleFraction accuracy a b z n alpha beta gamma result result' =
        //Expressing the Taylor expansion (above) as a single fraction copes better when b < 1
        //let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * complex (float(n) * (b + float(n) - 1.)) 0.
        let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * (complex (float(n)) 0.) * (b + (complex (float(n)-1.) 0.) )  //complex (float(n) * (b + float(n) - 1.)) 0.
        let beta'  = if n = 0 then complex 1. 0. else beta * (a + complex (float(n) - 1.) 0.) * z
        let gamma' = if n = 0 then complex 1. 0. else gamma * (complex (float(n)) 0.) * (b + (complex (float(n) - 1.) 0.))

        let result'' = (alpha' + beta')/gamma'
        //printf "alpha %A beta %A gamma %A result %A\nalpha' %A beta' %A gamma' %A result'%A\n" alpha beta gamma result alpha' beta' gamma' result'
        
        //result can diverge for complex z so fallback to taylor
        if System.Double.IsNaN(result''.r) || System.Double.IsNaN(result''.i) 
            then taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.)
        else if (((result'' - result')/result').Magnitude < accuracy) && (((result' - result)/result).Magnitude < accuracy) then ignore (if debug then printf "Single fraction in %A steps\na=%A b=%A c=%A result=%A\nhypergeom(%A+%Ai,%A+%Ai,%A+%Ai)\n" n a b z result'' a.r a.i b.r b.i z.r z.i); result'' else
            if n+1 < 500 then singleFraction accuracy a b z (n+1) alpha' beta' gamma' result' result'' else  printf "a %A b %A c %A r %A\n" a b z result'; failwith("Kummer M function failed to converge (Fraction)")
    let rec buchholz accuracy a b z n d d' d'' result =
        //Initialise system
        //First three steps are skipped as we know the d coefficients prior to calculation
        let unchangingCoefficient = (complex (Gamma.realGamma(b) * 2. ** (b-1.)) 0.) * exp(z/(complex 2. 0.))
        let sqrt_z_x_4a_minus_2b = sqrt( z*(complex (4.*a-2.*b) 0.) )
        let result = if n = 0 then unchangingCoefficient * ( ( besselJ (b-1.)    ( sqrt_z_x_4a_minus_2b ) )/(sqrt_z_x_4a_minus_2b**( b-1.   )) +
                                                             ( besselJ (b-1.+2.) ( sqrt_z_x_4a_minus_2b ) )/(sqrt_z_x_4a_minus_2b**( b-1.+2. )) * (complex (b/2.) 0.) * z ** 2.  ) else result
        let n = if n = 0 then 3 else n
        let d''' = ( (complex (float(n) - 2. + b) 0.) * d' + (complex (2.*a - b) 0.) * d )/(complex (float(n)) 0.)
        let delta = unchangingCoefficient * d''' * (besselJ (b-1.+float(n)) (sqrt_z_x_4a_minus_2b ) ) / (sqrt_z_x_4a_minus_2b**( (b-1.+float(n)) ))
        let result' = result + delta
        match (((delta/result').Magnitude < accuracy) ,n<500) with 
        | (true,_) -> ignore (if debug then printf "Buchholz in %A steps\n" n); result'
        | (false,true) -> buchholz accuracy a b z (n+1) d' d'' d''' result'
        | (false,false) -> printf "a %A b %A c %A\n" a b z ; failwith("Kummer M function failed to converge (Buchholz)")

    
    //Decide which approach to use
    //Temporarily disabling Buchholz method until I've confirmed that I can reproduce the original matlab results
    //Need a test which copes with NaN sign(a)=sign(z.Real)
    //if debug then printf "a %A b %A z %A\n" a b z

    match (b.Magnitude>=1.,false) with 
    | (true, _) -> taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.)
    | (false, _) -> singleFraction accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 1. 0.) (complex 0. 0.) (complex 1. 0.)
    //| (_, false)    -> buchholz accuracy a b z 0 (complex 1. 0.) (complex 0. 0.) (complex 0. 0.) (complex 0. 0.)

let U a (b:complex) z = 
    //undefined for integer b, so we make small perturbations to integer
    let b = if b.r%1. = 0. then b + complex 0.000000000001 0. else b
    //(M a b z)* (complex (gamma(1.-b)/gamma(1.+a-b)) 0. ) + (M (1.+a-b) (2.-b) z) * (complex (gamma(b-1.)/gamma(a)) 0.) * z**(1.-b)
    //(M a b z)* ( cGamma ((complex 1. 0.)-b)/cGamma ((complex 1. 0.)+a-b) ) + (M ((complex 1. 0.)+a-b) ((complex 2. 0.)-b) z) * ( cGamma (b-(complex 1. 0.))/ cGamma a ) * z**((complex 1. 0.)-b)
    
    let firstTermGammaRatio =   exp( Gamma.logLanczosGodfrey ((complex 1. 0.)-b) - Gamma.logLanczosGodfrey ((complex 1. 0.)+a-b) )
                                |> fun i -> if System.Double.IsNaN(i.r) then (Gamma.lanczosGodfrey ((complex 1. 0.)-b))/(Gamma.lanczosGodfrey ((complex 1. 0.)+a-b)) else i
    let secondTermGammaRatio =  exp( Gamma.logLanczosGodfrey (b-(complex 1. 0.)) - Gamma.logLanczosGodfrey a )
                                |> fun i -> if System.Double.IsNaN(i.r) then (Gamma.lanczosGodfrey ((complex 1. 0.)-b))/(Gamma.lanczosGodfrey ((complex 1. 0.)+a-b)) else i
    let result = (M a b z)* firstTermGammaRatio + (M ((complex 1. 0.)+a-b) ((complex 2. 0.)-b) z) * z**((complex 1. 0.)-b) * secondTermGammaRatio
    if debug then printf "U(%A,%A,%A)\nResult=%A\n" a b z result
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*b) * (M a b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a - b) + MathNet.Numerics.SpecialFunctions.Gamma(2.*b) * (M a -b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a + b)
    result 

let U' a (b:complex) z =
    let rec core a (b:complex) z attempt =
        if attempt > 10 then failwith "stuck in a loop"
        printf "b=%A\n" b
        if b.r%1. <> 0. || b.r > 1. then
                let result = (M a b z)* exp( Gamma.logLanczosGodfrey ((complex 1. 0.)-b) - Gamma.logLanczosGodfrey ((complex 1. 0.)+a-b) ) + (M ((complex 1. 0.)+a-b) ((complex 2. 0.)-b) z) * exp( Gamma.logLanczosGodfrey (b-(complex 1. 0.)) - Gamma.logLanczosGodfrey a ) * z**((complex 1. 0.)-b)
                if debug then printf "U(%A,%A,%A)\nResult=%A\n" a b z result
                result 
        else
            let c1 = complex 1. 0.
            let c2 = complex 2. 0.
            let a' = a+c1-b
            let b' = c2-b
            printf "b'=%A\n" b'
            let result = z**(c1-b) * (core a' b' z (attempt+1) )
            if debug then printf "U(%A,%A,%A)\nResult=%A\n" a b z result
            result
    core a b z 0
