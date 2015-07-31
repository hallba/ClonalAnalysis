module Kummer

open MathNet.Numerics

let debug = true

let gamma = MathNet.Numerics.SpecialFunctions.Gamma

let factorial n =
    let rec core n acc =
        if n=0 then acc else core (n-1) (float(n)*acc)
    core n 1.

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
    (x/(complex 2. 0.))**(complex v 0.) / (complex (gamma(v+1.)) 0.) * (hyperGeometric0F1 (v+1.) (-x*x/(complex 4. 0.)) )

let M a b (z:complex) = 
    //Different ways to calculate the hypergeometric function, based on their strengths. 
    //Current known weaknesses
    //-Anything with complex a or b (due to limitations of gamma function)
    //--because of the U function below- which requires a gamma function which handles complex numbers
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
        let numerator' = if n = 0 then complex 1. 0. else numerator * (complex (a + float(n-1)) 0.) * z
        let denominator' = if n = 0 then complex 1. 0. else denominator * complex (float(n) * (b + float(n-1))) 0.
        let delta = numerator' / denominator'
        let result' = result + delta
        if (delta/result).Magnitude < accuracy then ignore (if debug then printf "Taylor expansion in %A steps\n" n); result' else 
            if n+1 < 500 then taylorExpansion accuracy a b z (n+1) numerator' denominator' result' else printf "a %A b %A c %A r %A\n" a b z result' ; failwith("Kummer M function failed to converge (Taylor)")
    let rec singleFraction accuracy a b z n alpha beta gamma result result' =
        //Expressing the Taylor expansion (above) as a single fraction copes better when b < 1
        let alpha' = if n = 0 then complex 0. 0. else (alpha + beta) * complex (float(n) * (b + float(n) - 1.)) 0.
        let beta'  = if n = 0 then complex 1. 0. else beta * (complex (a + float(n) - 1.) 0.) * z
        let gamma' = if n = 0 then complex 1. 0. else gamma * (complex (float(n) * (b + float(n) - 1.)) 0.)

        let result'' = (alpha' + beta')/gamma'

        if (((result'' - result')/result').Magnitude < accuracy) && (((result' - result)/result).Magnitude < accuracy) then ignore (if debug then printf "Single fraction in %A steps\n" n); result'' else
            if n+1 < 500 then singleFraction accuracy a b z (n+1) alpha' beta' gamma' result' result'' else  printf "a %A b %A c %A r %A\n" a b z result'; failwith("Kummer M function failed to converge (Fraction)")
    let rec buchholz accuracy a b z n d d' d'' result =
        //Initialise system
        //First three steps are skipped as we know the d coefficients prior to calculation
        let unchangingCoefficient = (complex (gamma(b) * 2. ** (b-1.)) 0.) * exp(z/(complex 2. 0.))
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
    if debug then printf "a %A b %A z %A\n" a b z

    match (b>=1.,sign(a)=sign(z.Real)) with 
    | (true, _) -> taylorExpansion accuracy a b z 0 (complex 1. 0.) (complex 1. 0.) (complex 0. 0.)
    | (false, _) -> singleFraction accuracy a b z 0 (complex 0. 0.) (complex 0. 0.) (complex 1. 0.) (complex 0. 0.) (complex 1. 0.)
    | (_, false)    -> buchholz accuracy a b z 0 (complex 1. 0.) (complex 0. 0.) (complex 0. 0.) (complex 0. 0.)

let U a b z = 
    //undefined for integer b, so we make small perturbations to integer
    let b = if b%1. = 0. then b + 0.00000001 else b
    (M a b z)* (complex (gamma(1.-b)/gamma(1.+a-b)) 0. ) + (M (1.+a-b) (2.-b) z) * (complex (gamma(b-1.)/gamma(a)) 0.) * z**(1.-b)
    //MathNet.Numerics.SpecialFunctions.Gamma(-2.*b) * (M a b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a - b) + MathNet.Numerics.SpecialFunctions.Gamma(2.*b) * (M a -b z) / MathNet.Numerics.SpecialFunctions.Gamma(0.5 - a + b)