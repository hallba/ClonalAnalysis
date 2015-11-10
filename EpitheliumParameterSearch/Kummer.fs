﻿module Kummer

open MathNet.Numerics

let debug = Gamma.debug

let factorial n =
    let rec core n acc =
        if n<=0 then acc else core (n-1) (float(n)*acc)
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

let cPown (c:complex) n = 
    let rec core c n acc =
        let acc' = c * acc
        if n > 0 then core c (n-1) acc' else acc
    if n=0 then (complex 1. 0.) else core c n (complex 1. 0.)

let ci a = 
    //Convert an int to a complex number
    complex (float(a)) 0.

let cf a = 
    //Convert a float to a complex number
    complex a 0.
        

//       SUBROUTINE CHGUBI(A,B,X,HU,ID)
//C
//C       ======================================================
//C       Purpose: Compute confluent hypergeometric function
//C                U(a,b,x) with integer b ( b = ±1,±2,... )
//C       Input  : a  --- Parameter
//C                b  --- Parameter
//C                x  --- Argument
//C       Output:  HU --- U(a,b,x)
//C                ID --- Estimated number of significant digits
//C       Routines called:
//C            (1) GAMMA2 for computing gamma function Г(x)
//C            (2) PSI_SPEC for computing psi function
//C       ======================================================
//C
//        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

//Initialise key variables

//        ID=-100
//        EL=0.5772156649015329D0
//        N=ABS(B-1)
//        RN1=1.0D0
//        RN=1.0D0
//        DO 10 J=1,N
//           RN=RN*J
//           IF (J.EQ.N-1) RN1=RN
//10      CONTINUE
//        CALL PSI_SPEC(A,PS)
//        CALL GAMMA2(A,GA)

//Initialise variables dependent on sign of B

//        IF (B.GT.0.0) THEN
//           A0=A
//           A1=A-N
//           A2=A1
//           CALL GAMMA2(A1,GA1)
//           UA=(-1)**(N-1)/(RN*GA1)
//           UB=RN1/GA*X**(-N)
//        ELSE
//           A0=A+N
//           A1=A0
//           A2=A
//           CALL GAMMA2(A1,GA1)
//           UA=(-1)**(N-1)/(RN*GA)*X**N
//           UB=RN1/GA1
//        ENDIF

//Calculate hm1

//        HM1=1.0D0
//        R=1.0D0
//        HMAX=0.0D0
//        HMIN=1.0D+300
//        H0=0D0
//        DO 15 K=1,150
//           R=R*(A0+K-1.0D0)*X/((N+K)*K)
//           HM1=HM1+R
//           HU1=DABS(HM1)
//           IF (HU1.GT.HMAX) HMAX=HU1
//           IF (HU1.LT.HMIN) HMIN=HU1
//           IF (DABS(HM1-H0).LT.DABS(HM1)*1.0D-15) GO TO 20
//15         H0=HM1
//20      DA1=LOG10(HMAX)
//        DA2=0.0D0
//        IF (HMIN.NE.0.0) DA2=LOG10(HMIN)
//        ID=15-ABS(DA1-DA2)
//        HM1=HM1*DLOG(X)

//Calculate HM2

//        S0=0.0D0
//        DO 25 M=1,N
//           IF (B.GE.0.0) S0=S0-1.0D0/M
//25         IF (B.LT.0.0) S0=S0+(1.0D0-A)/(M*(A+M-1.0D0))
//        HM2=PS+2.0D0*EL+S0
//        R=1.0D0
//        HMAX=0.0D0
//        HMIN=1.0D+300
//        DO 50 K=1,150
//           S1=0.0D0
//           S2=0.0D0
//           IF (B.GT.0.0) THEN
//              DO 30 M=1,K
//30               S1=S1-(M+2.0D0*A-2.0D0)/(M*(M+A-1.0D0))
//              DO 35 M=1,N
//35               S2=S2+1.0D0/(K+M)
//           ELSE
//              DO 40 M=1,K+N
//40               S1=S1+(1.0D0-A)/(M*(M+A-1.0D0))
//              DO 45 M=1,K
//45               S2=S2+1.0D0/M
//           ENDIF
//           HW=2.0D0*EL+PS+S1-S2
//           R=R*(A0+K-1.0D0)*X/((N+K)*K)
//           HM2=HM2+R*HW
//           HU2=DABS(HM2)
//           IF (HU2.GT.HMAX) HMAX=HU2
//           IF (HU2.LT.HMIN) HMIN=HU2
//           IF (DABS((HM2-H0)/HM2).LT.1.0D-15) GO TO 55
//50         H0=HM2

//Some precision stuff- don't care

//55      DB1=LOG10(HMAX)
//        DB2=0.0D0
//        IF (HMIN.NE.0.0) DB2=LOG10(HMIN)
//        ID1=15-ABS(DB1-DB2)
//        IF (ID1.LT.ID) ID=ID1

//Calculate HM3

//        HM3=1.0D0
//        IF (N.EQ.0) HM3=0.0D0
//        R=1.0D0
//        DO 60 K=1,N-1
//           R=R*(A2+K-1.0D0)/((K-N)*K)*X
//60         HM3=HM3+R

//Calculate the output from UA, UB, HM1-3

//        SA=UA*(HM1+HM2)
//        SB=UB*HM3
//        HU=SA+SB
//        ID2=0.0D0
//        IF (SA.NE.0.0) ID1=INT(LOG10(ABS(SA)))
//        IF (HU.NE.0.0) ID2=INT(LOG10(ABS(HU)))
//        IF (SA*SB.LT.0.0) ID=ID-ABS(ID1-ID2)
//        RETURN
//        END

//Functional reading- backwards from fortran
//Only looking at HU output for now
//HU = SA+SB
//SB = UB*HM3
//SA = UA*(HM1+HM2)
//HM3 = SUM from k=1 to N-1 of f(k-1)*f(k) where f(0) = 1 and f(k) = (A2+K-1.0)/((K-N)*K)*X
//HM2 = digamma(a) + 2.0 * euler constant + S0 + SUM from k=1 to 150 of HW* f(k-1)*f(k) where f(0)=1 and f(k) = (A0+K-1.0D0)*X/((N+K)*K)
//      S0 = if B.r >= 0. then SUM from k=1 to n of -1/k
//           else SUM from k=1 to n of (1-a)/(k*(a+k-1))
//      HW(k) = 2.0 * euler constant + digamma(a)+S1-S2
//          S1(k) = if B.r >= 0. the SUM from m=1 to k of -(m + 2. * a - 2.)/(m*(m+a-1))
//                  else the SUM from m=1 to k+n of (1-a)/(m*(m+a-1))
//          S2(k) = if B.r >= 0. the SUM from m=1 to n of 1/(k+m)
//                  else the SUM from m=1 to k of 1/m
//HM1 = log(x) * SUM from k=1 to 150 of f(k-1)*f(k) where f(0) = 1 and f(k) = (A0+K-1.0D0)*X/((N+K)*K)

let seriesSum n f = 
    //Complex  numbers do not support get_zero so rewritten to be exclusive
    List.init n (fun n -> f n)
    |> List.fold (fun acc f -> acc + f) (complex 0. 0.)

let uInt a (b:complex) x =
    //Some convienience variables
    let c2 = complex 2. 0.
    let c1 = complex 1. 0.
    let c0 = complex 0. 0.
    
    //Define a tolerance for testing convergence
    let accuracy = 1.e-15
    
    //Initialise key variables
    let n = int(abs(b.r-1.))
    let cn = complex (float(n)) 0.
    let rn_1 = complex (factorial (n-1)) 0.
    let rn = if n > 0 then rn_1 * (complex (float(n)) 0.) else rn_1
    let ps = Gamma.diGammaComplex a
    let ga = Gamma.lanczosGodfrey a

    //Initialise variables dependent on the sign of b
    let (a0,a1,a2,ga1,ua,ub) = if b.r > 0. then                                                 //NB- this is some kind of reflextion
                                            let a0=a
                                            let a1= a- (complex (float(n)) 0.)
                                            let a2 = a1
                                            let ga1 = Gamma.lanczosGodfrey a1 
                                            let ua = if not (System.Double.IsInfinity(ga1.r)) then (complex -1. 0.)**(float(n)-1.)/(rn*ga1) else c0
                                            let ub = if not (System.Double.IsInfinity(ga.r)) then rn_1/ga*(x**(float(-n))) else c0
                                            (a0,a1,a2,ga1,ua,ub)
                                        else
                                            let a0 = a + cn
                                            let a1 = a0
                                            let a2 = a
                                            let ga1 = Gamma.lanczosGodfrey a1
                                            let ua = if not (System.Double.IsInfinity(ga1.r)) then (complex -1. 0.)**(float(n-1)) /(rn*ga) *(x**float(n)) else c0
                                            let ub = if not (System.Double.IsInfinity(ga.r)) then rn_1/ga1 else c0
                                            (a0,a1,a2,ga1,ua,ub)

    //printf "Initialisation complete\n"

    //Calculate hm1
    let rec calculatehm1 step tol x n a0 acc rAcc = 
        if step = 151 then acc*log(x) else
            let k = complex (float(step)) 0.
            let cn = complex (float(n)) 0.
            let rAcc' = rAcc * (a0+k-c1)*x/((cn+k)*k)
            let acc' = acc+rAcc'
            //printf "HM1 tick\n"
            //If change is below tol then the current result is good enough- return it
            if ((rAcc'.r/acc'.r)<tol) then acc*log(x)
            else calculatehm1 (step+1) tol x n a0 acc' rAcc'
    let hm1 = calculatehm1 1 accuracy x n a0 c1 c1
    //printf "HM1 complete %A\n" hm1
    //Calculate hm2
    let calculatehm2 a a0 (b:complex) ps n =
        //What happens when b=1 and n=0?
        let s0 = if b.r >= 0. then seriesSum (n) ( fun k-> -c1/( complex (float(k+1)) 0. ) ) //SUM from k=1 to n of -1/k
                 else seriesSum (n) (fun k ->   let ck = complex (float(k+1)) 0.
                                                (c1-a)/(ck*(a+ck-c1))               ) //SUM from k=1 to n of (1-a)/(k*(a+k-1))
        let calculates1 k (b:complex) a = 
            if b.r >= 0. then seriesSum (k) (fun m ->   let cm = complex (float(m+1)) 0.
                                                        -(cm+c2*a-c2)/(cm*(cm+a-c1))           )
            else seriesSum (k+n) (fun m ->      let cm = complex (float(m+1)) 0.
                                                (c1-a)/(cm*(cm+a-c1))                          )
            //if B.r >= 0. the SUM from m=1 to k of -(m + 2. * a - 2.)/(m*(m+a-1))
            //else the SUM from m=1 to k+n of (1-a)/(m*(m+a-1))
        let calculates2 k (b:complex) = 
            if b.r >= 0. then seriesSum n (fun m -> complex (1./(float(k+m+1))) 0. )
            else seriesSum k (fun m -> complex (1./(float(m+1))) 0. )
            //if B.r >= 0. the SUM from m=1 to n of 1/(k+m)
            //        else the SUM from m=1 to k of 1/m
        let rec corecalculatehm2 a a0 b ps n x step acc rAcc = 
            if step = 0 then corecalculatehm2 a a0 b ps n x 1 (ps + c2 * Gamma.eulerComplex + s0 ) c1 //initial value
            else if step = 150 then acc
            else    let k = complex (float(step)) 0.
                    let cn = complex (float(n)) 0.
                    let rAcc' = rAcc * (a0+k-c1)*x/((cn+k)*k)
                    let s1 = (calculates1 step b a)
                    let s2 = (calculates2 step b) 
                    let hw = c2 * Gamma.eulerComplex + ps + s1 - s2
                    //printf "HM2-R' %A HM2-hw' %A S1 %A S2 %A Cons %A\n" rAcc' hw s1 s2 (ps+c2 * Gamma.eulerComplex)
                    let acc' = hw * rAcc' + acc
                    //printf "acc %A diff %A acc' %A\n" acc (hw*rAcc') acc'
                    corecalculatehm2 a a0 b ps n x (step+1) acc' rAcc'
        //printf "S0 %A\n" s0
        corecalculatehm2 a a0 b ps n x 0 c0 c0
    let hm2 = calculatehm2 a a0 b ps n
    //printf "HM2 complete %A\n" hm2
    //Calculate hm3
    let rec calculatehm3 step a2 n x acc rAcc =
        if step=0 && n <> 0 then calculatehm3 (step+1) a2 n x c1 c1
        else if step = 0 then c0
        else if step = n then acc
        else 
            let k = complex (float(step)) 0.
            let cn = complex (float(n)) 0.
            let rAcc' = rAcc*(a2+k-c1)/((k-cn)*k)*x
            let acc' = acc+rAcc'
            //printf "HM3 tick\n"
            calculatehm3 (step+1) a2 n x acc' rAcc'
    let hm3 = calculatehm3 0 a2 n x c0 c0
    printf "HM1 %A HM2 %A HM3 %A\n" hm1 hm2 hm3
    printf "UA %A UB %A\n" ua ub
    //Calculate output based on hm1, hm2, hm3, ua and ub
    let sa=ua*(hm1+hm2)
    let sb=ub*hm3
    sa+sb