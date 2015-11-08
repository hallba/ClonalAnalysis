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

//let U' a (b:complex) z =
//    let rec core a (b:complex) z attempt =
//        if attempt > 10 then failwith "stuck in a loop"
//        printf "b=%A\n" b
//        if b.r%1. <> 0. || b.r > 1. then
//                let result = (M a b z)* exp( Gamma.logLanczosGodfrey ((complex 1. 0.)-b) - Gamma.logLanczosGodfrey ((complex 1. 0.)+a-b) ) + (M ((complex 1. 0.)+a-b) ((complex 2. 0.)-b) z) * exp( Gamma.logLanczosGodfrey (b-(complex 1. 0.)) - Gamma.logLanczosGodfrey a ) * z**((complex 1. 0.)-b)
//                if debug then printf "U(%A,%A,%A)\nResult=%A\n" a b z result
//                result 
//        else
//            let c1 = complex 1. 0.
//            let c2 = complex 2. 0.
//            let a' = a+c1-b
//            let b' = c2-b
//            printf "b'=%A\n" b'
//            let result = z**(c1-b) * (core a' b' z (attempt+1) )
//            if debug then printf "U(%A,%A,%A)\nResult=%A\n" a b z result
//            result
//    core a b z 0

let cPown (c:complex) n = 
    let rec core c n acc =
        let acc' = c * acc
        if n > 0 then core c (n-1) acc' else acc
    if n=0 then (complex 1. 0.) else core c n (complex 1. 0.)

let ci a = 
    //Convert a float or an int to a complex number
    complex (float(a)) 0.

let cf a = 
    //Convert a float or an int to a complex number
    complex a 0.

let UInt a (b:complex) (x:complex) =
    let accuracy = 1e-15
    let n = int(abs(b.r-1.))
    let cn = complex (float(n)) 0.
    let c2 = complex 2. 0.
    let c1 = complex 1. 0.
    let c0 = complex 0. 0.
    let rn_1 = complex (factorial (n-1)) 0.
    let rn = rn_1 * (complex (float(n)) 0.)
    let ps = Gamma.diGammaComplex a
    let ga = Gamma.lanczosGodfrey a
    let (a0,a1,a2,ga1,ua,ub) = if b.r > 0. then                                                 //NB- this is some kind of reflextion
                                               let a0=a
                                               let a1= a- (complex (float(n)) 0.)
                                               let a2 = a1
                                               let ga1 = Gamma.lanczosGodfrey a1
                                               let ua = (complex -1. 0.)**(float(n)-1.)/(rn*ga1)
                                               //UA=(-1)**(N-1)/(RN*GA1)
                                               let ub = rn_1/ga*(x**(float(-n)))
                                               //UB=RN1/GA*X**(-N)
                                               (a0,a1,a2,ga1,ua,ub)
                                            else
                                                //           A0=A+N
                                                let a0 = a + cn
                                                //           A1=A0
                                                let a1 = a0
                                                //           A2=A
                                                let a2 = a
                                                //           CALL GAMMA2(A1,GA1)
                                                let ga1 = Gamma.lanczosGodfrey a1
                                                //           UA=(-1)**(N-1)/(RN*GA)*X**N
                                                let ua = (complex -1. 0.)**(float(n-1)) /(rn*ga) *(x**float(n))
                                                //           UB=RN1/GA1
                                                let ub = rn_1/ga1
                                                (a0,a1,a2,ga1,ua,ub)
    let hMax,hMin,hM1,h0 =  List.init 150 (fun k -> let ck = complex (float(k)) 0.
                                                    cPown ((a0+ck-c1)*x/((ck+cn)*ck)) (k+1)        )
                            |> List.fold (fun ((hMax:complex),(hMin:complex),hM1,h0) r ->   let hM1' = hM1 + r
                                                                                            let hu1 = complex (abs(hM1'.r)) hM1'.i
                                                                                            let hMax',hMin' =   match (hu1.r>hMax.r,hu1.r<hMin.r) with
                                                                                                                | (false,false) -> hMax, hMin
                                                                                                                | (true,false)  -> hu1, hMin
                                                                                                                | (false,true)  -> hMax, hu1
                                                                                                                | (true,true)   -> hu1, hu1
                                                                                            //Some sort of break in the loop missing here based on accuracy
                                                                                            let h0' = hM1'
                                                                                            (hMax',hMin',hM1',h0')
                                                                                        ) (c0,(complex 1.e300 0.),c1,c0)

                                                
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
    let da1 = log10(hMax)
    let da2 = if hMin.r<>0. then log10(hMin) else c0
    let id = 15. - abs(da1.r-da2.r)
    let hM1 = hM1*log(x)
    let s0 = List.init (n-1) (fun m ->  let cm = complex (float(m+1)) 0.
                                        if b.r >= 0. then -c1/cm else (c1-a)/(cm*(a+cm-c1) ) )
             |> List.fold (fun acc ms -> acc + ms) c0
//20      DA1=LOG10(HMAX)
//        DA2=0.0D0
//        IF (HMIN.NE.0.0) DA2=LOG10(HMIN)
//        ID=15-ABS(DA1-DA2)
//        HM1=HM1*DLOG(X)
//        S0=0.0D0
//        DO 25 M=1,N
//           IF (B.GE.0.0) S0=S0-1.0D0/M
//25         IF (B.LT.0.0) S0=S0+(1.0D0-A)/(M*(A+M-1.0D0))
    let hM2 = ps + (complex 2. 0.)*Gamma.eulerComplex+s0
//        HM2=PS+2.0D0*EL+S0
    let calculateS1S2 vk (vb:complex) va =    
        List.init (vk-1) (fun k -> k+1)
        |> List.fold (fun (s1,s2) vm -> if vb.r > 0. then
                                            ( (s1-((ci vm)+c2*va-c2)) , (s2+c1/((ci vk)+(ci vm) )) ) 
                                        else 
                                            ( (s1+(c1-va)/((ci vm)*((ci vm)+va-c1))) , (s2+c1/(ci vm) )             )       ) (c0,c0)
    
    let (hM2,hu2) =     List.init 149 (fun k -> calculateS1S2 (k+1) b a )
                        |> List.fold () (hM2,hu2)
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
    let db1 = log10 hMax
    let db2 = if hMin.r <> 0. then log10 hMin else c0
    let id1 = 15. - abs((db1-db2).r)
    let id = if id1 < -100. then (complex id1 (db1-db2).i ) else (complex -100. 0.)
    //let hM3 = (if n=0 then c0 else c1)
    let hM3 =   List.init (n-1) (fun k ->   let ck = complex (float(k)) 0.
                                            let cn = complex (float(n)) 0.
                                            cPown ((a2+ck-c1)/((ck-cn)*ck)*x) k 
                                                                                )
                |> List.fold (fun acc r -> acc + r ) (if n=0 then c0 else c1)
    let sa = ua*(hM1+hM2)
//55      DB1=LOG10(HMAX)
//        DB2=0.0D0
//        IF (HMIN.NE.0.0) DB2=LOG10(HMIN)
//        ID1=15-ABS(DB1-DB2)
//        IF (ID1.LT.ID) ID=ID1
//        HM3=1.0D0
//        IF (N.EQ.0) HM3=0.0D0
//        R=1.0D0
//        DO 60 K=1,N-1
//           R=R*(A2+K-1.0D0)/((K-N)*K)*X
//60         HM3=HM3+R
//        SA=UA*(HM1+HM2)
//        SB=UB*HM3
//        HU=SA+SB
//        ID2=0.0D0
//        IF (SA.NE.0.0) ID1=INT(LOG10(ABS(SA)))
//        IF (HU.NE.0.0) ID2=INT(LOG10(ABS(HU)))
//        IF (SA*SB.LT.0.0) ID=ID-ABS(ID1-ID2)


    if false then complex 1. 0. else failwith "Function not complete"

        
        

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
//55      DB1=LOG10(HMAX)
//        DB2=0.0D0
//        IF (HMIN.NE.0.0) DB2=LOG10(HMIN)
//        ID1=15-ABS(DB1-DB2)
//        IF (ID1.LT.ID) ID=ID1
//        HM3=1.0D0
//        IF (N.EQ.0) HM3=0.0D0
//        R=1.0D0
//        DO 60 K=1,N-1
//           R=R*(A2+K-1.0D0)/((K-N)*K)*X
//60         HM3=HM3+R
//        SA=UA*(HM1+HM2)
//        SB=UB*HM3
//        HU=SA+SB
//        ID2=0.0D0
//        IF (SA.NE.0.0) ID1=INT(LOG10(ABS(SA)))
//        IF (HU.NE.0.0) ID2=INT(LOG10(ABS(HU)))
//        IF (SA*SB.LT.0.0) ID=ID-ABS(ID1-ID2)
//        RETURN
//        END