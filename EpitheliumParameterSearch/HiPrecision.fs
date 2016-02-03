module HiPrecision

open MathNet.Numerics

let a = 4N

type System.Decimal with
    static member accuracy=1e-30m
    static member log d = failwith ""
    static member exp d = failwith ""

type complexDecimal = {     real:       decimal
                            imaginary:  decimal
                            } with
                        member this.r = this.real
                        member this.i = this.imaginary
                        member this.Magnitude = (this.i*this.i + this.r*this.r)
                        static member (+) (c1,c2) = {c1 with real=(c1.real+c2.real);imaginary=(c1.imaginary+c2.imaginary)}
                        static member (-) (c1,c2) = {c1 with real=(c1.real-c2.real);imaginary=(c1.imaginary-c2.imaginary)}
                        //static member (-) (c1)    = {c1 with real= -c1.real;imaginary= -c1.imaginary}
                        
                        static member (*) (c1,c2) = {c1 with real=(c1.real*c2.real-c1.imaginary*c2.imaginary);imaginary=(c1.imaginary*c2.real+c2.imaginary*c1.real) }
                        static member (/) (c1,c2) = let complexConjugate = (c2.real*c2.real)+(c2.imaginary*c2.imaginary)
                                                    {c1 with real=(c1.real*c2.real+c1.imaginary*c2.imaginary)/complexConjugate;imaginary=(c1.imaginary*c2.real-c1.real-c2.imaginary)/complexConjugate}
                        
                        static member log c1  = failwith ""  
                        static member exp c1  = failwith ""

                        static member (/**) (c1,c2) = failwith ""
                        static member sqrt c1 = 

let complex (r:float) (i:float) = 
    {real= decimal(r); imaginary=decimal(i)}