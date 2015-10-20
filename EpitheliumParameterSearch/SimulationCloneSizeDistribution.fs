module SimulationCloneSizeDistribution

type division = {   AA  :   float;
                    BB  :   float;
                    AB  :   float;
                    B_C :   float;
                    }

type clone = {  numberA : int<Types.cell>;
                numberB : int<Types.cell>;
                numberC : int<Types.cell>;
                lambda  : float<Types.cell/Types.week>;
                rho     : float;
                r       : float;
                delta   : float;
                time    : float<Types.week>;
                rng     : System.Random;
                reportFrequency  : float<Types.week>;
                lastReport       : float<Types.week>;
                }
                with
                member this.gamma = this.lambda * this.rho / (1. - this.rho)
                member this.R = let R = ref None
                                match !R with
                                | Some(n) -> n
                                | None ->   let n = this.lambda * float(this.numberA) * 1.<Types.cell> + this.gamma * float(this.numberB) * 1.<Types.cell>
                                            R := Some(n)
                                            n
                member this.report = ()
                member this.pAA =   this.r*(1.+this.delta)*this.lambda*float(this.numberA)*1.<Types.cell>/this.R
                member this.pAB =   (1.-2.*this.r)*this.lambda*float(this.numberA)*1.<Types.cell>/this.R
                member this.pBB =   this.r*(1.-this.delta)*this.lambda*float(this.numberA)*1.<Types.cell>/this.R
                member this.pB2C =  float(this.numberB)*1.<Types.cell>*this.gamma/this.R
                member this.update =    let dt = - 1.<Types.week> * log (this.rng.NextDouble()/ (this.R*1.<Types.week Types.cell^-2> ) )
                                        let time'  = this.time + dt;
                                        let (lastReport',ignore) = if time' - this.lastReport > this.reportFrequency then (time',this.report) else (this.lastReport,())
                                        
                                        {this with time = time'; lastReport = lastReport' }

let initClone = {   numberA = 1<Types.cell>
                    numberB = 0<Types.cell>
                    numberC = 0<Types.cell>;
                    lambda  = 2.<Types.cell/Types.week>
                    rho = 0.85
                    r = 0.15
                    delta = 0.
                    time = 0.<Types.week>
                    rng = System.Random()
                    reportFrequency = 4.<Types.week>
                    lastReport = 0.<Types.week> }