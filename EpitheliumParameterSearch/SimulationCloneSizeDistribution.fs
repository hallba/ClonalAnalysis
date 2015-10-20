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
                time    : float<Types.week>;
                rng     : System.Random;
                reportFrequency  : float<Types.week>;
                lastReport       : float<Types.week>}
                with
                member this.gamma = this.lambda * this.rho / (1. - this.rho)
                member this.R = this.lambda * float(this.numberA) * 1.<Types.cell> + this.gamma * float(this.numberB) * 1.<Types.cell>
                member this.report = ()
                member this.update =    let dt = - 1.<Types.week> * log (this.rng.NextDouble()/ (this.R*1.<Types.week Types.cell^-2> ) )
                                        let time'  = this.time + dt;
                                        let (lastReport',ignore) = if time' - this.lastReport > this.reportFrequency then (time',this.report) else (this.lastReport,())
                                        
                                        {this with time = time'; lastReport = lastReport' }