module Test

let kasumiParameters =  {
                            Types.timePoints    =   [|(10.<Types.week>/7.); (11.<Types.week>/7.); (21.<Types.week>/7.) ; (42.<Types.week>/7.) ; (90.<Types.week>/7.) ; (180.<Types.week>/7.) ; (365.<Types.week>/7.); (545.<Types.week>/7.)|]
                            Types.rhoRange      =   Array.init 46 (fun i -> float(i)*0.02 + 0.05 )
                            Types.rRange        =   Array.init 20 (fun i -> float(i)*0.025<Types.probability>  + 0.01<Types.probability>)
                            Types.lambdaRange   =   Array.init 50 (fun i -> float(i)*0.1<Types.cell/Types.week> + 0.1<Types.cell/Types.week>)
                            Types.maxN          =   195
                            Types.deltaRange    =   Types.Zero
                            Types.results       =   None
                            }

let kasumiCG = {
                            Types.timePoints    =   [|(10.<Types.week>/7.); (11.<Types.week>/7.); (21.<Types.week>/7.) ; (42.<Types.week>/7.) ; (90.<Types.week>/7.) ; (180.<Types.week>/7.) ; (365.<Types.week>/7.); (545.<Types.week>/7.)|]
                            Types.rhoRange      =   Array.init 10 (fun i -> float(i)*0.1 + 0.05 )
                            Types.rRange        =   Array.init 4 (fun i -> float(i)*0.1<Types.probability>  + 0.05<Types.probability>)
                            Types.lambdaRange   =   Array.init 5 (fun i -> float(i)*1.0<Types.cell/Types.week> + 0.1<Types.cell/Types.week>)
                            Types.maxN          =   195
                            Types.deltaRange    =   Types.Zero
                            Types.results       =   None
                            }