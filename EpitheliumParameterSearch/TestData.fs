module Test

open CalculateLikelihood

let kasumiParameters =  {
                            Types.timePoints    =   [|(11.<Types.week>/7.); (21.<Types.week>/7.) ; (42.<Types.week>/7.) ; (90.<Types.week>/7.) ; (180.<Types.week>/7.) ; (365.<Types.week>/7.); (545.<Types.week>/7.)|]
                            Types.rhoRange      =   Array.init 46 (fun i -> float(i)*0.02 + 0.05 )
                            Types.rRange        =   Array.init 20 (fun i -> float(i)*0.025<Types.probability>  + 0.01<Types.probability>)
                            Types.lambdaRange   =   Array.init 50 (fun i -> float(i)*0.1<Types.cell/Types.week> + 0.1<Types.cell/Types.week>)
                            Types.maxN          =   195
                            Types.deltaRange    =   Types.Zero
                            Types.results       =   None
                            Types.excludeOnes   =   true
                            }

let kasumiCG = {
                            Types.timePoints    =   [|(11.<Types.week>/7.); (21.<Types.week>/7.) ; (42.<Types.week>/7.) ; (90.<Types.week>/7.) ; (180.<Types.week>/7.) ; (365.<Types.week>/7.); (545.<Types.week>/7.)|]
                            Types.rhoRange      =   Array.init 15 (fun i -> float(i)*0.05 + 0.15 )
                            Types.rRange        =   Array.init 9 (fun i -> float(i)*0.05<Types.probability>  + 0.1<Types.probability>)
                            Types.lambdaRange   =   Array.init 21 (fun i -> float(i)*0.2<Types.cell/Types.week> + 1.<Types.cell/Types.week>)
                            Types.maxN          =   195
                            Types.deltaRange    =   Types.Zero
                            Types.results       =   None
                            Types.excludeOnes   =   true
                            }//200 =~ 30 mins Simulation

let kasumiLocal =  {
                            Types.timePoints    =   [|(11.<Types.week>/7.); (21.<Types.week>/7.) ; (42.<Types.week>/7.) ; (90.<Types.week>/7.) ; (180.<Types.week>/7.) ; (365.<Types.week>/7.); (545.<Types.week>/7.)|]
                            Types.rhoRange      =   Array.init 18 (fun i -> float(i)*0.05 + 0.05 )
                            Types.rRange        =   Array.init 10 (fun i -> float(i)*0.025<Types.probability>  + 0.01<Types.probability>)
                            Types.lambdaRange   =   Array.init 20 (fun i -> float(i)*0.1<Types.cell/Types.week> + 0.1<Types.cell/Types.week>)
                            Types.maxN          =   195
                            Types.deltaRange    =   Types.Zero
                            Types.results       =   None
                            Types.excludeOnes   =   true
                            }

let tinyScan = {
                            Types.timePoints    =   [|(11.<Types.week>/7.); (21.<Types.week>/7.)|]
                            Types.rhoRange      =   [|0.15;0.20;0.25|]
                            Types.rRange        =   [|0.1<Types.probability>;0.15<Types.probability>|]
                            Types.lambdaRange   =   [|1.<Types.cell/Types.week>;1.2<Types.cell/Types.week>;1.4<Types.cell/Types.week>;1.6<Types.cell/Types.week>|]
                            Types.maxN          =   195
                            Types.deltaRange    =   Types.Zero
                            Types.excludeOnes   =   true
                            Types.results       =   None
                }

let testClone = Types.testSystem |> SimulationCloneSizeDistribution.parameterSetToClone (SimulationCloneSizeDistribution.Specified([1.428571429<Types.week>;])) ;;

let kasumiWTData = [
                    {
                    time=(11.<Types.week>/7.)
                    cloneSize=[|50;2|]
                    }
                    {
                    time=(11.<Types.week>/7.)
                    cloneSize=[|83;12;1;1|]
                    }
                    {
                    time=(11.<Types.week>/7.)
                    cloneSize=[|17;1|]
                    }
                    {
                    time=(11.<Types.week>/7.)
                    cloneSize=[|49;5;3|]
                    }
                    
                    {
                    time=(21.<Types.week>/7.)
                    cloneSize=[|40;8;2;1;0;0;1|]
                    }
                    {
                    time=(21.<Types.week>/7.)
                    cloneSize=[|39;11;4|]
                    }
                    {
                    time=(21.<Types.week>/7.)
                    cloneSize=[|45;7|]
                    }
                    {
                    time=(21.<Types.week>/7.)
                    cloneSize=[|30;2;1|]
                    }

                    {
                    time=(42.<Types.week>/7.)
                    cloneSize=[|40;18;1;1;1|]
                    }
                    {
                    time=(42.<Types.week>/7.)
                    cloneSize=[|32;13;1;4;2;1|]
                    }
                    {
                    time=(42.<Types.week>/7.)
                    cloneSize=[|77;35;12;1;3;0;1|]
                    }

                    {
                    time=(90.<Types.week>/7.)
                    cloneSize=[|6;3;1;1;0;0;0;1|]
                    }
                    {
                    time=(90.<Types.week>/7.)
                    cloneSize=[|4;2;2;1;0;0;1|]
                    }
                    {
                    time=(90.<Types.week>/7.)
                    cloneSize=[|22;11;1;1;3;2|]
                    }
                    {
                    time=(90.<Types.week>/7.)
                    cloneSize=[|27;14;10;3;2;1;0;1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;1|]
                    }
                    {
                    time=(90.<Types.week>/7.)
                    cloneSize=[|35;23;11;6;3;1;2;0;1|]
                    }

                    {
                    time=(180.<Types.week>/7.)
                    cloneSize=[|16;6;1;1;0;0;1;1;0;1;0;0;1|]
                    }
                    {
                    time=(180.<Types.week>/7.)
                    cloneSize=[|13;10;4;5;2;0;2;1;1;1|]
                    }
                    {
                    time=(180.<Types.week>/7.)
                    cloneSize=[|25;7;2;2;4;0;1;5;0;1;0;0;1;0;0;1;0;0;0;1;0;1|]
                    }
                    {
                    time=(180.<Types.week>/7.)
                    cloneSize=[|29;19;10;6;1;2;2;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|]
                    }
                    {
                    time=(180.<Types.week>/7.)
                    cloneSize=[|24;12;7;0;2;2;0;1;0;0;0;0;0;0;0;0;0;0;0;1|]
                    }
                    
                    {
                    time=(365.<Types.week>/7.)
                    cloneSize=[|29;17;11;9;4;0;1;0;1;0;1;2;1;1;0;0;0;0;0;0;0;0;0;0;0;0;3;1|]
                    }
                    {
                    time=(365.<Types.week>/7.)
                    cloneSize=[|22;14;9;5;6;5;0;5;1;1;1;3;1;0;0;3;1;0;0;1;0;0;0;0;0;1;0;1;0;1;0;0;0;0;0;0;0;0;0;0;1|]
                    }
                    {
                    time=(365.<Types.week>/7.)
                    cloneSize=[|25;15;9;3;5;3;5;0;0;3;2;1;0;1;1;0;2;0;0;0;0;1;0;0;0;0;0;0;0;1;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|]
                    }
                    {
                    time=(365.<Types.week>/7.)
                    cloneSize=[|28;5;5;7;1;0;0;1;1;1;1;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|]
                    }

                    {
                    time=(545.<Types.week>/7.)
                    cloneSize=[|31;10;7;5;3;5;1;2;0;0;3;1;0;0;1;1;0;1;0;0;0;0;1;0;1;0;0;0;0;1;0;0;0;0;2;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1|]
                    }
                    {
                    time=(545.<Types.week>/7.)
                    cloneSize=[|23;9;9;4;2;1;2;1;2;1;2;2;1;1;0;0;0;1;0;1;0;0;0;0;0;0;0;0;0;0;1;1;0;0;1;0;0;0;0;0;1;0;0;1;1|]
                    }
                    {
                    time=(545.<Types.week>/7.)
                    cloneSize=[|6;6;2;4;0;0;0;0;1;0;1;0;0;1;0;0;0;0;1|]
                    }
                    ]