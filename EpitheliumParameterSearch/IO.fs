module IO

open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.Data.Matlab

let storeSearch (filename: string) (search: Types.parameterSearch) = 
        use file = new System.IO.FileStream(filename, System.IO.FileMode.Create)
        let bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter()
        use mem = new System.IO.MemoryStream()
        bf.Serialize(mem,search)
        let byteArray = mem.ToArray()
        file.Write(byteArray,0,byteArray.Length)
        file.Close()
        ()

let retrieveSearch (filename: string) =
    use fileStream = new System.IO.FileStream(filename, System.IO.FileMode.Open)
    let bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter()
    let result = bf.Deserialize(fileStream)
    fileStream.Close()
    (result :?> Types.parameterSearch)

//Get data from a matlab generated matrix
let importMatlab (filename: string) = 
    //use fileStream = new System.IO.FileStream(filename, System.IO.FileMode.Open)
    let m = (MatlabReader.ReadAll<float>(filename,"rhoRange","lambdaRange","rRange","PScanPP","PSurvScanPP","t","tExp","timePoints") ) // :?> System.Collections.Generic.Dictionary< string, Matrix<float> > 
    //fileStream.Close()
    m