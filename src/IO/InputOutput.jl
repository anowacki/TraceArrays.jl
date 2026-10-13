"""
    InputOutput

Module for reading and writing data in array-based formats, including
DAS data
"""
module InputOutput

import Dates
import HDF5
import Seis
import TraceArrays

export
    read_optdas

include("optodas.jl")

end # module
