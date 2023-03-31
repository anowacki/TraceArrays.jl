module TraceArrays

export
    # Types
    DASArray,
    TraceArray,
    # Getters
    distances,
    # IO
    read_febus,
    # Operations
    cut_distance!,
    cut_distance,
    remove_median_trace!,
    remove_median_trace

import FebusTools
import Seis

using Dates: DateTime
using Seis: AbstractTrace, Event, Station
using Statistics: median

include("abstract.jl")
include("tracearray.jl")
include("dasarray.jl")


end
