module TraceArrays

export
    # Types
    DASArray,
    FourierDASArray,
    TraceArray,
    # Getters
    distances,
    wavenumbers,
    # IO
    read_febus,
    # Operations
    cut_distance!,
    cut_distance,
    decimate_distance!,
    decimate_distance,
    differentiate_distance!,
    differentiate_distance,
    integrate_distance,
    remove_median_trace!,
    remove_median_trace,
    # Reexported Seis things
    AbstractTrace,
    Event,
    Station,
    CartTrace,
    FourierTrace,
    cut,
    cut!,
    differentiate,
    differentiate!,
    endtime,
    fft,
    frequencies,
    ifft,
    integrate,
    integrate!,
    nfrequencies,
    normalise,
    normalise!,
    nsamples,
    remove_trend,
    remove_trend!,
    starttime,
    taper,
    taper!,
    trace

import DSP
import FFTW
import FebusTools
import Seis
# Re-exports
import Seis: AbstractTrace,
    Event,
    Station,
    CartTrace,
    FourierTrace,
    cut,
    cut!,
    differentiate,
    differentiate!,
    endtime,
    fft,
    frequencies,
    ifft,
    integrate,
    integrate!,
    nfrequencies,
    normalise,
    normalise!,
    nsamples,
    remove_trend,
    remove_trend!,
    starttime,
    taper,
    taper!,
    trace

using Dates: Dates, DateTime
using Statistics: median

# Types and Base/Seis methods
include("abstract.jl")
include("tracearray.jl")
include("dasarray.jl")
include("fourierdasarray.jl")
include("show.jl")
# Methods
include("filtering.jl")
include("operations.jl")

end
