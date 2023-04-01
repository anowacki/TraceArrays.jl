function Seis.cut!(t::AbstractTraceArray, b, e; warn=true, allowempty=false)
    ib, ie = _cut_time_indices(t, b, e; warn, allowempty)
    t.data = Seis.trace(t)[ib:ie,:]
    t.b = ib in eachindex(Seis.times(t)) ? Seis.times(t)[ib] : b
    t
end

function Seis.cut(t::AbstractTraceArray, b, e; warn=true, allowempty=false)
    ib, ie = _cut_time_indices(t, b, e; warn, allowempty)
    t′ = empty(t)
    t′.data = Seis.trace(t)[ib:ie,:]
    t′.b = ib in eachindex(Seis.times(t′)) ? Seis.times(t′)[ib] : b
    t′
end

"""
    _cut_time_indices(t::AbstractTraceArray, b, e; warn, allowempty) -> ib, ie

Compute the start index `ib` and end index `ie` which cuts the trace array
`t` in time between `b` s and `e` s.
"""
function _cut_time_indices(t::AbstractTraceArray, b, e; warn=true, allowempty=false)
    (b === missing || e === missing) && throw(ArgumentError("Start or end cut time is `missing`"))
    e < b && throw(ArgumentError("End cut time ($e) is before start cut ($b)"))
    if b > Seis.endtime(t) || e < Seis.starttime(t)
        if !allowempty
            b > Seis.endtime(t) &&
                throw(ArgumentError("Beginning cut time $b is later than end of trace ($(Seis.endtime(t)))."))
            e < Seis.starttime(t) &&
                throw(ArgumentError("End cut time $e is earlier than start of trace (t.b)."))
        end
        empty!(t.t)
        t.b = b
        return t
    end
    if b < t.b
        warn && @warn("Beginning cut time $b is before start of trace.  Setting to $(t.b).")
        b = t.b
    end
    if e > Seis.endtime(t)
        warn && @warn("End cut time $e is after end of trace.  Setting to $(Seis.endtime(t)).")
        e = Seis.endtime(t)
    end
    ib = round(Int, (b - t.b)/t.delta) + 1
    ie = Seis.nsamples(t) - round(Int, (Seis.endtime(t) - e)/t.delta)
    ib, ie
end
