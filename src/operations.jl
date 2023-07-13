function Seis.cut!(t::AbstractTraceArray, b::Union{Real,DateTime}, e::Union{Real,DateTime}; warn=true, allowempty=false)
    ib, ie, newb, isempty = Seis._cut_time_indices(t, b, e; warn, allowempty)
    t.data = isempty ? similar(Seis.trace(t), 0, length(t)) : Seis.trace(t)[ib:ie,:]
    t.b = newb
    t
end

function Seis.cut(t::AbstractTraceArray, b::Union{Real,DateTime}, e::Union{Real,DateTime}; warn=true, allowempty=false)
    ib, ie, newb, isempty = Seis._cut_time_indices(t, b, e; warn, allowempty)
    t′ = empty(t)
    t′.data = isempty ? similar(Seis.trace(t), 0, length(t)) : Seis.trace(t)[ib:ie,:]
    t′.b = newb
    t′
end

"""
    Seis.integrate!(t::AbstractTraceArray, method=:trapezium) -> t

Integrate each channel in `t`, replacing the underlying data with a copy.
"""
function Seis.integrate!(t::AbstractTraceArray, method::Symbol=:trapezium)
    if method === :trapezium
        # Note we lose one point and shift all points forward half a sample
        t.data = _integrate_trapezium(Seis.trace(t), t.delta)
        t.b += t.delta/2
    elseif method === :rectangle
        data = Seis.trace(t)
        _integrate_rectangle!(data, t.delta)
    else
        throw(ArgumentError("`method` must be one of `:trapezium` or `:rectangle`"))
    end
    t
end

function Seis.integrate(t::AbstractTraceArray, method::Symbol=:trapezium)
    t′ = empty(t)
    if method === :trapezium
        # Note we lose one point and shift all points forward half a sample
        t′.data = _integrate_trapezium(Seis.trace(t), t.delta)
        t′.b += t.delta/2
    elseif method === :rectangle
        data = copy(Seis.trace(t))
        _integrate_rectangle!(data, t.delta)
        t′.data = data
    else
        throw(ArgumentError("`method` must be one of `:trapezium` or `:rectangle`"))
    end
    t′
end    

"""
    _integrate_trapezium(m::AbstractMatrix, delta) -> ∫m

For an input matrix `m` whose columns represent some function evaluated at
evenly-spaced points `delta` apart, return an updated matrix `∫m`
whose columns contain the integrated columns of `m`.  Note therefore
that `∫m` has one less row than `m`.
"""
function _integrate_trapezium(m::AbstractMatrix, delta)
    size(m, 1) >= 2 ||
        throw(ArgumentError("need 2 or more points for trapezium integration"))
    ∫m = m[begin:end-1,:]
    h = delta/2
    @inbounds for icol in axes(∫m, 2)
        total = zero(eltype(m))
        for irow in axes(∫m, 1)
            total += h*(m[irow,icol] + m[irow+1,icol])
            ∫m[irow,icol] = total
        end
    end
    ∫m
end

"""
    _integrate_rectangle!(m::AbstractMatrix, delta) -> m

For an input matrix `m` whose columns represent some function evaluated
at evenly-spaced points `delta` apart, replace the columns in `m` with
the integrated values.
"""
function _integrate_rectangle!(m::AbstractMatrix, delta)
    @inbounds for icol in axes(m, 2)
        for irow in axes(m, 1)[2:end]
            m[irow,icol] = delta*m[irow,icol] + m[irow-1,icol]
        end
    end
    m
end

"""
    normalise!(t::AbstractTraceArray, val=1; all=false) -> t

Normalise all channels in `t` to have maximum absolute value `val`.

By default, this is done on a channel-by-channel basis.
To normalise all channels simultaneously, preserving the inter-channel
relative amplitudes, use `all=true`.  In this case, the single
maximum absolute amplitude across all the data will be normalised
to `val` instead.
"""
function Seis.normalise!(t::AbstractTraceArray, val=1; all=false)
    if all
        _normalise_all!(Seis.trace(t), val)
    else
        _normalise_columns!(Seis.trace(t), val)
    end
    t
end

"""
    _normalise_columns!(data::AbstractMatrix, val)

Normalise all columns in `data` to have the maximum absolute value
`val`.
"""
function _normalise_columns!(data::AbstractMatrix, val)
    @inbounds for icol in axes(data, 2)
        maxval = maximum(abs, @view(data[:,icol]))
        if !iszero(maxval)
            data[:,icol] .= @view(data[:,icol]) .* (val/maxval)
        end
    end
end

"""
    _normalise_all!(data, val)

Normalise all columns in `data` such that the maximum absolute value
in the whole matrix is `val`.
"""
function _normalise_all!(data::AbstractMatrix, val)
    maxval = maximum(abs, data)
    if !iszero(maxval)
        data .= data .* val/maxval
    end
    nothing
end

"""
    reverse!(t::AbstractTraceArray) -> t

Reverse the order of the channels in `t`.
"""
function Base.reverse!(t::AbstractTraceArray)
    reverse!(Seis.trace(t); dims=1)
    reverse!(t.sta)
    t
end

Base.reverse(t::AbstractTraceArray) = reverse!(deepcopy(t))
