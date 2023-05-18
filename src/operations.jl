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
        data[:,icol] .= @view(data[:,icol]) .* (val/maxval)
    end
end

"""
    _normalise_all!(data, val)

Normalise all columns in `data` such that the maximum absolute value
in the whole matrix is `val`.
"""
function _normalise_all!(data::AbstractMatrix, val)
    maxval = maximum(abs, data)
    data .= data .* val/maxval
    nothing
end
