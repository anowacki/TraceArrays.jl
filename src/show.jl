function Base.show(io::IO, ::MIME"text/plain", t::T) where {T<:AbstractTraceArray}
    indent = 2
    hdr_string_len = maximum(length.(string.(fieldnames(T)))) + indent

    function padded_print(name, vals...)
        print(io, "\n", lpad(string(name), hdr_string_len + indent) * ": ", vals...)
    end

    print(io, Seis.nsamples(t), " samples ⨯ ", length(t), " channels ", T, ":")
    for field in (:b, :delta)
        padded_print(field, getproperty(t, field))
    end

    print(io, "\n", " "^(indent÷2), typeof(t.evt), ":")
    for field in propertynames(t.evt)
        field === :pos && continue
        value = getproperty(t.evt, field)
        if !ismissing(value)
            field_name = string("evt.", field)
            if value isa Dict
                isempty(value) && continue
                padded_print(field_name, "")
                kindex = 0
                for (k, v) in value
                    kindex += 1
                    if kindex > 1
                        print(io, "\n", " "^(hdr_string_len + indent + 2), k," => ", v)
                    else
                        print(io, k, " => ", v)
                    end
                end
            else
                padded_print(field_name, value)
            end
        end
    end

    print("\n", " "^(indent÷2))
    print(io, length(t), "-element ", typeof(t.sta), ":")
    if length(t) > 0
        for field in propertynames(first(t.sta))
            field === :pos && continue
            values = getproperty(t.sta, field)
            if !all(ismissing, values)
                field_name = string("sta.", field)
                if field === :meta
                    padded_print(field_name, "Unique keys: ",
                        unique(Iterators.flatten(keys(meta) for meta in values))
                    )
                else
                    padded_print(field_name, values)
                end
            end
        end
    end

    for field in propertynames(t)
        if field in (:b, :delta, :sta, :evt, :picks, :meta, :data)
            continue
        else
            padded_print(field, getproperty(t, field))
        end
    end

    print(io, "\n Trace:")
    padded_print("picks", "")
    print(io, length(t.picks))
    padded_print("meta", "")
    Seis.show_dict(io, t.meta, hdr_string_len, indent)
end
