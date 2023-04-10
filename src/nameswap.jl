function swap_id_name(k, v; prop = :name)
    new_prop = getproperty(v, prop)
    @set! v.$prop = k
    new_prop, v
end

function replace_math_ident(id_name_dict, node)
    if typeof(node) == SBML.MathIdent
        return SBML.MathIdent(id_name_dict[node.id])
    elseif typeof(node) == SBML.MathApply
        return SBML.MathApply(node.fn, replace_math_ident(id_name_dict, node.args))
    elseif isa(node, AbstractArray)
        return map(x -> replace_math_ident(id_name_dict, x), node)
    else
        return node
    end
end

"SimBiology makes the names of everything a hashed ID, but keeps the ID in the name"
function fix_simbio_names(m::SBML.Model)
    cid_name_d = Dict(keys(m.compartments) .=> getproperty.(values(m.compartments), :name))
    pid_name_d = Dict(keys(m.parameters) .=> getproperty.(values(m.parameters), :name))
    rid_name_d = Dict(keys(m.reactions) .=> getproperty.(values(m.reactions), :name))

    for prop in [:compartments, :parameters, :reactions]
        xs = getproperty(m, prop)
        new_xs = Dict()
        for (k, v) in xs
            nk, nv = swap_id_name(k, v)
            new_xs[nk] = nv
        end
        @set! m.$prop = new_xs
    end

    new_ss = Dict()
    for (k, v) in m.species
        nk, nv = swap_id_name(k, v)
        cname = cid_name_d[v.compartment]
        # for species, we want to concat the name and compartment name like in simbiology and MTK
        sname = string(cname, "₊", nk)
        @set! nv.compartment = cname
        new_ss[sname] = nv
    end
    @set! m.species = new_ss

    sid_name_d = Dict(reverse.(keys(new_ss) .=> getproperty.(values(new_ss), :name)))

    ds = [cid_name_d, pid_name_d, rid_name_d, sid_name_d]
    slen = sum(length.(ds))
    id_name_dict = merge(ds...)
    @assert slen == length(id_name_dict)

    for (k, v) in m.reactions
        for (i, sr) in enumerate(v.reactants)
            @set! sr.species = id_name_dict[sr.species]
            m.reactions[k].reactants[i] = sr
        end
        for (i, sr) in enumerate(v.products)
            @set! sr.species = id_name_dict[sr.species]
            m.reactions[k].products[i] = sr
        end
        @set! m.reactions[k].kinetic_math = replace_math_ident(id_name_dict, v.kinetic_math)
    end

    for (i, r) in enumerate(m.rules)
        @set! m.rules[i].variable = id_name_dict[r.variable]
        new_tree = replace_math_ident(id_name_dict, r.math)
        @set! m.rules[i].math = new_tree
    end

    m
end
