# Conversion to symbolics
const IV = Catalyst.DEFAULT_IV
const D = Differential(IV)
symbolicsRateOf(x) = D(x)
const symbolics_mapping = Dict(SBML.default_function_mapping...,
                               "rateOf" => symbolicsRateOf)

# const SUBSDICT = get_substitutions(model)

map_symbolics_ident(x) = begin
    sym = Symbol(x.id)
    first(@species $sym)
end

function interpret_as_num(x::SBML.Math)
    SBML.interpret_math(x;
                        map_apply = (x::SBML.MathApply, interpret::Function) -> Num(symbolics_mapping[x.fn](interpret.(x.args)...)),
                        map_const = (x::SBML.MathConst) -> Num(SBML.default_constants[x.id]),
                        map_ident = map_symbolics_ident,
                        map_lambda = (_, _) -> throw(ErrorException("Symbolics.jl does not support lambda functions")),
                        map_time = (x::SBML.MathTime) -> IV,
                        map_value = (x::SBML.MathVal) -> x.val,
                        map_avogadro = (x::SBML.MathAvogadro) -> SBML.default_constants["avogadro"])
end

"""
Get dictonary to change types in kineticLaw
"""
function get_substitutions(model)
    u0map, parammap = get_mappings(model)
    subsdict = Dict()
    for item in first.(u0map)
        k = create_var(string(item.f.name))
        subsdict[k] = item
    end
    for item in first.(parammap)
        k = create_var(string(item.name))
        subsdict[k] = item
    end
    subsdict
end

function create_var(x; isbcspecies = false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@species $sym [isbcspecies = isbcspecies]))
end
function create_var(x, iv; isbcspecies = false, irreducible = false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@species $sym(iv) [
                               isbcspecies = isbcspecies,
                               irreducible = irreducible,
                           ]))
end
function create_param(x; isconstantspecies = false)
    sym = Symbol(x)
    Symbolics.unwrap(first(@parameters $sym [isconstantspecies = isconstantspecies]))
end

function has_rule_type(id::String, m::SBML.Model, T::Type{<:SBML.Rule})
    T == SBML.AlgebraicRule &&
        return any(SBML.isfreein(id, r.math) for r in m.rules if r isa SBML.AlgebraicRule)
    any(r.variable == id for r in m.rules if r isa T)
end
